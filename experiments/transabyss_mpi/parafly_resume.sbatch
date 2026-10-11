#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J parafly_resume
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem 250G
#SBATCH --output=parafly_resume_%j.log
# Step 1 of DEV26_PLAN.md: when a Trinity phase 2 is stopped part way and
# started again on more cores, does it pick up where it left off? dev26 rests
# on ParaFly skipping the commands listed in recursive_trinity.cmds.completed.
#
#     sbatch --export=ALL,R1=<cor_1.fq.gz>,R2=<cor_2.fq.gz> \
#         experiments/transabyss_mpi/parafly_resume.sbatch
#
# R1/R2: the corrected pair an ORP run assembled from,
#   <dir>/rcorr/<run>.TRIM_{1,2}P.cor.fq(.gz).
# Optional: NAME (default SRR1138704), OUT (default ./parafly_resume_<jobid>),
#   STRAND=RF|FR (as oyster.py --strand), NORMALIZE=0 (as oyster.py without
#   --normalize-reads; default 1, as the dev24 runs), SMALL_CPU (phase 2's
#   cores before the stop, default 8), KILL_AFTER (seconds of that phase 2
#   before the stop, default 900), GRACE (seconds from SIGTERM to SIGKILL,
#   default 30), PHASE1=<OUT>/phase1/<NAME>.trinity from an earlier run of
#   this script, to skip phase 1.
#
# What it does, all on this one node:
#   1. phase 1 as run_trinity_phase1 runs it (--no_distributed_trinity_exec)
#      into OUT/phase1/, kept untouched;
#   2. copies it to OUT/resume/ and OUT/straight/ (Trinity writes absolute
#      paths into its command list; the copies are pointed at themselves);
#   3. resume/: phase 2 at --CPU SMALL_CPU in its own session, stopped after
#      KILL_AFTER s (SIGTERM to the session, SIGKILL after GRACE s), then
#      relaunched at --CPU <all cores>;
#   4. straight/: phase 2 at --CPU <all cores>, never stopped;
#   5. prints a summary, also written to OUT/summary.txt.
# Phase 2 runs without --full_cleanup so both working dirs can be inspected.
# Takes ~3.5 h for SRR1138704 (phase 1 ~40 min, phase 2 ~1h20m twice).

set -euo pipefail
: "${R1:?set R1}" "${R2:?set R2}"
NAME="${NAME:-SRR1138704}"
OUT="${OUT:-$PWD/parafly_resume_${SLURM_JOB_ID:-$$}}"
STRAND="${STRAND:-}"
NORMALIZE="${NORMALIZE:-1}"
SMALL_CPU="${SMALL_CPU:-8}"
KILL_AFTER="${KILL_AFTER:-900}"
GRACE="${GRACE:-30}"
PHASE1="${PHASE1:-}"
CPU="${SLURM_CPUS_PER_TASK:-40}"
MEM=$(( ${SLURM_MEM_PER_NODE:-256000} / 1024 ))
MAXMEM=$(( MEM * 9 / 10 ))

for f in "$R1" "$R2"; do [ -f "$f" ] || { echo "missing read file: $f" >&2; exit 1; }; done
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"
SUMMARY="$OUT/summary.txt"

module load anaconda/colsa

now() { date +%s; }
hms() { local s=$1; printf '%dh%02dm%02ds' $((s / 3600)) $((s % 3600 / 60)) $((s % 60)); }
log() { echo "[$(date '+%F %T')] $*"; }

# The command _trinity_base_cmd builds, for an output dir and a core count.
trinity_cmd() {
    local outdir=$1 cpu=$2
    shift 2
    TCMD=(conda run --no-capture-output -n orp_trinity Trinity)
    case "$STRAND" in RF|FR) TCMD+=(--SS_lib_type "$STRAND") ;; esac
    TCMD+=(--no_version_check --bypass_java_version_check)
    [ "$NORMALIZE" = 1 ] || TCMD+=(--no_normalize_reads)
    TCMD+=(--seqType fq --output "$outdir" --max_memory "${MAXMEM}G"
           --left "$R1" --right "$R2" --CPU "$cpu" --inchworm_cpu 10 "$@")
}

# -- sessions ----------------------------------------------------------------
# Each phase 2 runs as `setsid <cmd> &`: bash without job control leaves a
# background child outside any group it leads, so setsid makes it a session
# leader without forking, and $! is the session id -- the same arrangement
# as Popen(start_new_session=True) in the dev26 plan.

SESSIONS=()
cleanup() {
    local s
    for s in ${SESSIONS[@]+"${SESSIONS[@]}"}; do pkill -KILL -s "$s" 2>/dev/null || true; done
}
trap cleanup EXIT
trap 'exit 143' TERM INT

start_session() {   # logfile cmd... -> SID
    local logf=$1
    shift
    setsid "$@" > "$logf" 2>&1 < /dev/null &
    SID=$!
    SESSIONS+=("$SID")
    sleep 1
    local sid
    sid=$(ps -o sid= -p "$SID" 2>/dev/null | tr -d ' ')
    if [ "$sid" != "$SID" ]; then
        echo "pid $SID is not a session leader (sid '$sid'); log:" >&2
        tail -20 "$logf" >&2
        exit 1
    fi
}

# Live (non-zombie) processes in a session: the leader stays a zombie until
# it is waited for, and pgrep would count it.
alive_in() { ps -o pid=,stat=,comm=,args= -s "$1" 2>/dev/null | awk '$2 !~ /^Z/' || true; }
parafly_up() { ps -o stat=,comm= -s "$1" 2>/dev/null | awk '$1 !~ /^Z/ && $2 == "ParaFly" {f = 1} END {exit !f}'; }

stop_session() {    # sid -> STOP_HOW, STOP_SECS
    local sid=$1 t0
    t0=$(now)
    pkill -TERM -s "$sid" 2>/dev/null || true
    STOP_HOW="SIGTERM"
    while [ -n "$(alive_in "$sid")" ]; do
        if (( $(now) - t0 >= GRACE )); then
            STOP_HOW="SIGKILL after ${GRACE}s (SIGTERM left: $(alive_in "$sid" | awk '{print $3}' | sort | uniq -c | xargs))"
            pkill -KILL -s "$sid" 2>/dev/null || true
            sleep 2
            break
        fi
        sleep 1
    done
    STOP_SECS=$(( $(now) - t0 ))
    wait "$sid" 2>/dev/null || true
}

# Trinity/ParaFly/Butterfly/Chrysalis-family processes still running anywhere
# in this job (scontrol listpids), or, if that is unavailable, as this user
# on this node (which may include other jobs').
survivors() {
    local pids
    pids=$(scontrol listpids "${SLURM_JOB_ID:-0}" 2>/dev/null | awk 'NR > 1 && $1 ~ /^[0-9]+$/ {print $1}' | paste -sd, -)
    if [ -n "$pids" ]; then
        ps -o pid=,sid=,stat=,comm=,args= -p "$pids" 2>/dev/null
    else
        ps -u "$USER" -o pid=,sid=,stat=,comm=,args=
    fi | awk '$3 !~ /^Z/' \
       | grep -E 'Trinity|ParaFly|java|Butterfly|Chrysalis|ReadsToTranscripts|GraphFromFasta|QuantifyGraph|FastaToDeBruijn|inchworm|jellyfish|bowtie2|salmon' \
       | grep -vE 'grep|ps -o' || true
}

# Watch a session until it ends (or until `limit` s after t0): when ParaFly
# starts and stops, and when the .completed file first grows past `base`.
# Sets W_PF_START, W_PF_END, W_FIRST (s after t0; empty if not seen),
# W_DONE (1 if the session ended) and W_RC.
watch_run() {       # sid t0 completed_file base [limit]
    local sid=$1 t0=$2 cf=$3 base=$4 limit=${5:-0} t
    W_PF_START=""; W_PF_END=""; W_FIRST=""; W_DONE=0; W_RC=""
    while [ -n "$(alive_in "$sid")" ]; do
        t=$(( $(now) - t0 ))
        if parafly_up "$sid"; then
            [ -n "$W_PF_START" ] || W_PF_START=$t
        elif [ -n "$W_PF_START" ] && [ -z "$W_PF_END" ]; then
            W_PF_END=$t
        fi
        if [ -z "$W_FIRST" ] && [ -f "$cf" ] && (( $(wc -l < "$cf") > base )); then
            W_FIRST=$t
        fi
        if (( limit > 0 && t >= limit )); then return 0; fi
        if (( t < 180 )); then sleep 1; else sleep 5; fi
    done
    t=$(( $(now) - t0 ))
    [ -z "$W_PF_START" ] || [ -n "$W_PF_END" ] || W_PF_END=$t
    W_DONE=1
    W_RC=0
    wait "$sid" || W_RC=$?
}

lines() { if [ -f "$1" ]; then wc -l < "$1" | tr -d ' '; else echo 0; fi; }
out_dirs() { { find "$1/read_partitions" -type d -name '*.out' 2>/dev/null || true; } | wc -l | tr -d ' '; }
# finished components' outputs (each component's Trinity runs --full_cleanup,
# leaving <reads>.out.Trinity.fasta), with inode and mtime
component_outputs() { { find "$1/read_partitions" -type f -name '*.Trinity.fasta' -printf '%p\t%i\t%T@\n' 2>/dev/null || true; } | sort; }

fasta_stats() {     # fasta -> "transcripts genes bases"
    if [ ! -s "$1" ]; then echo "- - -"; return; fi
    awk '/^>/ {n++; g = $1; sub(/^>/, "", g); sub(/_i[0-9]+$/, "", g); genes[g] = 1; next}
         {b += length($0)}
         END {printf "%d %d %d\n", n, length(genes), b}' "$1"
}

copy_phase1() {     # src dst -> COPY_FIXED
    local src=$1 dst=$2 f
    rm -rf "$dst"
    mkdir -p "$(dirname "$dst")"
    cp -a "$src" "$dst"
    # Trinity writes absolute paths (the per-component read files, inside
    # its output dir) into recursive_trinity.cmds and
    # partitioned_reads.files.list. Left alone, the copy's phase 2 would
    # assemble into the original. Checkpoints are named files, not
    # command strings, so nothing else needs to change.
    COPY_FIXED=0
    while IFS= read -r -d '' f; do
        if grep -qF "$src" "$f"; then
            sed -i "s#$src#$dst#g" "$f"
            COPY_FIXED=$(( COPY_FIXED + 1 ))
        fi
    done < <(find "$dst" -path "$dst/read_partitions" -prune -o -type f -size -500M -print0)
    if find "$dst" -path "$dst/read_partitions" -prune -o -type f -size -500M -print0 \
            | xargs -0 grep -lF "$src" 2>/dev/null; then
        echo "copy $dst still refers to $src (files above)" >&2
        exit 1
    fi
}

# -- environment -------------------------------------------------------------

log "=== ParaFly resume test: $NAME on $(hostname), ${CPU} cpus, ${MEM}G (Trinity --max_memory ${MAXMEM}G) ==="
log "R1=$R1"
log "R2=$R2"
log "OUT=$OUT  SMALL_CPU=$SMALL_CPU  KILL_AFTER=${KILL_AFTER}s  GRACE=${GRACE}s  NORMALIZE=$NORMALIZE  STRAND=${STRAND:-none}"

TRINITY_BIN=$(conda run -n orp_trinity bash -c 'readlink -f "$(command -v Trinity)"')
TRINITY_HOME=$(dirname "$TRINITY_BIN")
PARAFLY="$TRINITY_HOME/trinity-plugins/BIN/ParaFly"
log "Trinity: $TRINITY_BIN ($(git -C "$TRINITY_HOME" describe --always --dirty 2>/dev/null || echo 'not a git checkout'))"
[ -x "$PARAFLY" ] || { echo "no ParaFly at $PARAFLY" >&2; exit 1; }
gomp=$(nm -D "$PARAFLY" 2>/dev/null | grep -c GOMP_ || true)
pfcheck="$OUT/parafly_check"
mkdir -p "$pfcheck"
printf 'sleep 2\nsleep 2\n' > "$pfcheck/cmds"
rm -f "$pfcheck/cmds.completed"
cat > "$pfcheck/time.sh" <<'EOF2'
s=$(date +%s%N); "$1" -c "$2" -CPU 2 > /dev/null 2>&1; e=$(date +%s%N)
echo $(( (e - s) / 1000000 ))
EOF2
pf_ms=$(conda run -n orp_trinity bash "$pfcheck/time.sh" "$PARAFLY" "$pfcheck/cmds" | tail -1)
log "ParaFly: $PARAFLY; GOMP symbols: ${gomp:-?}; two 'sleep 2' at -CPU 2: ${pf_ms} ms"
if ! [ "${pf_ms:-x}" -eq "${pf_ms:-x}" ] 2>/dev/null || [ "$pf_ms" -gt 3600 ]; then
    echo "ParaFly runs its commands one at a time (no OpenMP): every timing here would be wrong" >&2
    exit 1
fi

# -- phase 1 -----------------------------------------------------------------

P1="$OUT/phase1/$NAME.trinity"
if [ -n "$PHASE1" ]; then
    [ -f "$PHASE1/recursive_trinity.cmds.ok" ] || { echo "$PHASE1: no recursive_trinity.cmds.ok, not a finished phase 1" >&2; exit 1; }
    [ ! -e "$PHASE1/recursive_trinity.cmds.completed" ] || { echo "$PHASE1 has run phase 2; give a phase-1-only dir" >&2; exit 1; }
    P1="$PHASE1"
    P1_SECS=0
    log "phase 1: reusing $P1"
else
    mkdir -p "$OUT/phase1"
    trinity_cmd "$P1" "$CPU" --no_distributed_trinity_exec
    log "phase 1: ${TCMD[*]}"
    t0=$(now)
    "${TCMD[@]}" > "$OUT/phase1/phase1.log" 2>&1 || { echo "phase 1 failed; see $OUT/phase1/phase1.log" >&2; exit 1; }
    P1_SECS=$(( $(now) - t0 ))
    [ -f "$P1/recursive_trinity.cmds.ok" ] || { echo "phase 1 left no recursive_trinity.cmds.ok" >&2; exit 1; }
    log "phase 1 done in $(hms "$P1_SECS")"
fi
TOTAL=$(lines "$P1/recursive_trinity.cmds")
log "recursive_trinity.cmds: $TOTAL commands"

A="$OUT/resume/$NAME.trinity"
B="$OUT/straight/$NAME.trinity"
t0=$(now)
copy_phase1 "$P1" "$A"; FIXED_A=$COPY_FIXED
copy_phase1 "$P1" "$B"; FIXED_B=$COPY_FIXED
COPY_SECS=$(( $(now) - t0 ))
touch "$OUT/copies.done"
log "copied phase 1 to resume/ and straight/ in $(hms "$COPY_SECS") (paths rewritten in $FIXED_A and $FIXED_B files)"
cd "$OUT"

# -- resume/: phase 2 at SMALL_CPU, stopped --------------------------------------

CF="$A/recursive_trinity.cmds.completed"
trinity_cmd "$A" "$SMALL_CPU"
log "phase 2 at --CPU $SMALL_CPU: ${TCMD[*]}"
t_small=$(now)
start_session "$OUT/resume/phase2_small.log" "${TCMD[@]}"
SID_SMALL=$SID
watch_run "$SID_SMALL" "$t_small" "$CF" 0 "$KILL_AFTER"
SMALL_PF_START=$W_PF_START
if [ "$W_DONE" = 1 ]; then
    echo "phase 2 at --CPU $SMALL_CPU ended (exit $W_RC) before the stop at ${KILL_AFTER}s: nothing to test" >&2
    tail -30 "$OUT/resume/phase2_small.log" >&2
    exit 1
fi
[ -n "$SMALL_PF_START" ] || log "WARNING: ParaFly had not started by the stop"
SMALL_SECS=$(( $(now) - t_small ))
DONE_BEFORE=$(lines "$CF")
INFLIGHT_BEFORE=$(out_dirs "$A")
log "stopping session $SID_SMALL after $(hms "$SMALL_SECS"): $DONE_BEFORE/$TOTAL completed, $INFLIGHT_BEFORE component dirs in flight"
alive_in "$SID_SMALL" | awk '{print $3}' | sort | uniq -c | sed 's/^/    running: /'
stop_session "$SID_SMALL"
DONE_KILL=$(lines "$CF")
INFLIGHT_KILL=$(out_dirs "$A")
SURV_SESSION=$(alive_in "$SID_SMALL")
sleep 5
SURV_JOB=$(survivors)
component_outputs "$A" > "$OUT/resume/outputs_at_kill.tsv"
OUTPUTS_KILL=$(lines "$OUT/resume/outputs_at_kill.tsv")
log "stopped by $STOP_HOW in ${STOP_SECS}s: $DONE_KILL/$TOTAL completed, $INFLIGHT_KILL partial component dirs, $OUTPUTS_KILL component outputs"
if [ -n "$SURV_SESSION$SURV_JOB" ]; then
    log "SURVIVORS after the stop:"
    printf '%s\n%s\n' "$SURV_SESSION" "$SURV_JOB" | sed '/^$/d; s/^/    /'
else
    log "no Trinity/ParaFly/java/Chrysalis processes left"
fi
cp "$CF" "$OUT/resume/completed_at_kill.txt" 2>/dev/null || : > "$OUT/resume/completed_at_kill.txt"

# -- resume/: relaunch on all cores --------------------------------------------

trinity_cmd "$A" "$CPU"
log "relaunch at --CPU $CPU: ${TCMD[*]}"
t_rl=$(now)
start_session "$OUT/resume/phase2_relaunch.log" "${TCMD[@]}"
watch_run "$SID" "$t_rl" "$CF" "$DONE_KILL"
RL_RC=$W_RC; RL_PF_START=$W_PF_START; RL_PF_END=$W_PF_END; RL_FIRST=$W_FIRST
RL_SECS=$(( $(now) - t_rl ))
log "relaunch exited $RL_RC after $(hms "$RL_SECS")"

# -- straight/: uninterrupted phase 2 on all cores ------------------------------

trinity_cmd "$B" "$CPU"
log "uninterrupted phase 2 at --CPU $CPU: ${TCMD[*]}"
t_st=$(now)
start_session "$OUT/straight/phase2.log" "${TCMD[@]}"
watch_run "$SID" "$t_st" "$B/recursive_trinity.cmds.completed" 0
ST_RC=$W_RC; ST_PF_START=$W_PF_START; ST_PF_END=$W_PF_END
ST_SECS=$(( $(now) - t_st ))
log "uninterrupted phase 2 exited $ST_RC after $(hms "$ST_SECS")"

# -- checks --------------------------------------------------------------------

RL_LOG="$OUT/resume/phase2_relaunch.log"
RL_NCMDS=$({ grep -ao 'Number of Commands: [0-9]*' "$RL_LOG" || true; } | tail -1 | awk '{print $4}')
RL_SKIPMSG=$(grep -ac 'has successfully completed from a previous run' "$RL_LOG" || true)
RL_SKIPMSG=${RL_SKIPMSG:-0}
EXPECT_NCMDS=$(( TOTAL - DONE_KILL ))
DONE_FINAL=$(lines "$CF")
DUPS=$({ sort "$CF" 2>/dev/null || true; } | uniq -d | wc -l | tr -d ' ')
PREFIX_OK=no
if [ "$DONE_KILL" = 0 ] || cmp -s "$OUT/resume/completed_at_kill.txt" <(head -n "$DONE_KILL" "$CF"); then PREFIX_OK=yes; fi
component_outputs "$A" > "$OUT/resume/outputs_final.tsv"
# a finished component that ran again has a new inode or mtime, or is gone
REDONE=$(awk -F'\t' 'FILENAME == ARGV[1] {after[$1] = $2 "\t" $3; next}
                     !($1 in after) || after[$1] != $2 "\t" $3 {n++}
                     END {print n + 0}' "$OUT/resume/outputs_final.tsv" "$OUT/resume/outputs_at_kill.tsv")
LEFT_A=$(out_dirs "$A")
LEFT_B=$(out_dirs "$B")
B_DONE=$(lines "$B/recursive_trinity.cmds.completed")
FA_A="$A.Trinity.fasta"
FA_B="$B.Trinity.fasta"
read -r TX_A GN_A BP_A <<< "$(fasta_stats "$FA_A")"
read -r TX_B GN_B BP_B <<< "$(fasta_stats "$FA_B")"
OUTN_A=$(lines "$OUT/resume/outputs_final.tsv")
OUTN_B=$(component_outputs "$B" | wc -l | tr -d ' ')
P1_LEAK=$({ find "$P1" -newer "$OUT/copies.done" 2>/dev/null || true; } | wc -l | tr -d ' ')

pass() { if [ "$1" = 1 ]; then echo PASS; else echo FAIL; fi; }
pct() { awk -v a="$1" -v b="$2" 'BEGIN {if (b > 0) printf "%.1f%%", 100 * a / b; else print "-"}'; }
diffpct() { awk -v a="$1" -v b="$2" 'BEGIN {if (b > 0 && a != "-") printf "%+d (%+.2f%%)", a - b, 100 * (a - b) / b; else print "-"}'; }
secs() { if [ -n "$1" ]; then hms "$1"; else echo "not seen"; fi; }

c_clean=$([ -z "$SURV_SESSION$SURV_JOB" ] && echo 1 || echo 0)
c_skip=$([ "${RL_NCMDS:-x}" = "$EXPECT_NCMDS" ] && { [ "$DONE_KILL" = 0 ] || [ "$RL_SKIPMSG" -gt 0 ]; } && echo 1 || echo 0)
c_once=$([ "$DUPS" = 0 ] && [ "$REDONE" = 0 ] && [ "$PREFIX_OK" = yes ] && echo 1 || echo 0)
c_finish=$([ "$RL_RC" = 0 ] && [ -s "$FA_A" ] && [ "$DONE_FINAL" = "$TOTAL" ] && echo 1 || echo 0)
c_partial=$([ "$LEFT_A" = "$LEFT_B" ] && echo 1 || echo 0)
c_straight=$([ "$ST_RC" = 0 ] && [ -s "$FA_B" ] && echo 1 || echo 0)

# ParaFly's own time on the remaining commands, against the uninterrupted
# ParaFly's scaled to the same share of the commands
REMAIN_FRAC=$(awk -v t="$TOTAL" -v d="$DONE_KILL" 'BEGIN {printf "%.4f", (t - d) / t}')
RL_PF=""; ST_PF=""; EXPECT_PF=""
[ -z "$RL_PF_START" ] || [ -z "$RL_PF_END" ] || RL_PF=$(( RL_PF_END - RL_PF_START ))
[ -z "$ST_PF_START" ] || [ -z "$ST_PF_END" ] || ST_PF=$(( ST_PF_END - ST_PF_START ))
[ -z "$ST_PF" ] || EXPECT_PF=$(awk -v s="$ST_PF" -v f="$REMAIN_FRAC" 'BEGIN {printf "%d", s * f}')

{
echo "=================== ParaFly resume test: $NAME ==================="
echo "node $(hostname), $CPU cpus; Trinity $TRINITY_BIN; ParaFly GOMP symbols ${gomp:-?}"
echo "phase 1: $(if [ "$P1_SECS" = 0 ]; then echo "reused $P1"; else hms "$P1_SECS"; fi); copies $(hms "$COPY_SECS")"
echo
echo "-- at the stop (phase 2 at --CPU $SMALL_CPU, stopped after $(hms "$SMALL_SECS"); ParaFly started at $(secs "$SMALL_PF_START"))"
echo "commands total                   $TOTAL"
echo "completed when signalled         $DONE_BEFORE ($(pct "$DONE_BEFORE" "$TOTAL"))"
echo "completed once all had exited    $DONE_KILL ($(pct "$DONE_KILL" "$TOTAL"))"
echo "partial component dirs           $INFLIGHT_KILL (rerun from scratch)"
echo "stopped by                       $STOP_HOW, all gone ${STOP_SECS}s after SIGTERM"
echo "survivors                        $(if [ "$c_clean" = 1 ]; then echo none; else echo "SEE LOG"; fi)"
echo
echo "-- relaunch at --CPU $CPU"
echo "ParaFly 'Number of Commands'     ${RL_NCMDS:-not printed} (expected $TOTAL - $DONE_KILL = $EXPECT_NCMDS)"
echo "ParaFly skip warnings            $RL_SKIPMSG (it prints at most 9)"
echo ".completed at the end            $DONE_FINAL of $TOTAL, duplicate lines $DUPS, starts with the pre-stop list: $PREFIX_OK"
echo "pre-stop outputs rewritten       $REDONE of $OUTPUTS_KILL"
echo "partial dirs left at the end     $LEFT_A (uninterrupted run: $LEFT_B)"
echo "exit                             $RL_RC; $(if [ -s "$FA_A" ]; then echo "wrote $(basename "$FA_A")"; else echo "NO $(basename "$FA_A")"; fi)"
echo
echo "-- restart overhead"
echo "relaunch -> ParaFly running      $(secs "$RL_PF_START") (uninterrupted run: launch -> ParaFly $(secs "$ST_PF_START"))"
echo "relaunch -> first new completion $(secs "$RL_FIRST")"
echo "ParaFly on the remaining $(pct "$((TOTAL - DONE_KILL))" "$TOTAL")   $(secs "$RL_PF") (uninterrupted ParaFly $(secs "$ST_PF") x $REMAIN_FRAC = $(secs "$EXPECT_PF"))"
echo "relaunch wall                    $(hms "$RL_SECS"); uninterrupted wall $(hms "$ST_SECS"); stopped run $(hms "$SMALL_SECS") at --CPU $SMALL_CPU"
echo
echo "-- assemblies (transcripts / genes / bases)"
echo "resumed                          $TX_A / $GN_A / $BP_A  ($OUTN_A component outputs)"
echo "uninterrupted                    $TX_B / $GN_B / $BP_B  ($OUTN_B component outputs, $B_DONE/$TOTAL completed, exit $ST_RC)"
echo "resumed - uninterrupted          transcripts $(diffpct "$TX_A" "$TX_B"), genes $(diffpct "$GN_A" "$GN_B")"
echo "(byte identity is not expected: Trinity output depends on run dir and JVM hash, NOTES.md 2026-10-05)"
echo
echo "-- checks"
echo "$(pass "$c_clean")  stop leaves no Trinity/ParaFly/java/Chrysalis processes"
echo "$(pass "$c_skip")  relaunch skips the completed commands"
echo "$(pass "$c_once")  no command ran twice (no duplicate .completed lines, pre-stop outputs untouched)"
echo "$(pass "$c_finish")  relaunch finishes: exit 0, Trinity.fasta, every command completed"
echo "$(pass "$c_partial")  partial component dirs don't survive the relaunch"
echo "$(pass "$c_straight")  uninterrupted run finishes"
[ "$P1_LEAK" = 0 ] || echo "FAIL  $P1_LEAK files in the pristine phase-1 dir changed after the copies: a copy ran into it"
echo
echo "fastas: $FA_A"
echo "        $FA_B"
echo "logs:   $OUT/resume/phase2_small.log, $RL_LOG, $OUT/straight/phase2.log"
} | tee "$SUMMARY"
