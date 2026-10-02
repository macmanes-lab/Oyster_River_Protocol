#!/usr/bin/env bash
# Small-file and streaming I/O benchmark for one directory.
#
# usage: io_test.sh DIR [LABEL]
#
# Creates a scratch subdirectory under DIR, then times:
#   create 2k    2,000 small files in a fresh directory
#   create 20k   20,000 small files in one directory (OrthoFinder writes one
#                file per orthogroup, 50k-300k of them, into a single dir)
#   list 20k     ls of that directory
#   delete 20k   rm -rf of it
#   write 1G     1 GiB sequential write, fsync'd
#   read 1G      the same file read back with O_DIRECT, so not from page cache
# Each step has a time limit, and a step that runs out is reported as
# ">LIMIT s" with how far it got. The scratch directory is removed at the end.
set -u

DIR=${1:?usage: io_test.sh DIR [LABEL]}
LABEL=${2:-$DIR}
LIMIT=${LIMIT:-60}
T=$(mktemp -d -p "$DIR" .iotest.XXXX) || { echo "$LABEL: cannot create scratch dir"; exit 1; }
trap 'rm -rf "$T"' EXIT

now() { date +%s.%N; }
secs() { printf '%.2f' "$(echo "$2 - $1" | bc)"; }

create() {  # create N files in $T/$2, within LIMIT seconds
    local n=$1 d=$T/$2 i=0 s e
    mkdir -p "$d"; s=$(now)
    while (( i < n )); do
        echo x > "$d/f$i" || break
        i=$((i + 1))
        if (( i % 500 == 0 )) && (( $(echo "$(now) - $s > $LIMIT" | bc) )); then
            echo ">$LIMIT s ($i of $n)"; return
        fi
    done
    e=$(now); echo "$(secs "$s" "$e") s"
}

r_c2=$(create 2000 c2)
r_c20=$(create 20000 c20)
s=$(now); timeout "$LIMIT" ls -f "$T/c20" > /dev/null && r_ls="$(secs "$s" "$(now)") s" || r_ls=">$LIMIT s"
s=$(now); timeout "$LIMIT" rm -rf "$T/c20" && r_rm="$(secs "$s" "$(now)") s" || r_rm=">$LIMIT s"
s=$(now)
if timeout "$LIMIT" dd if=/dev/zero of="$T/big" bs=4M count=256 conv=fsync status=none; then
    e=$(now); r_w="$(printf '%.0f' "$(echo "1024 / ($e - $s)" | bc -l)") MB/s"
else r_w=">$LIMIT s"; fi
s=$(now)
if timeout "$LIMIT" dd if="$T/big" of=/dev/null bs=4M iflag=direct status=none; then
    e=$(now); r_r="$(printf '%.0f' "$(echo "1024 / ($e - $s)" | bc -l)") MB/s"
else r_r=">$LIMIT s"; fi

printf '%-34s create2k %-18s create20k %-20s ls20k %-9s rm20k %-9s write %-9s read %s\n' \
    "$LABEL" "$r_c2" "$r_c20" "$r_ls" "$r_rm" "$r_w" "$r_r"
