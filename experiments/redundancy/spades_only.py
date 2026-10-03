#!/usr/bin/env python3
"""Run only oyster.py's two SPAdes assemblies (spadesauto, spadeshigh).

Takes oyster.py's own command line (--read1/--read2/--trimmed-corrected-reads/
--spades1-kmer/--runout/--dir/--cpu/--mem ...), so the k values resolve exactly
as in a full run and print as the "[kmer]" lines. Writes
<dir>/assemblies/<runout>.spades{auto,high}.fasta.

usage: spades_only.py --oyster CODE_DIR/oyster.py <oyster.py arguments>
"""
import importlib.util
import sys
from pathlib import Path

i = sys.argv.index("--oyster")
oyster_path = Path(sys.argv[i + 1]).resolve()
del sys.argv[i:i + 2]
sys.path.insert(0, str(oyster_path.parent))
spec = importlib.util.spec_from_file_location("oyster", oyster_path)
oyster = importlib.util.module_from_spec(spec)
spec.loader.exec_module(oyster)

p = oyster.Pipeline(oyster.parse_args())
p.setup()
p.readcheck()
p.use_corrected_reads()
p.run_spadesauto()
p.run_spadeshigh()
print("[spades_only] done")
