#!/bin/bash
# build_bfly.sh OUTJAR : compile Butterfly sources in the submodule, overlay onto the shipped fat jar
set -euo pipefail
E=~/orp_envs/orp_trinity; export PATH=$E/bin:$PATH
B=~/trinity_eng/src/Butterfly/Butterfly/src
OUT=$1
W=$(mktemp -d /tmp/bfbuild.XXXX)
mkdir -p $W/fat $W/classes
(cd $W/fat && unzip -q ~/trinity_eng/src/Butterfly/Butterfly.jar)
CP=$(ls $B/lib/*.jar | tr "\n" ":")
find $B/src -name "*.java" > $W/sources.txt
javac -nowarn --release 8 -encoding UTF-8 -cp "$CP" -d $W/classes @$W/sources.txt 2>&1 | grep -v -E "^Note:|warning" | head -20 || true
cp -r $W/classes/* $W/fat/
(cd $W/fat && jar --create --file $OUT --manifest META-INF/MANIFEST.MF -C . . 2>/dev/null || jar cfm $OUT META-INF/MANIFEST.MF -C . .)
ls -la $OUT; rm -rf $W
