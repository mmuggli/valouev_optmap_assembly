#!/usr/bin/python3

import sys
import re

if len(sys.argv) < 4:
    print("Usage: ", sys.argv[0], " <maps> <ovlps> <skipped>")
    sys.exit(1)

map_names = set()
for lno,line in enumerate(open(sys.argv[1])):
    if lno % 3 == 0: map_names.add(line.strip())

ovlp_map_names = set()
for lno,line in enumerate(open(sys.argv[2])):
    if lno % 3 == 0:
        fields = line.split()
        ovlp_map_names.add(fields[0])
        ovlp_map_names.add(fields[1])

skipped = open(sys.argv[3],"w")
for unaligned_name in map_names - ovlp_map_names:
    skipped.write(unaligned_name + "\n")
skipped.close()

