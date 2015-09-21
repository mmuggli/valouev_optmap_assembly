#!/usr/bin/python3

import sys
import re

if len(sys.argv) < 4:
    print("Usage: ", sys.argv[0], " <maps> <ovlps> <patched_ovlps>")
    sys.exit(1)

names = {}
for lno,line in enumerate(open(sys.argv[1])):
    if lno % 3 == 0: names[line.strip()] = int(lno / 3)

def lookup(name):
    return str(names[name])
    
patched = open(sys.argv[3], "w")
for lno,line in enumerate(open(sys.argv[2])):
    if lno % 3 == 0:
        fields = line.split()
        patched.write(lookup(fields[0]) + " " + lookup(fields[1]) + " " + line)
    else:
        patched.write(line)
patched.close()
