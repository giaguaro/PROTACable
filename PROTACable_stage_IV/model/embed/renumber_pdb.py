"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""

import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file) as infile, open(output_file, "w") as outfile:
    c = 0
    ter_found = False
    for line in infile:
        if line.startswith("ATOM"):
            c = int(line[6:11])
        if line.startswith("TER"):
            ter_found = True
            c += 1
            f = int(line[6:11])
            if f == c:
                c += 1
        if line.startswith("HETATM"):
            if not ter_found: c += 1
            line = line[:6] + f"{c:5d}" + line[11:]
            outfile.write(line)
            c += 1
            continue
        outfile.write(line)
