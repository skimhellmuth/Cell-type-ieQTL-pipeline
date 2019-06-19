#!/usr/bin/env python3

import sys

def main():
    """Extract SNP ID, and sum the alleles"""
    for line in sys.stdin:
        xcols = line.rstrip().split("\t")
        if line.startswith("#CHR"):
            xout = ['ID'] + xcols[9:len(xcols)]
            print("\t".join(xout))
        elif line.startswith("#") == False:
            gt_index = xcols[8].split(":").index("GT") # there can be additional genotype files separated by :, i.e. GQ - genotype quality,
            xout = [xcols[2]]
            for i in range(9,len(xcols)):
                xgt = xcols[i].split(":")[gt_index].split("|")
                if "." not in xgt:
                    xout.append(str(sum(map(int,xgt))))
                else:
                    xout.append("NA")
            print("\t".join(xout))


if __name__=='__main__':
    main();
