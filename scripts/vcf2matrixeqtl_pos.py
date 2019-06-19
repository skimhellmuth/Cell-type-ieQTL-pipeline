#!/usr/bin/env python3

import sys

def main():
    """Extract SNP ID, chr, pos from the VCF file"""
    for line in sys.stdin:
        xcols = line.rstrip().split("\t")
        if line.startswith("#CHR"):
            xout = ['snp','chr_snp','pos']
            print("\t".join(xout))
        elif line.startswith("#") == False:
            xout = [xcols[2],xcols[0],xcols[1]]
       	    print("\t".join(xout))


if __name__=='__main__':
    main()
