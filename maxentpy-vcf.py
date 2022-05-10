#! /usr/bin/env python3

import argparse
from maxentpy_vcf import run_maxentpy


def main(args):
    run_maxentpy(args.vcf, args.genome, args.refgene, args.output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Batch run AutoPVS1')
    parser.add_argument('--vcf', '-i', required=True, help='VCF file')
    parser.add_argument('--output', '-o', required=True, help='output file')
    parser.add_argument('--genome', '-g', required=True, help='reference fasta file')
    parser.add_argument('--refgene', '-r', required=True, help='refGene')
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
