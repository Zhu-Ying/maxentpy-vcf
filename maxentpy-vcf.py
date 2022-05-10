#! /usr/bin/env python3

import argparse
from maxentpy_annovar import run_maxentpy


def main(args):
    gene_details = args.gene_detail or ["GeneDetail.refGeneWithVer", "AAChange.refGeneWithVer"]
    run_maxentpy(args.annovar, args.genome, args.refgene, gene_details, args.output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Batch run AutoPVS1')
    parser.add_argument('--annovar', '-i', required=True, help='input ANNOVAR result for SNV')
    parser.add_argument('--output', '-o', required=True, help='output file')
    parser.add_argument('--genome', '-g', required=True, help='reference fasta file')
    parser.add_argument('--refgene', '-r', required=True, help='refGene')
    parser.add_argument('--gene_detail', '-c', action='append', help='refGene')
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
