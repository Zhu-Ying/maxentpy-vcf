import csv
import vcf as pyvcf
from pyfaidx import Fasta
from .transcripts import read_refgene, make_transcript
from .splicing import SplicingMaxEnt
from .utils import vcf_to_av


def run_maxentpy(vcf_file: str, genome_file: str, refgene_file: str, outfile: str):
    genome = Fasta(genome_file)
    refgenes = read_refgene(refgene_file)
    vcf_reader = pyvcf.Reader(filename=vcf_file)
    with open(outfile, 'w') as fo:
        writer = csv.DictWriter(fo, fieldnames=[
            'Chr', 'Start', 'End', 'Ref', 'Alt',
            'Maxent_type', 'Maxent_pred',
            'Maxent_score_ref', 'Maxent_score_alt',
            'Maxent_score_var', 'Maxent_foldchange'
        ], delimiter='\t')
        writer.writeheader()
        for record in vcf_reader:
            vcf_chrom = str(record.CHROM)
            vcf_pos = int(record.POS)
            vcf_ref = str(record.REF)
            for vcf_alt in record.ALT:
                vcf_alt = str(vcf_alt)
                chrom, start, end, ref, alt = vcf_to_av(vcf_chrom, vcf_pos, vcf_ref, vcf_alt)
                vcf_chrom = 'chrM' if chrom == 'MT' else f'chr{chrom}'
                result = None
                for refgene in refgenes.query(f'Chrom == "{vcf_chrom}" and TxStart <= {end} and TxEnd >= {start}').iloc:
                    transcript = make_transcript(refgene)
                    splicing = SplicingMaxEnt(vcf_chrom, vcf_pos, vcf_ref, vcf_alt, transcript, genome)
                    if result:
                        if abs(splicing.maxentscore_var) > abs(result.maxentscore_var):
                            result = splicing
                        elif abs(splicing.maxentscore_var) == abs(result.maxentscore_var):
                            if splicing.maxentscore_alt > result.maxentscore_alt:
                                result = splicing
                    else:
                        result = splicing
                if result:
                    writer.writerow({
                        'Chr': chrom, 'Start': start, 'End': end, 'Ref': ref, 'Alt': alt,
                        'Maxent_type': result.splice_type,
                        'Maxent_pred': result.maxentpred,
                        'Maxent_score_ref': result.maxentscore_ref,
                        'Maxent_score_alt': result.maxentscore_alt,
                        'Maxent_score_var': result.maxentscore_var,
                        'Maxent_foldchange': result.maxentscore_foldchange
                    })
