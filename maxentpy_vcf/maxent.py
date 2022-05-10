import re
import csv
import vcf as pyvcf
import pandas
from pyfaidx import Fasta
from .transcripts import read_refgene, make_transcript
from .splicing import SplicingMaxEnt

def vcf_to_av(chrom: str, pos, ref: str, alt: str) -> [str, int, int, str, str]:
    start, ref, alt = int(pos), ref.upper(), alt.upper()
    if len(ref) > 1 or len(alt) > 1 and ref != alt:
        if ref.startswith(alt) or ref.endswith(alt):
            if ref.startswith(alt):
                start = start + len(alt)
            ref = ref.replace(alt, '', 1)
            alt = ''
        elif alt.startswith(ref) or alt.endswith(ref):
            start = start + len(ref) - 1 if alt.startswith(ref) else start - len(alt) + len(ref)
            alt = alt.replace(ref, '', 1)
            ref = ''
        else:
            ref_rev, alt_rev, substr, stop, index = ref[::-1], alt[::-1], '', False, 0
            while index < len(ref) and index < len(alt):
                if ref_rev[index] != alt_rev[index]:
                    stop = True
                if ref_rev[index] == alt_rev[index] and not stop:
                    substr = ref_rev[index] + substr
                index += 1
            ref = re.sub(r'%s$' % substr, '', ref)
            alt = re.sub(r'%s$' % substr, '', alt)
            substr, stop, index = '', False, 0
            while index < len(ref) and index < len(alt):
                if ref[index] != alt[index]:
                    stop = True
                if ref[index] == alt[index] and not stop:
                    substr += ref[index]
                index += 1
            ref = re.sub(r'^%s' % substr, '', ref)
            alt = re.sub(r'^%s' % substr, '', alt)
            start += len(substr) - 1 if len(substr) and not ref else len(substr)
    end = start + len(ref) - 1 if ref else start
    return re.sub(r'[Cc][Hh][Cc]', '', chrom), start, end, ref if ref else '-', alt if alt else '-'

def run_maxentpy(vcf_file: str, genome_file: str, refgene_file: str, outfile: str):
    genome = Fasta(genome_file)
    refgenes = read_refgene(refgene_file)
    vcf_reader = pyvcf.Reader(filename=vcf_file)
    rows = list()
    for record in vcf_reader:
        for alt in record.ALT:
            chrom, start, end, ref, alt2 = vcf_to_av(record.CHROM, int(record.POS), str(record.REF), str(alt))
            result = None
            transcript_names = set()
            chrom2 = 'chrM' if chrom == 'MT' else f'chr{chrom}'
            for refgene in refgenes.query(f'Chrom == "{chrom2}" and TxStart <= {end} and TxEnd >= {start}').iloc:
                transcript = make_transcript(refgene)
                splicing = SplicingMaxEnt(chrom2, str(record.POS), str(record.REF), str(alt), transcript, genome)
                if result:
                    if abs(splicing.maxentscore_var) > abs(result.maxentscore_var):
                        result = splicing
                    elif abs(splicing.maxentscore_var) == abs(result.maxentscore_var):
                        if splicing.maxentscore_alt > result.maxentscore_alt:
                            result = splicing
                else:
                    result = splicing
            if result:
                rows.append({
                    'Chr': chrom, 'Start': start, 'End': end, 'Ref': ref, 'Alt': alt2,
                    'Maxent_type': result.splice_type,
                    'Maxent_pred': result.maxentpred,
                    'Maxent_score_ref': result.maxentscore_ref,
                    'Maxent_score_alt': result.maxentscore_alt,
                    'Maxent_score_var': result.maxentscore_var,
                    'Maxent_foldchange': result.maxentscore_foldchange
                })
                print(chrom, start, end, ref, alt, result.splice_type, result.maxentpred)
                
    with open(outfile, 'w') as fo:
        writer = csv.DictWriter(fo, fieldnames=[
            'Chr',
            'Start',
            'End',
            'Ref',
            'Alt',
            'Maxent_type',
            'Maxent_pred',
            'Maxent_score_ref',
            'Maxent_score_alt',
            'Maxent_score_var',
            'Maxent_foldchange'
        ], delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)