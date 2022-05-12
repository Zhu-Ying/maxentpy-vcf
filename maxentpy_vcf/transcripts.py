from collections import namedtuple
import pandas as pd

TRANSCRIPT = namedtuple('TRANSCRIPT', ['name', 'chrom', 'tx_start', 'tx_end', 'strand', 'exons', 'introns'])


def read_refgene(infile: str) -> pd.DataFrame:
    dataframe = pd.read_csv(infile, header=None, delimiter="\t")
    dataframe.columns = [
        'Bin', 'Name', 'Chrom', 'Strand', 'TxStart', 'TxEnd', 'CdsStart', 'CdsEnd',
        'ExonCount', 'ExonStarts', 'ExonEnds', 'Score', 'Gene', 'CdsStartStat', 'CdsEndStat', 'ExonFrames'
    ]
    return dataframe


def make_transcript(refgene) -> TRANSCRIPT:
    """
    Iterate through a refGene file.

    GenePred extension format:
    http://genome.ucsc.edu/FAQ/FAQformat.html#GenePredExt

    Column definitions:
    0. uint undocumented id
    1. string name;                 "Name of gene (usually transcript_id from GTF)"
    2. string chrom;                "Chromosome name"
    3. char[1] strand;              "+ or - for strand"
    4. uint txStart;                "Transcription start position"
    5. uint txEnd;                  "Transcription end position"
    6. uint cdsStart;               "Coding region start"
    7. uint cdsEnd;                 "Coding region end"
    8. uint exonCount;              "Number of exons"
    9. uint[exonCount] exonStarts;  "Exon start positions"
    10. uint[exonCount] exonEnds;   "Exon end positions"
    11. int score;                  "Score"
    12. string name2;               "Alternate name (e.g. gene_id from GTF)"
    13. string cdsStartStat;        "enum('none','unk','incmpl','cmpl')"
    14. string cdsEndStat;          "enum('none','unk','incmpl','cmpl')"
    15. lstring exonFrames;         "Exon frame offsets {0,1,2}"
    """
    exon_starts = map(int, str(refgene.ExonStarts).strip(',').split(','))
    exon_ends = map(int, str(refgene.ExonEnds).strip(',').split(','))
    exons = list(zip(exon_starts, exon_ends))
    introns = [(exons[i][1], exons[i+1][0]) for i in range(int(refgene.ExonCount) - 1)]
    return TRANSCRIPT(
        name=str(refgene.Name),
        chrom=str(refgene.Chrom),
        tx_start=int(refgene.TxStart),
        tx_end=int(refgene.TxEnd),
        strand=str(refgene.Strand),
        exons=exons,
        introns=introns,
    )
