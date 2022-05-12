from pyfaidx import Fasta
from .transcripts import TRANSCRIPT
from .maxentpy import maxent
from .maxentpy.maxent import load_matrix5, load_matrix3


class SplicingMaxEnt:
    MATRIX5 = load_matrix5()
    MATRIX3 = load_matrix3()
    DONOR_EXON = 3
    DONOR_INTRON = 6
    ACCEPTOR_EXON = 3
    ACCEPTOR_INTRON = 20

    def __init__(self, chrom: str, pos: int, ref: str, alt: str, transcript: str, genome: Fasta):
        self.chrom, self.pos, self.ref, self.alt = chrom, int(pos), ref, alt
        self.splice_type, self.refseq, self.altseq = self.parse_splicing(self.chrom, self.pos, self.ref, self.alt, transcript, genome)
        self.maxentscore_ref = -1.0
        self.maxentscore_alt = -1.0
        self.maxentscore_foldchange = 1
        self.maxentscore_var = 0
        self.calculate_maxentscore()

    @staticmethod
    def reverse_complement(seq):
        """Retrun a reverse complementary seq"""
        nt_complement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C',
                         'a': 't', 'c': 'g', 't': 'a', 'g': 'c'}
        reverse_seq = list(reversed(seq))
        rev_comp_seq = [nt_complement[k] for k in reverse_seq]
        return ''.join(rev_comp_seq)

    @classmethod
    def format_donor(cls, raw_seq):
        cls.DONOR_EXON = 3
        format_seq = raw_seq[:cls.DONOR_EXON].lower() + \
            raw_seq[cls.DONOR_EXON:cls.DONOR_EXON + 2].upper() + \
            raw_seq[cls.DONOR_EXON + 2:].lower()
        return format_seq

    @classmethod
    def format_acceptor(cls, raw_seq):
        cls.ACCEPTOR_INTRON = 20
        format_seq = raw_seq[:cls.ACCEPTOR_INTRON - 2].lower() + \
            raw_seq[cls.ACCEPTOR_INTRON - 2:cls.ACCEPTOR_INTRON].upper() + \
            raw_seq[cls.ACCEPTOR_INTRON:].lower()
        return format_seq

    @classmethod
    def parse_splicing(cls, chrom: str, offset: int, ref: str, alt: str, transcript: TRANSCRIPT, genome: Fasta):
        """
        Get the refseq and altseq around splice sites
        """
        refseq, altseq = '', ''
        refseq_start, refseq_end = None, None
        splice_type = 'NA'
        for i, (intron_start, intron_end) in enumerate(transcript.introns):
            for offset_i in range(offset, offset+len(ref)):
                to_start = offset_i - intron_start
                to_end = offset_i - intron_end
                altseq_exon_end = altseq_intron_end = ''
                if transcript.strand == '+':
                    if -cls.DONOR_EXON < to_start <= cls.DONOR_INTRON:
                        splice_type = 'donor'
                        refseq_start = intron_start - cls.DONOR_EXON
                        refseq_end = intron_start + cls.DONOR_INTRON
                        refseq = genome[chrom][refseq_start:refseq_end].seq
                        if intron_start - cls.DONOR_EXON < offset - 1:
                            altseq_exon_end = genome[chrom][refseq_start:offset - 1].seq
                        if offset + len(ref) - 1 < refseq_end + len(ref) - len(alt):
                            altseq_intron_end = genome[chrom][offset + len(ref) - 1:
                                                              refseq_end + len(ref) - len(alt)].seq
                        altseq = altseq_exon_end + alt + altseq_intron_end
                    if -cls.ACCEPTOR_INTRON < to_end <= cls.ACCEPTOR_EXON:
                        splice_type = 'acceptor'
                        refseq_start = intron_end - cls.ACCEPTOR_INTRON
                        refseq_end = intron_end + cls.ACCEPTOR_EXON
                        refseq = genome[chrom][refseq_start:refseq_end].seq
                        if offset + len(ref) - 1 < refseq_end:
                            altseq_exon_end = genome[chrom][offset + len(ref) - 1:refseq_end].seq
                        if refseq_start - len(ref) + len(alt) < offset - 1:
                            altseq_intron_end = genome[chrom][refseq_start - len(ref) + len(alt):offset - 1].seq
                        altseq = altseq_intron_end + alt + altseq_exon_end
                else:
                    if -cls.ACCEPTOR_EXON < to_start <= cls.ACCEPTOR_INTRON:
                        splice_type = 'acceptor'
                        refseq_start = intron_start - cls.ACCEPTOR_EXON
                        refseq_end = intron_start + cls.ACCEPTOR_INTRON
                        refseq = genome[chrom][refseq_start:refseq_end].reverse.complement.seq
                        if refseq_start < offset - 1:
                            altseq_exon_end = genome[chrom][refseq_start:offset - 1].seq
                        if offset + len(ref) - 1 < refseq_end + len(ref) - len(alt):
                            altseq_intron_end = genome[chrom][offset + len(ref) - 1:
                                                              refseq_end + len(ref) - len(alt)].seq
                        altseq = altseq_exon_end + alt + altseq_intron_end
                        altseq = cls.reverse_complement(altseq)
                    if -cls.DONOR_INTRON < to_end <= cls.DONOR_EXON:
                        splice_type = 'donor'
                        refseq_start = intron_end - cls.DONOR_INTRON
                        refseq_end = intron_end + cls.DONOR_EXON
                        refseq = genome[chrom][refseq_start:refseq_end].reverse.complement.seq
                        if offset + len(ref) - 1 < refseq_end:
                            altseq_exon_end = genome[chrom][offset + len(ref) - 1:refseq_end].seq
                        if refseq_start - len(ref) + len(alt) < offset - 1:
                            altseq_intron_end = genome[chrom][refseq_start - len(ref) + len(alt):
                                                              offset - 1].seq
                        altseq = altseq_intron_end + alt + altseq_exon_end
                        altseq = cls.reverse_complement(altseq)

            # Format upper and lower case for better demonstration
            if splice_type == 'donor':
                refseq = cls.format_donor(refseq)
                altseq = cls.format_donor(altseq)
            elif splice_type == 'acceptor':
                refseq = cls.format_acceptor(refseq)
                altseq = cls.format_acceptor(altseq)
        return splice_type, refseq, altseq

    def calculate_maxentscore(self):
        """
        --- Calculate the maxentscan socre ---
        When a mutation occurs, if the WT score is above the threshold and
        the score variation (between WT and Mutant) is under -10% for HSF (-30% for MaxEnt)
        we consider that the mutation breaks the splice site.
        In the other case, if the WT score is under the threshold and
        the score variation is above +10% for HSF (+30% for MaxEnt) we consider that
        the mutation creates a new splice site.
        """
        maxentscore_alt = maxentscore_ref = -1.00
        if self.splice_type == 'donor':
            if len(self.refseq) == 9:
                maxentscore_ref = maxent.score5(self.refseq, matrix=self.MATRIX5)
            if len(self.altseq) == 9:
                maxentscore_alt = maxent.score5(self.altseq, matrix=self.MATRIX5)
        elif self.splice_type == 'acceptor':
            if len(self.refseq) == 23:
                maxentscore_ref = maxent.score3(self.refseq, matrix=self.MATRIX3)
            if len(self.altseq) == 23:
                maxentscore_alt = maxent.score3(self.altseq, matrix=self.MATRIX3)
        maxentscore_foldchange = maxentscore_alt / maxentscore_ref
        maxentscore_var = (maxentscore_alt - maxentscore_ref)/maxentscore_ref
        self.maxentscore_ref = round(maxentscore_ref, 2)
        self.maxentscore_alt = round(maxentscore_alt, 2)
        self.maxentscore_foldchange = round(maxentscore_foldchange, 2)
        self.maxentscore_var = round(maxentscore_var, 2)

    @property
    def maxentpred(self):
        if self.splice_type in ['donor', 'acceptor']:
            ref, alt, var = self.maxentscore_ref, self.maxentscore_alt, self.maxentscore_var
            return 'D' if (ref > 3 and var < -0.3) or (ref < 3 and abs(var) > 0.3 and alt > 3) else 'T'
        return 'U'
