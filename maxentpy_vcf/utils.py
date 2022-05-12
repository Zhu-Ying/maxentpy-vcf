import re


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
