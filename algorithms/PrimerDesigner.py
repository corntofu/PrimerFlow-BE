import sqlite3
import numpy as np
import pysam
from typing import List, Dict, Tuple

############################################
# Basic utilities
############################################

RC_MAP = str.maketrans("ATCG", "TAGC")

def reverse_complement(seq: str) -> str:
    return seq.translate(RC_MAP)[::-1]

############################################
# Thermodynamics (SantaLucia 1998 simplified)
############################################

NN_PARAMS = {
    'AA': (-7.9, -22.2), 'TT': (-7.9, -22.2),
    'AT': (-7.2, -20.4), 'TA': (-7.2, -21.3),
    'CA': (-8.5, -22.7), 'TG': (-8.5, -22.7),
    'GT': (-8.4, -22.4), 'AC': (-8.4, -22.4),
    'CT': (-7.8, -21.0), 'AG': (-7.8, -21.0),
    'GA': (-8.2, -22.2), 'TC': (-8.2, -22.2),
    'CG': (-10.6, -27.2), 'GC': (-9.8, -24.4),
    'GG': (-8.0, -19.9), 'CC': (-8.0, -19.9)
}

def calc_tm_nn(seq: str, dna_nM=50.0, salt_mM=50.0) -> float:
    dh, ds = 0.0, 0.0
    for i in range(len(seq) - 1):
        h, s = NN_PARAMS.get(seq[i:i+2], (0, 0))
        dh += h
        ds += s

    ds += 0.368 * (len(seq) - 1) * np.log(salt_mM / 1000.0)
    R = 1.987
    c = dna_nM * 1e-9
    tm = (dh * 1000) / (ds + R * np.log(c / 4))
    return tm - 273.15

def calc_3prime_dG(seq: str, n: int = 5) -> float:
    dg = 0.0
    tail = seq[-n:]
    for i in range(len(tail) - 1):
        h, s = NN_PARAMS.get(tail[i:i+2], (0, 0))
        dg += h
    return dg

############################################
# Needleman–Wunsch
############################################

def nw_3prime_mismatch(a: str, b: str, k=5) -> int:
    return sum(x != y for x, y in zip(a[-k:], b[-k:]))

############################################
# Main Designer
############################################

class PrimerDesigner:
    def __init__(
        self,
        genome_fasta: str,
        annotation_db: str
    ):
        self.genome = pysam.FastaFile(genome_fasta)
        self.db = sqlite3.connect(annotation_db)
        self.cur = self.db.cursor()

    ########################################
    # Stage 2.1 – Sliding window + filtering
    ########################################

    def generate_candidates(
        self,
        template: str,
        k_min=18,
        k_max=25,
        tm_range=(57, 63)
    ) -> List[Dict]:

        candidates = []

        for k in range(k_min, k_max + 1):
            for i in range(len(template) - k + 1):
                seq = template[i:i+k]
                rc = reverse_complement(seq)

                for strand, s in [('+', seq), ('-', rc)]:
                    tm = calc_tm_nn(s)
                    if not (tm_range[0] <= tm <= tm_range[1]):
                        continue

                    dg3 = calc_3prime_dG(s)
                    if dg3 <= -10.0:
                        continue

                    if s[-1] not in "GC":
                        continue

                    if s[-5:].count('G') + s[-5:].count('C') > 4:
                        continue

                    candidates.append({
                        'seq': s,
                        'start': i,
                        'end': i + k - 1,
                        'strand': strand,
                        'tm': tm,
                        'dg3': dg3
                    })

        return candidates

    ########################################
    # Stage 2.2 – Local DB filtering
    ########################################

    def local_db_filter(
        self,
        chrom: str,
        primer: Dict
    ) -> bool:

        # SNP (3' end strict)
        s = primer['end'] - 4
        e = primer['end']
        self.cur.execute(
            "SELECT COUNT(*) FROM snp WHERE chrom=? AND pos BETWEEN ? AND ?",
            (chrom, s, e)
        )
        if self.cur.fetchone()[0] > 0:
            return False

        # Restriction enzyme overlap
        self.cur.execute(
            "SELECT COUNT(*) FROM restriction_site "
            "WHERE chrom=? AND NOT (end < ? OR start > ?)",
            (chrom, primer['start'], primer['end'])
        )
        if self.cur.fetchone()[0] > 0:
            return False

        # Exon junction span
        self.cur.execute(
            "SELECT start, end FROM exon WHERE chrom=?",
            (chrom,)
        )
        exons = self.cur.fetchall()
        for (a_end, b_start) in zip(exons[:-1], exons[1:]):
            if primer['start'] <= a_end[1] and primer['end'] >= b_start[0]:
                return False

        return True

    ########################################
    # Stage 2.3 – Genome specificity
    ########################################

    def specificity_check(
        self,
        chrom: str,
        primer: Dict,
        max_hits=50,
        mismatch_cutoff=2
    ) -> bool:

        hits = 0
        for ref in self.genome.references:
            seq = self.genome.fetch(ref)
            pos = seq.find(primer['seq'])
            while pos != -1:
                if ref == chrom:
                    pos = seq.find(primer['seq'], pos + 1)
                    continue

                off = seq[pos:pos+len(primer['seq'])]
                mm = nw_3prime_mismatch(primer['seq'], off)
                if mm < mismatch_cutoff:
                    return False

                hits += 1
                if hits > max_hits:
                    return False

                pos = seq.find(primer['seq'], pos + 1)

        return True

    ########################################
    # Stage 2.4 – Pairing & scoring
    ########################################

    def pair_primers(
        self,
        primers: List[Dict],
        product_range=(100, 300),
        max_tm_diff=3.0
    ) -> List[Dict]:

        fwd = [p for p in primers if p['strand'] == '+']
        rev = [p for p in primers if p['strand'] == '-']

        pairs = []

        for f in fwd:
            for r in rev:
                size = r['start'] - f['end']
                if not (product_range[0] <= size <= product_range[1]):
                    continue

                if abs(f['tm'] - r['tm']) > max_tm_diff:
                    continue

                penalty = (
                    abs(f['tm'] - 60) +
                    abs(r['tm'] - 60) +
                    abs(f['dg3']) +
                    abs(r['dg3'])
                )

                pairs.append({
                    'fwd': f,
                    'rev': r,
                    'product_size': size,
                    'penalty': penalty
                })

        return sorted(pairs, key=lambda x: x['penalty'])
