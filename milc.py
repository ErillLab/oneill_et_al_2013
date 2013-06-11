"""The purpose of this script is to compute MILC scores (see Supek and
Vlahovicek 2005) for protein-coding genes, given a reference set
computed from a set of ribosomal proteins containing at least 100aas.
The MILC score for a given amino acid a is given by:
M_a = 2*sum_c ln(f_c/g_c)

where f_c, g_c are respectively the observed frequency codon c in the
gene, and g_c the frequency in the reference set.  The sum is taken
over all synonymous codons.

The MILC score for a gene is just the sum over all amino acids,
i.e. all 20 amino acids, not all amino acids encoded by the sequence.
"""
from utils import *
from biochem import *
from collections import Counter
from Bio import SeqIO
from reciprocal_blast import get_cdss

class MILC(object):
    def __init__(self,genome):
        """Initialize a MILC instance from a genbank filename,
        extracting ribosomal proteins therefrom."""
        CDSs = get_cdss(genome)
        print "identifying ribosomal proteins"
        ribosomals = [cds for cds in CDSs
                      if 'product' in cds.qualifiers
                      and any(("ribosomal subunit" in product
                               or "ribosomal protein" in product)
                              for product in
                              cds.qualifiers['product'])]
        cutoff = 100 # throw away sequences < 100aa
        ribosomal_seqs = [str(r.extract(genome).seq).lower() for r in ribosomals
                          if len(r.extract(genome)) >= cutoff]
        assert all(len(r) % 3 == 0 for r in ribosomal_seqs)
        counts = Counter(group_codons("".join(ribosomal_seqs)))
        coding_counts = {k:v for k,v in counts.items() if not k in stop_codons}
        z = float(sum(coding_counts.values()))
        self.refset_freqs = {k:v/z for k,v in coding_counts.items()}
        self.num_ribosomals = len(ribosomals)

    def score(self,seq):
        """Score a coding sequence, assuming sum should be taken over
        amino acid alphabet"""
        seq = seq.lower()
        aa_seq = translate(seq)
        used_aas = set(aa_seq)
        counts = Counter(group_codons(seq))
        coding_counts = {k:v for k,v in counts.items() if not k in stop_codons}
        z = float(sum(coding_counts.values()))
        seq_freqs = {k:v/z for k,v in coding_counts.items()}
        def Ma(aa):
            return sum(2*counts[codon]*log(seq_freqs[codon]/
                                           self.refset_freqs[codon])
                       for codon in translation_table[aa]
                       if codon in seq_freqs and codon in self.refset_freqs)
        M = sum([Ma(aa) for aa in translation_table.keys()])
        L = float(len(seq)/3)
        C = sum(len(redundancy) - 1
                for redundancy in aa_seq)/L - 0.5
        return M/L - C

