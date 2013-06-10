"""Code for benchmarking/validating CUB algorithms against expression
data.  Expression data taken from canonical expression spreadsheet
(Table_S1.csv, see gnxp.py)"""

from gnxp import master_exp_dict
from orgs import validation_orgs
from reciprocal_blast import get_cdss,get_genome_filename,get_genome,head
from scipy.stats import pearsonr, spearmanr
from utils import mean,transpose,se
import warnings

from milc import MILC

def benchmark_milc():
    #for org in validation_orgs:
    for org in ["Escherichia_coli","Mycobacterium_smegmatis"]:
        print org
        gbk_filename = get_genome_filename(org,'gbk')
        genome = get_genome(org)
        cdss = get_cdss(genome)
        milc = MILC(genome)
        org_exp_dict = master_exp_dict[org]
        milcs = []
        expss = [] # list of lists; one for each replicates
        for i,cds in enumerate(cdss):
            if i % 100 == 0:
                print i
            locus_tag = head(cds.qualifiers['locus_tag'])
            try:
                gene = head(cds.qualifiers['gene'])
            except:
                gene = None
            if locus_tag in org_exp_dict:
                exp_scores = org_exp_dict[locus_tag]
            elif gene in org_exp_dict:
                exp_scores = org_exp_dict[gene]
            try:
                seq = str(cds.extract(genome).seq)
                milc_score = milc.score(seq)
                expss.append(exp_scores)
                milcs.append(milc_score)
            except:
                print "tag failed:",locus_tag
        print "recovered milcs for: ",len(milcs),"tags"
        print "distinct milc scores:",len(set(milcs))
        milcs = [-x for x in milcs] #flip scores so they correlate
                                    #positively w/ expression
        spearmans = [spearmanr(milcs,exps)[0] for exps in transpose(expss)]
        spearman_mean = mean(spearmans)
        spearman_se = se(spearmans)
        print "Correlation:",org,spearman_mean,spearman_se,milc.num_ribosomals

print "loaded benchmark"
