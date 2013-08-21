"""Code for benchmarking/validating CUB algorithms against expression
data.  Expression data taken from canonical expression spreadsheet
(Table_S1.csv, see gnxp.py)"""

from gnxp import master_exp_dict
from orgs import validation_orgs
from reciprocal_blast import get_cdss,get_genome_filename,get_genome,head,org2nc_id
from scipy.stats import pearsonr, spearmanr
from utils import mean,transpose,se
import warnings
import re,os

from milc import MILC
from delta import Delta


def benchmark_cat():
    """Compute correlation with expression for the CDC method of Zhang
    BMC Bioinformatics 2012"""
    for org in validation_orgs:
        print org
        try:
            gbk_filename = get_genome_filename(org,'gbk')
            genome = get_genome(org)
            cdss = get_cdss(genome)
            ncid = org2nc_id(org)
            cat_filename = os.path.join("index_results",ncid+"_CAT",ncid+".cat")
            with open(cat_filename) as f:
                lines = [line.split("\t") for line in f.readlines()[1:]]
            labels,cdcs = transpose([(fields[0],fields[10]) for fields in lines])
            matches = [re.search(r":(c?)(\d+)-(\d+)",label).groups() for label in labels]
            locations = [((int(start),int(stop))
                          if c == ''
                          else (int(stop) - 1,int(start)))
                         for (c,start,stop) in matches]
            cat_dict = {location:float(cdc) for location,cdc in zip(locations,cdcs)}
            org_exp_dict = master_exp_dict[org]
        # a dictionary of form {(start,stop):[locus tags]}
            location2lt = {(feature.location.start+1,feature.location.end):
                               feature.qualifiers['locus_tag'][0]
                           for feature in genome.features
                           if ('locus_tag' in feature.qualifiers)}
            correlates = [(cdc,org_exp_dict[location2lt[location]])
                          for location,cdc in cat_dict.items()
                          if location in location2lt
                          and location2lt[location] in org_exp_dict]
            cdcs,exps = transpose(correlates)
            
            rhos = [spearmanr(cdcs,map(lambda xs:xs[i],exps))[0]
                    for i in range(len(exps[0]))]
            print "num correlates:",len(correlates)
            print "Correlation:",org,mean(rhos),se(rhos)
        except:
            print "Failed on:",org
        

def benchmark_method(method="delta"):     # or MILC
    for org in ["Mycobacterium_smegmatis"]:#problem_species:#validation_orgs:
        print org
        gbk_filename = get_genome_filename(org,'gbk')
        genome = get_genome(org)
        cdss = get_cdss(genome)
        Method = MILC if method == "MILC" else Delta
        method = Method(org)
        org_exp_dict = master_exp_dict[org]
        scores = []
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
                method_score = method.score(seq)
                expss.append(exp_scores)
                scores.append(method_score)
            except:
                print "tag failed:",locus_tag
        print "recovered scores for: ",len(scores),"tags"
        print "distinct scores:",len(set(scores))
        if method == "MILC":
            scores = [-x for x in scores] #flip scores so they correlate
                                    #positively w/ expression
        spearmans = [spearmanr(scores,exps)[0] for exps in transpose(expss)]
        spearman_mean = mean(spearmans)
        spearman_se = se(spearmans)
        print "Correlation:",org,spearman_mean,spearman_se#,milc.num_ribosomals



print "loaded benchmark"
