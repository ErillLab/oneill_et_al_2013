#!/usr/bin/env python

# Author: Sefa Kilic
# Modifications: PON
"""
This module contains functions for common tasks related to the
analysis of homologous coding sequences through reciprocal BLAST.

A typical workflow begins by constructing localized BLAST databases
for a list of organisms of interest.  Assume that this list is
represented as a standard Python list of strings called 'orgs'.  We
can construct the databases with the populate_dbs function:

populate_dbs(orgs)

which creates a database for each organism consisting of all amino
acid sequences contained in its FASTA file.  NB: this function assumes
the NCBI blast suite to be installed in the directory specified by
BLAST_BIN_PATH in the source.  The FASTA file is assumed to live in
data/$orgdir/, where $orgdir is the directory name in
ftp://ftp.ncbi.nih.gov/genomes/Bacteria/ from which the genome was
downloaded.  $orgdir is assumed to contain both the genus and species
name, but is otherwise free to vary.  The resulting database is built
in /blast_db.

Next, we wish to BLAST each protein in each genome against every other
genome.  In practice, the massively parallel nature of this task makes
cluster computing a more natural and efficient setting for this task.
The primary resource for parallel BLASTing can be found on
tara.rs.umbc.edu in the folder erill_common/cub, where the Python
script setup_slurm.py will automatically setup and process BLAST jobs
on tara's compute nodes.  In the anecdotal experience of the author,
the task of computing pairwise reciprocal BLASTs for ~20 organisms can
be accelerated from ~25hrs to ~10min by parallelizing, although
compressing, scp-ing and uncompressing the results files adds a
roughly O(log(n)) constant to the wall-time.  Since a typical job runs
in approximately 2 min, a reasonable rule of thumb is to parallelize
if running more than 15 jobs.  Otherwise, one can run the jobs locally with:

reciprocal_blasts2(orgs)

In any case, assume that the desired BLAST results are present in the
directory blast_results.  The next step is to collate the one-way
blast hits into reciprocal results.  We do this by calling:

collate_reciprocal_blasts(orgs)

which generates, for each pair of organisms, a tab-separated file
whose lines consist of a pair of locus tags, one from each organism,
if the locus tags are reciprocal blast hits.

One often wishes to find proteins conserved between more than two
organisms.  We solve this problem efficiently by constructing a graph
whose vertices are all locus tags for all organisms.  The graph
contains an edge between two tags if their respective proteins have
been identified as reciprocal blast hits.  The problem of identifying
proteins conserved over multiple organisms is therefore reduced to the
problem of finding n-vertex cliques, where n = |orgs|.  Although the
general complexity of this problem is daunting, it is quite tractable
in practice.  We construct the graph as follows:

all_results = load_reciprocals(orgs)
g = make_graph(all_results)

and find n-cliques by calling the find_full_cliques function:

cliques = find_full_cliques(g,orgs)

Ultimately, we will wish to record the results on disk.  We construct
an annotation dictionary of the form d[org][locus_tag] = annotation by calling:

anno_dict = annotation_dict(orgs)

We then construct locus tag dictionary of the form
d[org][gi] = locus_tag by calling:

named_ltds = make_named_ltds(orgs)

Now we can generate annotations for a group of cliques by computing:

clique_dict = analyze_cliques(cliques,orgs,anno_dict,named_ltds)

of the form clique_dict[i] = [annotations], where i is the index of
the clique in cliques and the list of annotations contains the
annotation for each locus tag in the clique.

Finally, the results are stored with:

cliques2csv(cliques,orgs,named_ltds,clique_dict,filename)

which represents the cliques in a tab-separated file.  A typical line
of this file consists of:

org1_locus_tag org1_index ... orgn_locus_tag orgn_index annotation

and has a header row giving the organism names.
"""
from math import *
import os
import time
import sys
import re
import csv
import operator
import string
import scipy.stats
import networkx as nx
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import Seq
from Bio import pairwise2

from orgs import *
from utils import *
from biochem import *
import tai

from collections import Counter, defaultdict
from string import Template
from matplotlib import pyplot as plt
import pylab

pairwise2.MAX_ALIGNMENTS = 1
BLAST_BIN_PATH = "/home/poneill/ncbi-blast-2.2.26+/bin/"
DB_PATH = 'blastdb'
SITE_LENGTH = 82
ORG_PATH = "data"


def org2nc_id(org):
    """Get NCid associated with org"""
    dir_name = org2dirname(org)
    file_name = head(os.listdir(os.path.join("data", dir_name)))
    (NCid, _) = os.path.splitext(file_name)
    return NCid

def org2dirname(org):
    """Get directory name associated with org"""
    dir_contents = os.listdir("data")
    dir_name = head([dc for dc in dir_contents if org_matches_dir(org, dc)])
    return dir_name
    
def makedb(sp):
    """Make database using fasta files of species"""
    os.system(BLAST_BIN_PATH + 'makeblastdb \
            -in %s -title %s -dbtype prot -out %s -logfile %s'
            % (sp['genome_file'], sp['name'], sp['db'], sp['name']+'.log'))


def blastn(seq, sp):
    """Run blastn on given seq and db"""
    qfile = 'query.tmp' #query file
    ofile = 'my_blast.xml' # output file
    with open(qfile, 'w') as f:
        f.write('%s' % seq)
    os.system(BLAST_BIN_PATH + 'blastn \
            -query %s -task blastn -db %s -out %s -outfmt 5'
            % (qfile, sp['db'], ofile))
    # parse output file
    result_handle = open(ofile)
    blast_record = NCBIXML.read(result_handle)
    # if no hit, return None
    if len(blast_record.alignments) == 0:
        return None
    # get first alignment
    alignment = blast_record.alignments[0]
    # make sure that hsps are sorted
    alignment.hsps.sort(key=lambda x: x.score, reverse=True)
    os.remove(qfile) # delete query file
    os.remove(ofile) # delete output file
    return alignment.hsps[0] # return best hit

def extend_site(genome, sbjct_start, sbjct_end, query_start, query_end):
    '''Extend site'''
    assert query_start < query_end
    if sbjct_start < sbjct_end: # strand:1
        if sbjct_start < 50 and sbjct_end > len(genome.seq)-50:
            return None
        ext_start = sbjct_start - (query_start-1) # extend start to left
        ext_end = sbjct_end + (SITE_LENGTH - query_end) # extend end to right
        retseq = genome.seq[ext_start-1 : ext_end]
    else: # strand:0
        if sbjct_start > len(genome.seq)-50 and sbjct_end < 50:
            return None
        ext_start = sbjct_start + (query_start-1)
        ext_end = sbjct_end - (SITE_LENGTH - query_end)
        retseq = genome.seq[ext_end-1 : ext_start].reverse_complement()
    return retseq
#


def populate_dbs(orgs):
    """given a list of organism names,  find directories,  grab files
    and construct databases"""
    org_dirs = os.listdir(ORG_PATH)
    file_ext = ".faa"
    for org in orgs:
        print "constructing db for ",  org
        org_dir = head([od for od in org_dirs if org_matches_dir(org, od)])
        full_org_dir = os.path.join(ORG_PATH, org_dir)
        fasta_file = head([f for f in os.listdir(full_org_dir)
                           if f.endswith(file_ext)])
        print "using fasta file", fasta_file
        sp = {'name':org, 
              'genome_file':os.path.join(full_org_dir, fasta_file), 
              'db': os.path.join(DB_PATH, org)}
        makedb(sp)
    
def reciprocal_blasts2(orgs):
    """Compute reciprocal blast hits for orgs"""
    org_dirs = os.listdir(ORG_PATH)
    results_contents = os.listdir('blast_results')
    file_ext = '.faa'
    for org1 in orgs:
        for org2 in orgs:
            print "starting on: ", org1, org2, "at", time.ctime()
            out_file = "results_%s_%s.txt" % (org1, org2)
            if org1 == org2:
                print "skipping", org1, org2
                continue 
            org_dir = head([od for od in org_dirs if org_matches_dir(org1, od)])
            full_org_dir = os.path.join(ORG_PATH, org_dir)
            fasta_file = head([f for f in os.listdir(full_org_dir)
                                            if f.endswith(file_ext)])
            full_fasta_file = os.path.join(full_org_dir, fasta_file)
            full_db_path = os.path.join(DB_PATH, org2)
            if out_file in results_contents:
                print "skipping", org1, org2
                continue
            os.system(os.path.join(BLAST_BIN_PATH, 
                      ('blastp -query %s -task blastp -db %s -out %s -outfmt 5'
                       % (full_fasta_file, full_db_path, 
                          os.path.join("blast_results", out_file)))))
            print "finished", org1, org2, "at", time.ctime()

def find_missing_results_files(orgs):
    """Print a list of expected files not found in blast_results"""
    dir_contents = os.listdir('blast_results')
    for org1 in (orgs):
        for org2 in orgs:
            if org1 == org2:
                continue
            filename1 = ("results_%s_%s.txt" % (org1, org2))
            filename2 = ("results_%s_%s.txt" % (org2, org1))
            for fn in [filename1, filename2]:
                if not fn in dir_contents:
                    print fn
                    

def collate_reciprocal_blasts(orgs):
    full_path = lambda f: os.path.join("blast_results", f)
    out_path = lambda f: os.path.join("reciprocal_results", f)
    locus_tag_dicts = {}
    dir_contents = os.listdir('reciprocal_results')
    for i, org1 in enumerate(orgs):
        for org2 in orgs[i+1:]:
            out_filename = "reciprocals_%s_%s.csv" % (org1, org2)
            out_file = out_path(out_filename)
            if out_filename in dir_contents:
                print "found file, moving on"
                continue
            print "collating recips for %s,%s at %s" % (org1, org2, time.ctime())
            if not org1 in locus_tag_dicts:
                locus_tag_dicts[org1] = make_locus_tag_dict(get_genome_filename(org1,"gbk"))
            if not org2 in locus_tag_dicts:
                locus_tag_dicts[org2] = make_locus_tag_dict(get_genome_filename(org2,"gbk"))
            org1_ltd, org2_ltd = locus_tag_dicts[org1], locus_tag_dicts[org2]
            filename1 = full_path("results_%s_%s.txt" % (org1, org2))
            filename2 = full_path("results_%s_%s.txt" % (org2, org1))
            print "parsing results for %s, %s" % (org1, org2)
            print "using:", filename1, filename2
            hits1 = parse_results(filename1, org1_ltd, org2_ltd)
            print "parsing results for %s, %s" % (org2, org1)
            hits2 = parse_results(filename2, org2_ltd, org1_ltd)
            print "finding reciprocals"
            reciprocals = find_reciprocals(hits1, hits2)
            print "writing"
            with open(out_file, 'w') as f:
                text = "\n".join(["\t".join(tup) for tup in reciprocals]) + "\n"
                f.write(text)

def make_locus_tag_dict(gbk_filename):
    print gbk_filename
    genome = SeqIO.read(gbk_filename, 'genbank')
    print "finding CDSs"
    CDSs = ([feature for feature in genome.features if feature.type == 'CDS'])
    print "findings gis"
    gis = [cds.qualifiers['db_xref'][0][3:] for cds in CDSs]
    print "finding locus tags"
    locus_tags = [cds.qualifiers['locus_tag'][0] for cds in CDSs]
    return {gi:locus_tag for (gi, locus_tag) in zip(gis, locus_tags)}

def make_named_ltds(orgs):
    return {org:make_locus_tag_dict(get_genome_filename(org,"gbk")) for org in orgs}

master_lt2o_dict = {}

def locus_tag2org(tag, named_ltds):
    if not tag in master_lt2o_dict:
        master_lt2o_dict[tag] = locus_tag2org_inner(tag, named_ltds)
    return master_lt2o_dict[tag]

def locus_tag2org_inner(tag, named_ltds):
    prefix = re.match(r"^[A-Za-z]+", tag).group(0)
    for org in named_ltds:
        first_tag = named_ltds[org].values()[0]
        if prefix == locus_tag_prefix(first_tag):
            return org
    #otherwise, search thoroughly
    for org in named_ltds:
        if any(value.startswith(tag) for value in named_ltds[org].values()):
            return org
    raise Exception("Couldn't find locus tag %s" % tag)

def get_genome_filename(org_name,ext):
    dir_contents = os.listdir('data')
    dir_name = os.path.join("data", 
                            head(filter(lambda d: org_matches_dir(org_name, d), 
                                        dir_contents)))
    fn = head(filter(lambda f: f.endswith('.' + ext), os.listdir(dir_name)))
    return os.path.join(dir_name, fn)
    
def parse_results(filename, query_locus_dict, target_locus_dict):
    """Accept a file containing the results of blasting query_org against
    target_org and return a dictionary of the form {protein in query_org: first
    blast hit in target_org}"""
    hits = {}
    query_def_pattern = re.compile(r"""<Iteration_query-def>
                                       gi\|(.*)\|ref.*
                                       </Iteration_query-def>""", re.X)
    hit_def_pattern = re.compile(r"<Hit_def>gi\|(.*?)\|ref.*?</Hit_def>")
    evalue_pattern =  re.compile(r"<Hsp_evalue>(.*?)</Hsp_evalue>")
    cutoff = 1e-10
    found_best_hit = False
    with open(filename) as f:
        for line in f:
            query_def_match = query_def_pattern.search(line)
            hit_def_match = hit_def_pattern.search(line)
            evalue_match = evalue_pattern.search(line)
            if query_def_match:
                query_name = query_def_match.group(1)
                found_best_hit = False
#                print "found query_name: %s" % query_name
            elif hit_def_match and not found_best_hit: 
                hit_name = hit_def_match.group(1)
#                print "found hit_name: %s" % hit_name
            elif evalue_match and not found_best_hit:
                evalue = float(evalue_match.group(1))
 #               print "found evalue: %s" % str(evalue)
                if evalue < cutoff:
                    query_locus = query_locus_dict[query_name]
                    target_locus = target_locus_dict[hit_name]
                    hits[query_locus] = (target_locus, cutoff)
                    found_best_hit = True 
    return hits

def find_reciprocals(d1, d2):
    """Given two dictionaries d1 and d1 collecting the matches
    for org1 against org2 and vice versa, return a list of tuples
    [(prot1, prot2)] such that prot1:prot2 is in hits1 and prot2:prot1
    is in hits2"""
    hits1 = {k:v[0] for k, v in d1.iteritems()}
    hits2 = {k:v[0] for k, v in d2.iteritems()}
    reciprocals = [(h1, h2) for h1 in hits1 for h2 in hits2
                   if h2 in hits1[h1] and h1 in hits2[h2]]
    return reciprocals

def load_reciprocals(orgs):
    """Read the reciprocal blast results and return a dictionary of form
    {results_filename:[(geneA, geneB)]}"""
#    dir_contents = os.listdir("reciprocal_results")
    def get(org1,org2):
        filename = "reciprocals_%s_%s.csv" % (org1,org2)
        lines = open(os.path.join("reciprocal_results", filename)).readlines()
        results = [re.search(r"(.*)\t(.*)\n", line).groups() for line in lines]
        return results
#    dir_contents = ["reciprocals_%s_%s.csv" % tup for tup in choose2(orgs)]
    all_results = {}
    for org1,org2 in choose2(orgs):
        #for filename in dir_contents:
        try:
            results = get(org1,org2)
        except IOError:
            reversed_results = get(org2,org1)
            results = map(reverse,reversed_results)
        filename = "reciprocals_%s_%s.csv" % (org1,org2)
        all_results[filename] = results
    return all_results

def make_graph(all_results):
    """Given the a dictionary of form
    {results_filename:[(geneA, geneB)]}, construct a graph whose
    vertices are genes and whose edges are reciprocal blast hits"""
    G = nx.Graph()
    vertices = set([tup[0] for result in all_results.values()
                    for tup in result])
    for v in vertices:
        G.add_node(v)
    for result in all_results.values():
        for tup in result:
            G.add_edge(*tup)
    return G

def get_genome(org):
    return SeqIO.read(get_genome_filename(org,"gbk"), 'genbank')

def gene_name2locus_tag(gene, genome):
    print gene
    cdss = ([feature for feature in genome.features
                 if feature.type == 'CDS'])
    for cds in cdss:
        if ('gene' in cds.qualifiers
            and gene in cds.qualifiers['gene']):
            return cds.qualifiers['locus_tag'][0]
    
def annotation_dict(orgs):
    """return a dictionary of form anno_dict[org][locus_tag] == annotation"""
    anno_dict = {}
    for org in orgs:
        print org
        anno_dict[org] = {}
        genome = get_genome(org)
        CDSs = ([feature for feature in genome.features
                 if feature.type == 'CDS'])
        for cds in CDSs:
            description = (head(cds.qualifiers['product'])
                           if 'product' in cds.qualifiers
                           else 'None')
            anno_dict[org][head(cds.qualifiers['locus_tag'])] = description
    return anno_dict

def analyze_cliques(cliques, orgs, anno_dict, named_ltds):
#    clique_dict = {i:[] for i, cl in enumerate(cliques)}
    clique_dict = {}
    for i, cl in enumerate(cliques):
        clique_dict[i] = [anno_dict[locus_tag2org(tag, named_ltds)][tag]
                          for tag in cl]
    return clique_dict


def get_cdss(genome):
    return ([feature for feature in genome.features if feature.type == 'CDS'])

memoize_cdss_from_org = True
cdss = {}
def cdss_from_org(org):
    if not org in cdss:
        cdss[org] = get_cdss(get_genome(org))
    return cdss[org]
        
def analyze_cliques2(cliques, orgs): # deprecated
    clique_dict = {i:[] for i, cl in enumerate(cliques)}
    for n, org in enumerate(orgs):
        print n, org
#        genome = SeqIO.read(get_genome_filename(org,"gbk"), 'genbank')
        genome = get_genome(org)
        print "computing cdss"
        CDSs = get_cdss(genome)
        print "searching cdss"
        for cds in CDSs:
#            print cds.qualifiers['locus_tag']
            for i, cl in enumerate(cliques):
                for locus_tag in cl:
                    if locus_tag in cds.qualifiers['locus_tag']:
                        clique_dict[i].append(head(cds.qualifiers['product']))
                        #print "found cds"
                        break
    return clique_dict

def maj_annotation(annotations):
    """Given a list of annotations (i.e. clique_dict[i]), return the
    most common"""
    return Counter(annotations).most_common(1)[0][0]

def find_full_cliques(G, orgs):
    return [cl for cl in nx.find_cliques(G) if len(cl) == len(orgs)]

def locus_tag_prefix(tag):
    return re.match("[A-Za-z]+", tag).group()

def cliques2csv2(cliques, orgs, named_ltds, clique_dict, filename): #deprecated
    sorted_tags = sort_on(cliques[0], orgs, 
                          lambda tag: locus_tag2org(tag, named_ltds))
    row_header = map(locus_tag_prefix, sorted_tags)
    indices = index_dicts(orgs,index)
    sorted_cliques = [sort_on(clique, row_header, locus_tag_prefix)
                      for clique in cliques]
    rows = [zipWith(lambda tag, org:(tag, indices[org][tag]), sc, orgs)
            for sc in sorted_cliques]
    annotations = [maj_annotation(clique_dict[i]) for i in range(len(cliques))]
    formatted_indices = ["\t".join(map(str, sum(row, ()))) for row in rows]
    formatted_rows = "\n".join(zipWith(lambda row, anno: row + '\t' + anno, 
                                       formatted_indices, annotations))
    formatted_header = ("\t".join(sum(zip(orgs, 
                                        ["" for i in range(len(orgs))]), ())) +
                        "\tAnnotation\n")
    with open(filename, 'w') as f:
        f.write(formatted_header)
        f.write(formatted_rows)


def locus_tag_prefixes(org, named_ltds):
    return set(locus_tag_prefix(tag) for tag in named_ltds[org].itervalues())

def cliques2csv(cliques, orgs, named_ltds, clique_dict, filename,index):
    print "finding all_prefixes"
    all_prefixes = set([locus_tag_prefix(tag) for clique in cliques
                        for tag in clique])
    print "sorting prefixes"
    sorted_prefixes = sort_on(all_prefixes, orgs, 
                              lambda tag: locus_tag2org(tag, named_ltds))
#    row_header = map(locus_tag_prefix, sorted_tags)
    print "indices"
    indices = index_dicts(orgs,index)
    print "sorted cliques"
    sorted_cliques = [sort_on(clique, sorted_prefixes, locus_tag_prefix)
                      for clique in cliques]
    print "rows"
    print orgs
    rows = [zipWith(lambda tag, org:(tag, indices[org][tag]), sc, orgs)
            for sc in sorted_cliques]
    annotations = [maj_annotation(clique_dict[i]) for i in range(len(cliques))]
    formatted_indices = ["\t".join(map(str, sum(row, ()))) for row in rows]
    formatted_rows = "\n".join(zipWith(lambda row, anno: row + '\t' + anno, 
                                       formatted_indices, annotations))
    formatted_header = ("\t".join(sum(zip(orgs, 
                                        ["" for i in range(len(orgs))]), ())) +
                        "\tAnnotation\n")
    with open(filename, 'w') as f:
        f.write(formatted_header)
        f.write(formatted_rows)

def generate_pairwise_spreadsheets(orgs, named_ltds, anno_dict):
    for org1, org2 in choose2(orgs):
        print "starting", org1, org2
        filename = "pairwise_correlations_%s_%s.csv" % (org1, org2)
        full_dir = os.path.join("clique_csvs", "pairwise_correlations")
        full_path = os.path.join(full_dir, filename)
        if filename in os.listdir(full_dir):
            print "found", filename, "skipping"
            continue
        pair = [org1, org2]
        print "making graph"
        g = make_graph(load_reciprocals(pair))
        print "finding cliques"
        cliques = find_full_cliques(g, pair)
        print "making clique_dict"
        clique_dict = analyze_cliques(cliques, pair, anno_dict, named_ltds)
        print "writing csv"
        cliques2csv(cliques, pair, named_ltds, clique_dict, full_path)

def generate_pairwise_spreadsheets2(orgs, named_ltds, anno_dict):#deprecated
    for org1, org2 in choose2(orgs):
        print "starting", org1, org2
        pair = [org1, org2]
        print "making graph"
        g = make_graph(load_reciprocals(pair))
        print "finding cliques"
        cliques = find_full_cliques(g, pair)
        print "making clique_dict"
        clique_dict = analyze_cliques2(cliques, pair)
        filename = "pairwise_correlations_%s_%s2.csv" % (org1, org2)
        print "writing csv"
        cliques2csv(cliques, pair, named_ltds, clique_dict, filename)

    
def read_genome_info_txt(filename):
    lines = open(filename).readlines()
    return [line.split("\t") for line in lines]

def find_indices_over_cliques(genome_info, cliques):
    return [float(line[0]) for clique in cliques
            for line in genome_info if line[1] in clique]

def nc_template(index):
    "%s-1.0_rcc_RCA" if index == "RCA" else "%s-1.0-tAI"
    
def correlate_index_over_cliques(cliques, orgs,index):
    """Return a list of lists containing, for each organism, the index value of
    that organism's contribution to each clique"""
    NCids = [org2nc_id(org) for org in orgs]
    
    dirs = [nc_template % nc for nc in NCids]
    genome_info_filenames = [os.path.join('index_results', d, 'genome_info.txt')
                             for d in dirs]
    return [find_indices_over_cliques(read_genome_info_txt(gif), cliques)
            for gif in genome_info_filenames]

def genome_info_filename(org,index):
    "return the filepath for genome_info.txt belonging to org"
    if index == "RCA" or index == "CAI":
        d = "%s-1.0_rcc_%s" % (org2nc_id(org), index)
    else:
        d = "%s-%s" % (org2nc_id(org), index)
    return os.path.join("index_results", d, "genome_info.txt")
    
def index_dict(org,index):
    """return a dictionary of the form {locus_tag:index} over org"""
    genome_info = read_genome_info_txt(genome_info_filename(org,index))
    return {row[1]:float(row[0]) for row in genome_info}

def index_dicts(orgs,index):
    return {org:index_dict(org,index) for org in orgs}

def check_correlation_in_pairwise_cliques():
    dirname = os.path.join("clique_csvs", "pairwise_correlations")
    fs = [f for f in os.listdir(dirname)
          if f.startswith("pairwise_correlations") and f.endswith(".csv")]
    results = []
    for f in fs:
        filename = os.path.join(dirname, f)
        lines = [line for line in csv.reader(open(filename), delimiter='\t')]
        xs = [float(line[1]) for line in lines[1:]]
        ys = [float(line[3]) for line in lines[1:]]
        results.append((f, scipy.stats.pearsonr(xs, ys)))
    return sorted(results, key=lambda result: result[1][0], reverse=True)

conditional_log_trasnform=True
def exp_dict(filename):
#    print filename
    conditional_log_transform = False
    MAX_LOG_VALUE = 100
    lines = [line for line in csv.reader(open(filename), delimiter=',')]
    d = {}
    for line in lines:
        locus_tag = line[0]
        data = map(float, line[1:])
        if not locus_tag in d:
            d[locus_tag] = data
        else:
            d[locus_tag].extend(data)
    if conditional_log_transform:
        vals = sum(d.values(),[])
        _max = max(vals)
        _min = abs(min(vals))
        def safe_log(x):
#            print x, _min
            return log(x + _min + 1)
        if _max > MAX_LOG_VALUE:
            d = {k:map(safe_log,d[k]) for k in d}
    return d


def pairwise_correlations(org1, org2):
    dirname = os.path.join("clique_csvs", "pairwise_correlations")
    f = "pairwise_correlations_%s_%s.csv" % (org1, org2)
    try:
        file_handle = open(os.path.join(dirname, f))
        rev = False
    except IOError:
        f = "pairwise_correlations_%s_%s.csv" % (org2, org1)
        file_handle = open(os.path.join(dirname, f))
        rev = True
    lines = [line for line in csv.reader(file_handle, delimiter='\t')]
    if not rev:
        pairs =  [(line[0], float(line[1]), line[2], float(line[3]))
                  for line in lines[1:]]
    else:
        pairs =  [(line[2], float(line[3]),line[0], float(line[1]))
                  for line in lines[1:]]
    return pairs

def exp_filename(org):
    return os.path.join("exp_csvs", org+"_exp.csv")
    
def expression_report_locals(org1, org2, p=0, upper=True):
    correlations = pairwise_correlations(org1, org2)
    (raw_exp_dict1, raw_exp_dict2) = map(lambda org:exp_dict(exp_filename(org)), 
                                        [org1, org2])
    exp_dict1_cutoff = percentile(p, raw_exp_dict1.values(), upper=upper)
    exp_dict2_cutoff = percentile(p, raw_exp_dict2.values(), upper=upper)
    meets_cutoff = operator.ge if upper else operator.le
    exp_dict1 = {k:raw_exp_dict1[k] for k in raw_exp_dict1
                 if meets_cutoff(raw_exp_dict1[k], exp_dict1_cutoff)}
    exp_dict2 = {k:raw_exp_dict2[k] for k in raw_exp_dict2
                 if meets_cutoff(raw_exp_dict2[k], exp_dict2_cutoff)}
    exp_data_for = lambda line: (line[0] in exp_dict1
                                 and line[2] in exp_dict2)
    exp_pairs = [map(mean, (exp_dict1[line[0]], exp_dict2[line[2]]))
                 for line in correlations
                 if exp_data_for(line)]
    (org1_exp, org2_exp) = zip(*exp_pairs) if len(exp_pairs) else ([], [])
    num_cliques = len(correlations)
    index_pairs = [(line[1], line[3]) for line in correlations
                  if exp_data_for(line)]
    (org1_indices, org2_indices) = zip(*index_pairs)
    org1_index_vs_org2_index = pearson(org1_indices, org2_indices)
    conserved_in_org1 = {line[0]:line[1] for line in correlations}
    conserved_in_org2 = {line[2]:line[3] for line in correlations}
    org1_exp_vs_org2_exp = pearson(org1_exp, org2_exp)
    num_expressed = len(exp_pairs)
    org1_index_vs_org1_exp = pearson(org1_indices, org1_exp)
    org2_index_vs_org2_exp = pearson(org2_indices, org2_exp)
    org1_index_vs_org2_exp = pearson(org1_indices, org2_exp)
    org1_exp_vs_org2_index = pearson(org2_indices, org1_exp)
    return locals() #sketchy!

def expression_vs_percentile_single(org, ps, upper=True):
    rs = [pearson(*expression_vs_percentile_single_data(org, p, upper))
          for p in ps]
    plt.plot(ps, rs)
    plt.show()

def expression_vs_percentile_single_data(org, index,p, upper=True):
    indices = index_dict(org,index)
    exp_path = lambda f: os.path.join("exp_csvs", f+"_exp.csv")
    exps = exp_dict(exp_path(org))
    cutoff = percentile(p, exps.values(), upper=upper)
    indices_exps = [(indices[tag], mean(exps[tag])) for tag in indices
                  if tag in exps and exps[tag] >= cutoff]
    final_indices, final_exps = zip(*indices_exps)
    return final_indices, final_exps
    
def expression_vs_percentile(org1, org2, ps, upper=True):
    vals = [expression_report_locals(org1, org2, p, upper=upper) for p in ps]
    variables = ["org1_index_vs_org2_index", 
                 "org1_exp_vs_org2_exp", 
                 "num_expressed", 
                 "org1_index_vs_org1_exp", 
                 "org2_index_vs_org2_exp", 
                 "org1_index_vs_org2_exp", 
                 "org1_exp_vs_org2_index"]
    for i, var in enumerate(variables):
        exec("%s_list=[v['%s'] for v in vals]" % (var, var))
        exec("p%s, =plt.plot(ps, %s_list)" % (i, var))
    plt.xlabel("Percentile cutoff")
    plt.ylabel("Pearson correlation")
    plt.legend([eval("p%s" % i) for i in range(len(variables))], variables, loc=3)
    plt.show()
    
def expression_report(org1, org2, p=0):
    d = expression_report_locals(org1, org2, p)
    template = Template("""Comparing: $org1, $org2
Found $num_cliques conserved proteins
Using $num_expressed found in both expression datasets
org1 index vs org2 index:\t\t$org1_index_vs_org2_index
org1 exp vs org2 exp:\t\t$org1_exp_vs_org2_exp
org1 self-expression:\t\t$org1_index_vs_org1_exp
org2 self-expression:\t\t$org2_index_vs_org2_exp
org1 index vs. org2 expression\t$org1_index_vs_org2_exp
org2 index vs. org1 expression\t$org1_exp_vs_org2_index
NB: Pearson r-values, not r^2
""")
    return template.substitute(d)

def all_expression_reports(orgs):
    for org1 in orgs:
        for org2 in orgs:
            try:
                print expression_report(org1, org2)
            except:
                continue

def expression_report_self_cross_comparison(orgs):
    selfs = defaultdict(list)
    crosses = defaultdict(list)
    for org1 in orgs:
        for org2 in orgs:
            print org1, org2
            try:
                print "trying"
                d = expression_report_locals(org1, org2)
                print "got dict for ", org1, org2
                selfs[org1].append(d["org%s_index_vs_org%s_exp" % (1, 1)])
                selfs[org2].append(d["org%s_index_vs_org%s_exp" % (2, 2)])
                print "extended selves"
                crosses[(org1, org2)].append(d["org1_index_vs_org2_exp"])
                crosses[(org2, org1)].append(d["org1_exp_vs_org2_index"])
                print "extended crosses"
            except:
                continue
    return selfs, crosses

def expression_report_graph(orgs):
    g = nx.MultiGraph()
    for org1 in orgs:
        for org2 in orgs:
            try:
                d = expression_report_locals(org1, org2)
                locals().update(d)
                org1_index = org1+"_index"
                org1_exp = org1+"_exp"
                org2_index = org2+"_index"
                org2_exp = org2+"_exp"
                if not org1 in g:
                    print "adding node", org1_index
                    g.add_node(org1_index)
                    print "adding node", org1_exp
                    g.add_node(org1_exp)
                if not org1 in g:
                    print "adding node", org2_index
                    g.add_node(org2_index)
                    print "adding node", org2_exp
                    g.add_node(org2_exp)
                var_pairs = choose2(["".join(pair)
                                     for pair in
                                     cart_product(["org1_", "org2_"], 
                                                  ["index", "exp"])])
                for var1, var2 in var_pairs:
                    label = eval("%s_vs_%s" % (var1, var2))
                    node1 = var1.replace("org1", org1).replace("org2", org2)
                    node2 = var2.replace("org2", org2).replace("org1", org1)
                    print "adding edge", node1, node2, label
                    g.add_edge(node1, node2, weight=label)
            except:
                pass
    return g

def collate_cdss_with_expression(org):
    genome = get_genome(org)
    cdss = get_cdss(genome)
    exp_d = exp_dict(exp_filename(org))
    return {cds.extract(genome).seq:exp_d[head(cds.qualifiers['locus_tag'])]
            for cds in verbose_gen(cdss)
            if head(cds.qualifiers['locus_tag']) in exp_d}

def tag2aas(tag, cdss, genome):
    print tag
    cds = head([cds for cds in cdss if tag in cds.qualifiers['locus_tag']])
    #drop symbol for stop codon
    return cds.extract(genome).seq.translate().tostring()[:-1] 

def analyze_index_vs_conservation(org1, org2, n=None, length_normalize=False):
    print "reading genomes"
    genome1, genome2 = map(get_genome, [org1, org2])
    print "finding CDSs"
    cdss1, cdss2 = map(get_cdss, [genome1, genome2])
    print "reading index dicts"
    index_dict1, index_dict2 = map(lambda org:index_dict(org,index), [org1, org2])
    print "reading correlations"
    correlations = pairwise_correlations(org1, org2)[:n]
    tag2aas1 = lambda tag: tag2aas(tag, cdss1, genome1)
    tag2aas2 = lambda tag: tag2aas(tag, cdss2, genome2)
    def nw_score(seq1, seq2):
        denom = float(seq1 + seq2)/2 if length_normalize else 1
        return head(pairwise2.align.globaldx(seq1, seq2, BLOSUM62))[2]/denom
    print "computing alignments"
    scores = [nw_score(tag2aas1(line[0]), tag2aas2(line[2]))
              for line in verbose_gen(correlations)]
    return [(line[0], index_dict1[line[0]], line[2], index_dict2[line[2]], scores[i])
            for i, line in enumerate(correlations)]

def test_gc(k, n):
    aa_alphabet = 'ABCDEFGHIKLMNPQRSTVWXYZ'
    for i in range(n):
        print i
        seq1 = "".join([random.choice(aa_alphabet) for i in range(k)])
        seq2 = "".join([random.choice(aa_alphabet) for i in range(k)])
        pairwise2.align.globaldx(seq1, seq2, BLOSUM62)

def parse_iterations_verbose(org):
    """parse iterationsVerbose.txt for org and return three
    dictionaries, the first of the form {codon:genome-wide codon
    freq}, the second of the form {codon:refset codon freq}, and
    the third of the form {codon:w-value}"""
    ncid = org2nc_id(org)
    filename = os.path.join("index_results", ncid+"-1.0_rcc_RCA", 
                            "iterationVerbose.txt")
    with open(filename) as f:
        lines = [filter(iota, line) for line in csv.reader(f, delimiter="\t")]
    genomic_freq_lines = lines[3:11] #magic numbers by inspection
    refset_freq_lines = lines[len(lines)-36:len(lines)-28]
    w_lines = lines[len(lines)-9:]
    def lines2dict(lines):
        d = {}
        for line in lines:
            for field in line:
                codon, val_string = field.split(':')
                val = float(val_string)
                d[codon] = val
        return d
    results = map(lines2dict, [genomic_freq_lines, refset_freq_lines, w_lines])
    return results

def all_pairwise_indices(orgs):
    return sorted([(all_correlations(org1, org2)[0], org1, org2)
                   for (org1, org2) in verbose_gen(choose2(orgs))], 
                  key = lambda tup: tup[0])

def get_codon_dicts(orgs,genomic=True):
    """Return a dictionary of codon usage biases index by org"""
    f = genomic_codon_freqs if genomic else refset_codon_freqs
    return {org:f(org) for org in orgs}

#helper functions...
def genomic_codon_freqs(org):
    return parse_iterations_verbose(org)[0]

def refset_codon_freqs(org):
    return parse_iterations_verbose(org)[1]

def all_correlations(org1, org2, n=None,correlation_func=pearson):
    """Compare the correlation between index values to the correlation
    between RCA weights"""
#    genome1, genome2 = map(get_genome, [org1, org2])
#    cdss1, cdss2 = map(get_cdss, [genome1, genome2])
#    cdss1, cdss2 = map(cdss_from_org, [org1,org2])
    index_dict1, index_dict2 = map(index_dict, [org1, org2])
    correlations = pairwise_correlations(org1, org2)[:n]
    index_table1 = [index_dict1[line[0]] for line in correlations]
    index_table2 = [index_dict2[line[2]] for line in correlations]
    gen_freqs1, refset_freqs1, ws1 = parse_iterations_verbose(org1)
    gen_freqs2, refset_freqs2, ws2 = parse_iterations_verbose(org2)
    gen_codon_table1 = [gen_freqs1[codon] for codon in codons]
    gen_codon_table2 = [gen_freqs2[codon] for codon in codons]
    refset_codon_table1 = [refset_freqs1[codon] for codon in codons]
    refset_codon_table2 = [refset_freqs2[codon] for codon in codons]
    w_table1 = [ws1[codon] for codon in codons]
    w_table2 = [ws2[codon] for codon in codons]
    results = map(lambda (x,y): correlation_func(x,y),[(index_table1, index_table2),
                                    (gen_codon_table1, gen_codon_table2),
                                    (refset_codon_table1, refset_codon_table2),
                                    (w_table1, w_table2)])
    index_corr, gen_codon_corr, refset_codon_corr, w_corr = results
    return (index_corr, gen_codon_corr, refset_codon_corr, w_corr)

def all_correlations_all_orgs(orgs,cf=pearson):
    d = all_correlations_all_orgs_dict(orgs,cf)
    return dictmap(d,choose2(orgs))
    
def all_correlations_all_orgs_dict(orgs,cf=pearson):
    return {(org1,org2):all_correlations(org1,org2,correlation_func = cf)
            for (org1,org2) in verbose_gen(choose2(orgs))}


def parse_ribo_str_crit(filename):
    with open(filename) as f:
        lines = f.readlines()
    results = [float(line.split(": ")[1]) for line in lines]
    (ribosomal_crit, strength_crit, content_crit) = results
    return (ribosomal_crit, strength_crit, content_crit)

def init_analysis_summary(orgs, outfile="initial_analysis_summary.csv"):
    def get_file(org, index):
        return os.path.join("index_results", org2nc_id(org) + "-1.0_rcc_" + index, 
                            "ribo_strength_criterion.txt")
    header = ("Organism, " +
              "RCA Ribosomal Crit, RCA Strength Crit, RCA Content Crit, "
              + "CAI Ribosomal Crit, CAI Strength Crit, CAI Content Crit")
    lines = [", ".join(map(str, (org, ) +
              parse_ribo_str_crit(get_file(org, "RCA")) +
              parse_ribo_str_crit(get_file(org, "CAI")))) for org in orgs]
    with open(outfile, 'w') as f:
        f.write(header)
        f.write("\n".join(lines))

def pairwise_summary(orgs, outfile, f=pearson):
    """Summarize the pairwise_correlation files by providing csv
    tables that report the Pearson and Spearman correlations, and
    number of genes, in each pairwise_correlation file. """
    header = ", "+", ".join(orgs)
    def get_lines(org1, org2):
        filename = os.path.join("clique_csvs", "pairwise_correlations", 
                                "pairwise_correlations_%s_%s.csv" % (org1, org2))
        with open(filename) as fn:
            #skim off the header
            lines = [line for line in csv.reader(fn, delimiter="\t")][1:] 
        return lines 
    def extract(org1, org2, content):
        if org1 == org2:
            return float(content == pearson or content == spearman)
        try:
            lines = get_lines(org1, org2)
        except IOError:
            try:
                lines = get_lines(org2, org1)
            except IOError:
                raise IOError(org1, org2)
        vals = [map(float, (line[1], line[3])) for line in lines]
        org1_vals, org2_vals = zip(*vals)
        return f(org1_vals, org2_vals)
    
    outlines = [[org1] + [extract(org1, org2, f)
                       for org2 in orgs] for org1 in orgs]
    with open(outfile, 'w') as outf:
        outf.write(header+"\n")
        outf.write("\n".join([", ".join(map(str, line)) for line in outlines]))

def escape_spaces(path):
    return path.replace(' ', '\ ')
def run_trna_scan(orgs):
    home_dir = "/home/poneill"
    trna_scan_path = os.path.join(home_dir, "bin")
    relative_home = os.getcwd()
    esc_relative_home = escape_spaces(relative_home)
    home_path = lambda fn: os.path.join(esc_relative_home, fn)
    trnas_path = lambda org: os.path.join(esc_relative_home, 
                                          "trnas", org + "_trnas")
    outfiles = os.listdir("trnas")
    for org in orgs:
        print org, time.ctime()
        print os.getcwd()
        org_dir = os.path.join("data", org2dirname(org))
        fn = head(os.listdir(org_dir), lambda f: f.endswith(".fna"))
        if org + "_trnas" in outfiles:
            print "skipping:", org
            continue
        fna_filename = home_path(os.path.join(org_dir, fn)) 
        outfile = trnas_path(org)
        os.chdir(trna_scan_path)
        command = "tRNAscan-SE -BH %s -f %s" % (fna_filename, outfile)
        print command
        os.system(command)
        os.chdir(relative_home)

def reverse_complement(seq):
    table = string.maketrans("ACGTacgt", 'TGCAtgca')
    return seq.translate(table)[::-1]

def org2trna_filename(org):
        return os.path.join("trnas", "%s_trnas" % org)

def trna_dict2trna_list(trna_dict,codon_order = "lexicographic"):
    lookup = {"lexicographic":codons,
              "sorted":sorted_codons,
              "table":table_codons}
    cs = lookup[codon_order]
    return [trna_dict[c] for c in cs]

def write_tai_data_for_R_script(org):
    print org
    w_list = w_list_for_R_script(org)
    orf_freqs = orf_frequencies_for_R_script(org)
    data2csv([w_list],os.path.join("tai_script",org+"_ws.csv"))
    data2csv(orf_freqs,os.path.join("tai_script",org+"_orf_freqs.csv"))
                  
def w_list_for_R_script(org):
    return trna_dict2trna_list(genomic_codon_freqs(org),codon_order = "table")

def orf_frequencies_for_R_script(org):
    sequences = coding_sequences(org)
    m = []
    for sequence in verbose_gen(sequences):
        m.append(trna_dict2trna_list(Counter(group_codons(sequence)),
                                     codon_order="table"))
    return m
    

def get_trna_dicts(orgs):
    """Return a dictionary of trna counts.
    Format: {org:{codon:[scores]}}"""
    return {org:parse_trnas(org2trna_filename(org))
            for org in orgs}

def parse_trnas(filename,cutoff=None):
#    type_pat = re.compile(r"Type: ([A-Za-z]{3})")
    anticodon_pat = re.compile(r"Anticodon: ([ATGC]{3})")
    score_pat = re.compile(r"Score: ([0-9.]+)")
    with open(filename) as f:
        lines = f.readlines()
    trnas = ["".join(lines[7*i:7*(i+1)]) for i in range(len(lines)/7)]
    codon_keys = defaultdict(list)
    for trna in trnas:
        if "Undet" in trna or "???" in trna or "pseudo" in trna:
            continue
        anticodon = re.search(anticodon_pat, trna).groups(0)[0]
        score = float(re.search(score_pat, trna).groups(0)[0])
        codon = reverse_complement(anticodon).lower()
        codon_keys[codon].append(score)
    return defaultdict(int,{codon:len(codon_keys[codon])
                             for codon in verbose_gen(codon_keys)})
#    return codon_keys

def my_extract(cds,genome):
    if cds.sub_features:
        return "".join([my_extract(sf,genome) for sf in cds.sub_features])
    else:
        seq = str(genome[cds.location.start:cds.location.end].seq)
        return (seq if cds.strand == 1
                else reverse_complement(seq))

def coding_sequences(org):
    try:
        fn = get_genome_filename(org,"ffn")
        cdss = SeqIO.parse(fn,'fasta')
        print "found ffn"
        return [cds.lower() for cds in cdss]
    except IndexError:
        print "found gbk"
        genome = get_genome(org)
        cdss = get_cdss(genome)
    return [str(cds.extract(genome).seq).lower() for cds in verbose_gen(cdss)]

def codon_usage(genome):
    codons = defaultdict(int)
    cdss = get_cdss(genome)
    old_objects = len(gc.get_objects())
    for cds in verbose_gen(cdss):
#        print "extracting"
#        seq = cds.extract(genome)
        seq = my_extract(cds,genome)
#       print "finished extracting"
        current_objects = len(gc.get_objects())
        new_objects = current_objects - old_objects
        for codon in group_codons(seq):
            codons[codon] += 1
        print "This: %d, New: %d, Garbage: %d, Collection Counts: %s"\
            % (current_objects, new_objects, len(gc.garbage), gc.get_count())
        old_objects = current_objects
    return codons

def codon_trna_correlation(codon_dicts,trna_dicts):
    assert(codon_dicts.keys() == trna_dicts.keys())
    orgs = codon_dicts.keys()
    for codon in codons:
        ts = [trna_dicts[org][codon] for org in orgs]
        print "ts:",ts
        cs = [codon_dicts[org][codon] for org in all_orgs]
        print "cs:",cs
        print codon,pearson(ts,cs)

def trna_analysis(gen_codon_dicts,refset_codon_dicts,trna_dicts):
    for org in all_orgs:
 	for codon in codons:
            gen_freq = gen_codon_dicts[org][codon]
            refset_freq = refset_codon_dicts[org][codon]
            copy_num = trna_dicts[org][codon]
            print "%s,%s,%s" % (refset_freq, gen_freq, copy_num)

def synonymous_codon_analysis(codon_dict,normalized=True):
    """Given a codon frequency dict, return a normcd alized dictionary of
    synonymous codons"""
    cond_norm = lambda xs:map(truncate,normalize(xs)) if normalized else iota
    return {aa:cond_norm([codon_dict[codon] for codon in translation_table[aa]])
            for aa in aas}
    
def codon_behavior(gen_codon_dicts,refset_codon_dicts,trna_dicts):
    for codon in codons:
        gen_freqs = [gen_codon_dicts[org][codon] for org in all_orgs]
        refset_freqs = [refset_codon_dicts[org][codon] for org in all_orgs]
        avg_cn = mean([trna_dicts[org][codon] for org in all_orgs])
        avg_diff = mean(zipWith(operator.sub,refset_freqs,gen_freqs))
        print codon,avg_diff,scipy.stats.wilcoxon(refset_freqs,gen_freqs)[1]<.05,avg_cn


def generate_correlations(): 
    def write_correlations(org_name,method):
        return data2csv(all_correlations_all_orgs(eval(org_name),method),
                        "%s_%s_correlations.csv" % (org_name,method.func_name),
                        header=["indices","genfreqs","refsetfreqs","ws"])
    # org_names = ["all_orgs","actinos","firmicutes","gammas","pseudos",
    #              "psychros","pvp_orgs"]
    org_names = ["pseudos","psychros"]
    methods = [pearson,spearman,l2]
    for org_name,method in cart_product(org_names,methods):
        write_correlations(org_name,method)

def normalize_codon_dict(d):
    """Take a dictionary of the form {codon:numeric_val} and return a
    dictionary normalized by amino acid bias"""
    print "hello"
    d_copy = d.copy()
    def normalization(codon):
        print codon
        denom = float(sum([d_copy[c] for c in synonymous_codons(codon)
                           if c in d_copy]))
        return denom if denom else 1
    print "defined"
    return defaultdict(int,{codon:d[codon]/normalization(codon) for codon in d})
    
def write_codon_trna_spreadsheet(filename):
    tt = translation_table
    sorted_codons = [codon for aa in tt for codon in tt[aa]]
    print "trna_dicts"
    trna_dicts = map_over_vals(normalize_codon_dict,get_trna_dicts(pvp_orgs))
    print "genomic_dicts"
    genomic_dicts = map_over_vals(normalize_codon_dict,get_codon_dicts(pvp_orgs,genomic=True))
    print "refset_dicts"
    refset_dicts = map_over_vals(normalize_codon_dict,get_codon_dicts(pvp_orgs,genomic=False))
    num_fields = 3 #refset freqs, genomic freqs, trna copy number
    rearrange_comma = lambda xs: "," + xs[:len(xs) - 1]
    aa_header = rearrange_comma("".join([aa + "," * len(tt[aa]) * num_fields
                                         for aa in tt]))
    codon_header = rearrange_comma("".join(codon + "," * num_fields
                                           for codon in sorted_codons))
    field_header = rearrange_comma("".join(["refset,genomic,trna,"
                                            for codon in sorted_codons]))
    def line(org):
        print org
        return [org] + sum([[refset_dicts[org][codon],
                             genomic_dicts[org][codon],
                         trna_dicts[org][codon]] for codon in sorted_codons],[])
    pseudo_lines = map(line,pseudos)
    psychro_lines = map(line,psychros)
    def moment_lines(lines,name,f):
        return [name] + map(lambda xs: truncate(f(xs),9),transpose(lines)[1:])#drop
                                                                                 #name
                                                                                 #column
    pseudo_avg = moment_lines(pseudo_lines,"Pseudomonas average",mean)
    psychro_avg = moment_lines(psychro_lines,"Psychrobacter average",mean)
    pseudo_std = moment_lines(pseudo_lines,"Pseudomonas standard dev",sd)
    psychro_std = moment_lines(psychro_lines,"Psychrobacter standard dev",sd)
    with open(filename,'w') as f:
        f.write(aa_header+"\n")
        f.write(codon_header+"\n")
        f.write(field_header+"\n")
        f.write("\n".join([",".join(map(str,line)) for line in pseudo_lines]) + "\n")
        f.write(",".join(map(str,pseudo_avg)) + "\n")
        f.write(",".join(map(str,pseudo_std)) + "\n")
        f.write("\n".join([",".join(map(str,line)) for line in psychro_lines]) + "\n")
        f.write(",".join(map(str,psychro_avg)) + "\n")
        f.write(",".join(map(str,psychro_std)) + "\n")
print "loaded"

def codon_trna_barchart(filename):
    with open(filename) as f:
        lines = [line for line in csv.reader(f, delimiter=",")]
    pseudo_averages = map(float,lines[12][1:])
    psychro_averages = map(float,lines[16][1:])
    pseudo_refset = [field for (i,field) in enumerate(pseudo_averages) if i % 3 == 0]
    psychro_refset = [field for (i,field) in enumerate(psychro_averages) if i % 3 == 0]
    pseudo_genomic = [field for (i,field) in enumerate(pseudo_averages) if i % 3 == 1]
    psychro_genomic = [field for (i,field) in enumerate(psychro_averages) if i % 3 == 1]
    pseudo_trnas = [field for (i,field) in enumerate(pseudo_averages) if i % 3 == 2]
    psychro_trnas = [field for (i,field) in enumerate(psychro_averages) if i % 3 == 2]
    ys = list(sum(zip(pseudo_refset,psychro_refset,pseudo_genomic,
             psychro_genomic,pseudo_trnas,psychro_trnas),()))
    fig = pylab.figure()
    ax = fig.add_subplot(1,1,1)
    indices = range(len(ys))
    sorted_codons = [codon for aa in tt for codon in translation_table[aa]]
    color_labels = ["%s %s %s" % (codon,species,typ) for codon in sorted_codons
                    for typ in ["refset","genomic","trnas"]
                    for species in ["pseudos","psychros"]]
    labels = ["%s (%s)" % (codon,translate(codon)) for codon in sorted_codons
              for typ in ["refset","genomic","trnas"]
              for species in ["pseudos","psychros"]]
    def color(label):
        if "genomic" in label and "pseudos" in label:
            return "magenta"
        elif "genomic" in label and "psychros" in label:
            return "cyan"
        elif "refset" in label and "pseudos" in label:
            return "red"
        elif "refset" in label and "psychros" in label:
            return "blue"
        elif "trna" in label and "pseudos" in label:
            return "yellow"
        elif "trna" in label and "psychros" in label:
            return "green"
        else:
            print(label)
    for index,y,label,color_label in zip(indices,ys,labels,color_labels):
        ax.bar(index,y,color=color(color_label))
    ax.set_xticks(indices)
    ax.set_xticklabels(labels,rotation='vertical')
    pylab.show()

def run_tai(orgs):
    print "trna dicts"
    trna_dicts = get_trna_dicts(orgs)
    print "sorting"
    sorted_trna_dicts = {org:[trna_dicts[org][codon]
                              for codon in sorted_codons]
                         for org in orgs}
    print "annotations"
    anno_dict = annotation_dict(orgs)
    for org in orgs:
        print "taing ",org
        folder = org2nc_id(org) + "-tAI"
        dirname = os.path.join("index_results",folder)
        if not folder in os.listdir('index_results'):
            print "making folder for ",org
            os.mkdir(dirname)
        filename = "genome_info.txt"
        if filename in os.listdir(dirname):
            print "found ",filename, "for ", org
            continue
        org_tai = tai.TAI(sorted_trna_dicts[org])
        genome = get_genome(org)
        cdss = get_cdss(genome)
        lines = []
        for cds in cdss:
            sequence = str(cds.extract(genome).seq)
            locus_tag = head(cds.qualifiers['locus_tag'])
            description = anno_dict[org][locus_tag]
            tai_value = org_tai.tai(sequence)
            lines.append([tai_value,locus_tag,description])
            print locus_tag
        data2csv(lines,os.path.join(dirname,filename))
        
def refset_info(org,index):
    nc = org2nc_id(org)
    dirname = nc + "-1.0_rcc_%s" % index
    full_path = os.path.join("index_results",dirname,"refset_info.txt")
    with open (full_path) as f:
        lines = [line for line in csv.reader(f,delimiter="\t")]
    return [(line[1],line[4]) for line in lines]

def top_scoring_annotations(org,index,n):
    nc = org2nc_id(org)
    dirname = nc + "-1.0_rcc_%s" % index
    full_path = os.path.join("index_results",dirname,"genome_info.txt")
    with open (full_path) as f:
        lines = [line for line in csv.reader(f,delimiter="\t")]
    return [line[5] for line in lines][:n]

def refset_annotations(org,index):
    return [line[-1] for line in refset_info(org,index)]
        
def refset_summary(orgs,filename):
    data = []
    for org in orgs:
        data.append([org])
        data.append([])
        for annotation in refset_annotations(org,"RCA"):
            data.append([annotation])
    data2csv(data,filename)

def top_scoring_summary(orgs,filename,n):
    data = []
    for org in orgs:
        data.append([org])
        data.append([])
        for annotation in top_scoring_annotations(org,"RCA",n):
            data.append([annotation])
    data2csv(data,filename)
