"""This module contains various utility functions.  The rule of thumb
for inclusion in this module is that the function could conceivably be
used in an entirely unrelated project."""

import os
import scipy.stats
from math import *

iota = lambda x:x #it is sometimes useful to have the identity
                  #function lying around

def sort_on(xs, ys, f):
    """Given f: x -> y, sort xs by corresponding ys"""
    zs = [(x, index(ys, lambda y: y == f(x))) for x in xs]
    return [tup[0] for tup in sorted(zs, key=lambda tup:tup[1])]

def head(xs, p=iota):
    """Take first element of xs, optionally satisfying predicate p"""
    return filter(p, xs)[0]

def index(xs, p):
    """Return index of first x satisfying p, or None if none do"""
    winners = filter(lambda (i, x): p(x), zip(range(len(xs)), xs))
    return winners[0][0] if winners else None

def choose2(xs):
    """return list of choose(xs, 2) pairs, retaining ordering on xs"""
    return [(x1, x2) for i, x1 in enumerate(xs) for x2 in xs[i+1:]]


def zipWith(f, xs, ys):
    """Return list of form [f(x, y)]"""
    return map(lambda (x, y): f(x, y), zip(xs, ys))

def mean(xs):
    if xs:
        return float(sum(xs))/len(xs)
    else:
        return None
    
def verbose_gen(xs,n=None):
    for i, x in enumerate(xs):
        if not n or i % n == 0:
            print i
        yield x

def cart_product(xs, ys):
    return [(x, y) for x in xs for y in ys]

def normalize(xs):
    total = float(sum(xs))
    return [x/total for x in xs]

def percentile(p, xs, upper=True):
    """Find the n% percentile cutoff for xs, i.e. smallest y such that
    P(x > y) <= n%"""
    zs = sorted(xs, reverse=not upper)
    n = len(zs)
    i = int(n * p)
    return zs[i]

def truncate(x,places=2):
    return int(x * 10**places)/float(10**places)

def mapmap(f,xxs):
    return [map(f,xs) for xs in xxs]

def variance(xs):
    return mean([x**2 for x in xs]) - mean(xs)**2

def sd(xs):
    return sqrt(mean([x**2 for x in xs]) - mean(xs)**2)

def fano(xs):
    return variance(xs)/mean(xs)

def dictmap(d,xs):
    """Map d over xs"""
    return [d[x] for x in xs]

def map_over_dict(f,d):
    """[k:v] -> [k:f(v)]"""
    fd = {x:y for (x,y) in map(f,zip(d.keys(),d.values()))}
    if type(d) is dict:
        return fd
    elif type(d) is type(defaultdict(list)):
        return defaultdict(lambda :f(d.default_factory()),fd)

def map_over_vals(f,d):
    """[k:v] -> [x:y], f(k,v) = (x,y)"""
    fd = {k:f(d[k]) for k in d}
    if type(d) is dict:
        return fd
    elif type(d) is type(defaultdict(list)):
        return defaultdict(lambda :f(d.default_factory()),fd)

def pearson(x, y):
    return scipy.stats.pearsonr(x, y)[0]
def spearman(x, y):
    return scipy.stats.spearmanr(x, y)[0]
def l2(x,y):
    return sum(zipWith(lambda u,v: (u-v)**2,x,y))
def length(x, y):
    return len(x)

def data2csv(data, filename, sep=", ",header=None,overwrite=False):
    make_line = lambda row: sep.join([str(field) for field in row]) + "\n"
    if filename in os.listdir('.') and not overwrite:
        print "found ",filename
        pass
    with open(filename, 'w') as f:
        if header:
            f.write(make_line(header))
        f.write("".join([make_line(row) for row in data]))

def get_name(obj):
    return head(dir(),lambda name: eval(name) == obj)

print("loaded utils")

def transpose(xxs):
    """Transpose a list of the form [[a1,a2...],[b1,b2..]...] into a
    list of the form [[a1,b1,...],[a2,b2,...]...]"""
    return zip(*xxs)

def group_codons(seq):
    return [seq[i*3:(i+1)*3] for i in range(len(seq)/3)]

def product(xs):
    return reduce(lambda x,y: x*y,xs)

def reverse(xs):
    return xs[::-1]

def org_matches_dir(org, org_dir):
    """Returns whether org_dir is the org_dir of org"""
    return all(word.lower() in org_dir.lower() for word in org.split('_'))
