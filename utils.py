"""This module contains various utility functions.  The rule of thumb
for inclusion in this module is that the function could conceivably be
used in an entirely unrelated project."""

import os
import scipy.stats
from math import *
from collections import Counter
import random

iota = lambda x:x #it is sometimes useful to have the identity
                  #function lying around

def sort_on(xs, ys, f):
    """Given f: x -> y, sort xs by corresponding ys"""
    zs = [(x, index(ys, lambda y: y == f(x))) for x in xs]
    return [tup[0] for tup in sorted(zs, key=lambda tup:tup[1])]

def head(xs, p=iota):
    """Take first element of xs, optionally satisfying predicate p"""
    filtered_xs = filter(p, xs)
    return filtered_xs[0] if filtered_xs else []

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

def variance(xs,correct=True):
    n = len(xs)
    correction = n/float(n-1) if correct else 1
    mu = mean(xs)
    return correction * mean([(x-mu)**2 for x in xs])

def sd(xs,correct=True):
    return sqrt(variance(xs,correct=correct))

def se(xs,correct=True):
    return sd(xs,correct)/sqrt(len(xs))

def geo_mean(xs):
    if xs:
        return product(xs) ** (1/float(len(xs)))
    
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
    print xs
    if xs:
        return sqrt(mean([x**2 for x in xs]) - mean(xs)**2)
    else:
        return None

def se(xs):
    if xs:
        return sd(xs)/sqrt(len(xs))
    else:
        return None
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

def maxZip(xs,ys):
    """Zip, padding the shorter list with Nones"""
    xlen = len(xs)
    ylen = len(ys)
    return zip(xs + [None for i in range(ylen - xlen)],
               ys + [None for i in range(ylen - xlen)])
    
def transpose(xxs):
    """Transpose a list of the form [[a1,a2...],[b1,b2..]...] into a
    list of the form [[a1,b1,...],[a2,b2,...]...]"""
    return zip(*xxs)

def group_codons(seq):
    return [seq[i*3:(i+1)*3] for i in range(len(seq)/3)]

def group_by(xs,n):
    chunks = [xs[i:i+n] for i in range(0,len(xs),n)]
    assert(xs == concat(chunks))
    return chunks

def product(xs):
    return reduce(lambda x,y: x*y,xs)

def reverse(xs):
    return xs[::-1]

def org_matches_dir(org, org_dir):
    """Returns whether org_dir is the org_dir of org"""
    return all(word.lower() in org_dir.lower() for word in org.split('_'))

def concat(xxs):
    return sum(xxs,[])

def plurality(xs):
    """Return the most frequently appearing x in xs"""
    return head(head(Counter(xs).most_common()))

def sample(xs,probs):
    probs = normalize(probs)
    partials = [sum(probs[:i+1]) for i in range(len(probs))]
    r = random.random()
    i = index(partials,lambda x: x > r)
    return xs[i]

def permute(xs):
    ys = xs[:]
    random.shuffle(ys)
    return ys

def base(b,n): #used for 1-based indexing in R
    return choose2(range(1,b + 1))[n - 1]

def pairs(xs):
    return zip([None] + xs,xs+[None])[1:len(xs)]

def ceiling(x):
    if x == floor(x):
        return x
    else:
        return floor(x + 1)
    
def chop(xs,k):
    """chop xs into k (approximately equal) contiguous sublists"""
    n = len(xs)
    #n = dk + r; last r items get spread out into last sublists
    d = n // k
    r = n % k
    print r,d
    xxs = []
    for i in range(k):
        xxs.append([xs.pop() for j in range(d + (i < r))])
    return xxs
        
def rollmean(xs,n):
    return [mean(xs[i:i+n]) for i in range(len(xs)-n +1)]

def split_on(xs, pred):
    """Split xs into a list of lists each beginning with the next x
    satisfying pred, except possibly the first"""
    indices = [i for (i,v) in enumerate(xs) if pred(v)]
    return [xs[i:j] for (i,j) in zip([0]+indices,indices+[len(xs)]) if i != j]

def catch(f):
    """Courtesy of Bryan Head on SO, slightly modified.  Usage:
    eggs = (1,3,0,3,2)
    [catch(lambda : 1/egg) for egg in eggs]
    [1, 0, None, 0, 0]"""
    try:
        return f()
    except Exception as e:
        return None

def argmax(xs):
    i,x = max(enumerate(xs),key= lambda (i,x):x)
    return i

def argmin(xs):
    i,x = min(enumerate(xs),key= lambda (i,x):x)
    return i

def rslice(xs,js):
    return [xs[j] for j in js]

def sample(n,xs,replace=True):
    if replace:
        return [random.choice(xs) for i in range(n)]
    else:
        ys = list(xs[:])
        samp = []
        for i in range(n):
            y = random.choice(ys)
            samp.append(y)
            ys.remove(y)
        return samp

def bs(xs):
    return sample(len(xs),xs,replace=True)

def sorted_indices(xs):
    """Return a list of indices that puts xs in sorted order.
    E.G.: sorted_indices([40,10,30,20]) => [1,3,2,0]"""
    return [i for (i,v) in sorted(enumerate(xs),key=lambda(i,v):v)]

def sliding_window(seq,w,verbose=False):
    i = 0
    n = len(seq)
    while i < n - w + 1:
        if verbose:
            if i % verbose == 0:
                print i
        yield seq[i:i+w]
        i += 1
