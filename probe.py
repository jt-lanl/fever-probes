# routines for probesingle.py, probepair.py, probemulti.py

import re
import functools
import itertools

from hamming import hamming

@functools.lru_cache(maxsize=None)
def kmerify(k,s):
    ''' break a string s into a set of substrings of length k '''
    return set(s[i:i+k] for i in range(len(s)+1-k))

def coverage_offby(seqs,kmer_list,offby=0):
    if not kmer_list:
        return 0

    K = len(kmer_list[0])
    c = 0
    for s in seqs:
        kmer_set = kmerify(K,s.seq)
        if any( hamming(*xy)<= offby 
                for xy in itertools.product(kmer_list,kmer_set) ):
            c += 1
    return c

def kmer_seq_hamming(kmer,s):
    '''distance between a kmer and a sequence seq is 
    min hamming distance between kmer and all the k-mers in seq'''
    ## in case there's an exact match, test for that quickly
    if kmer in s.seq:
        return 0
        
    K = len(kmer)
    kmer_set =  kmerify(K,s.seq)
    return min( hamming(kmer,x) for x in kmer_set )

def uncovered_sequences_exact(seqs,kmer_list):
    '''list of sequences that are not covered by any of the kmers in kmer_list'''
    return [s for s in seqs if not any(k in s.seq for k in kmer_list)]

def max_distance(seqs,kmer_list):
    '''hamming distance of set kmer_list to furthest sequence'''
    subseqs = uncovered_sequences_exact(seqs,kmer_list)
    if len(subseqs)==0:
        return 0
    dist=[]
    for s in subseqs:
        dist_s = min( [kmer_seq_hamming(k,s) for k in kmer_list] )
        dist.append(dist_s)
    return max(dist)

def gc_content(s):
    ''' count the number of G's and C's in a probe string s '''
    return s.count("G") + s.count("C")


### Reverse Complement; eg, "AAACCCTG" -> "CAGGGTTT"

RC_Table = str.maketrans("ATCG","TAGC")
def reverse_complement(s):
    ''' return the reverse complement of string s:
    A <-> T, C <-> G, and order is revered '''
    ## Uses RC_Table which is defined globally, above
    return s.translate(RC_Table)[::-1]


### Hairpin Propensity -- find longest "stem" in a string 
### stem is a substring whose reverse complement is also in the string 

## Based on:
## https://stackoverflow.com/questions/43018230/how-to-find-the-longest-stem-for-a-stem-loop-in-dna

def longestStem(s):
    n = len(s)
    k = int(n/2) #length of longest possible stem
    candidate = ''
    i = 1

    while i <= k and len(candidate) == i - 1:
        for j in range(n-2*i+1):
            t = s[j:i+j]
            if reverse_complement(t) in s[i+j:]:
                candidate = t
                break
        i +=1
    return candidate

def lenStemLoop(s,stem):
    '''given a stem in a stem loop, what is the length of the loop'''
    ## ....stemxxxxxxxMETS....  (where METS is rev_complement of stem)
    ## return number of characters in xxxxx
    rcstem = reverse_complement(stem)
    ndx = s.index(stem)+len(stem)
    lenloop = s[ndx:].index(rcstem)
    #print(s.index(rcstem),s.index(stem),len(stem))
    return lenloop

def stemfound(s,length):
    ''' if a stem-loop found of given length, return the stem; else return blank '''
    n = len(s)
    candidate = ''
    for j in range(n-2*length+1):
        t = s[j:j+length]
        if reverse_complement(t) in s[j+length:]:
            candidate = t
            break

    return candidate



## _filter_kmers is a basic function, used by filter_kmerdict and filter_kmerset
## 
def _filter_kmers(kmers,remove_fcn,gc=0,stem=0):
    ## Remove kmers with low GC content
    if gc:
        k_low_gc=[k for k in kmers if gc_content(k) < gc]
        list(map( remove_fcn, k_low_gc ))
    ## Remove kmers with long self-complimentary patterns "stems"
    if stem:
        k_stem=[k for k in kmers if stemfound(k,length=1+stem)]
        list(map( remove_fcn, k_stem ))
    return kmers

def filter_kmerset(kmers,**kwfilter):
    ''' return set of kmers with those failing GC and HPP conditions removed '''
    return _filter_kmers(kmers,kmers.remove,**kwfilter)

def filter_kmerdict(kmers,**kwfilter):
    ''' return dict of kmers with those failing GC and HPP conditions removed '''
    return _filter_kmers(kmers,kmers.pop,**kwfilter)

def old_filter_kmerdict(kmers,gc=0,stem=0):
    ## Remove kmers with low GC content
    if gc:
        k_low_gc=[k for k in kmers if gc_content(k) < gc]
        list(map( kmers.pop, k_low_gc ))
    ## Remove kmers with long self-complimentary patterns "stems"
    if stem:
        k_stem=[k for k in kmers if stemfound(k,length=1+stem)]
        list(map( kmers.pop, k_stem ))
    return kmers


def filter_seq_names_by_pattern(seqs,pattern):
    ''' return a list of sequences whose names match the given pattern '''
    patt = re.compile(pattern)
    return [s for s in seqs if patt.search(s.name)]
    
