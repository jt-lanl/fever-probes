
import sys
from collections import Counter
from hamming import hamming

'''
Epitope = string of length E (typically E=9)

OFFBY:

Algorithms for finding (and counting) all the epitopes in a set
that differ by a fixed number of characters from a given epitope

These algorithms are based on the fact that if two strings of length E
differ by at most 'offby' characters, then they will share a common 
contiguous substring of length s.  The formula
    s=int((E-1)/offby) 
is based on considering wrap-around strings to be contiguous; eg
       ABCDEFGHI and ABxDEFxHI 
share the common contiguous substring "HIAB": AB.....HI
In this example, E=9, offby=2, and s=4
'''

def offby_struct_prepare(E,fullepiset,offby):
    '''
    Set up data structure ldict in time O(E*N) 
       where N is number of epitopes in fullepiset
       and E is the length of the epitope
    ldict is a list of E dict()'s
    ldict[k] is a dict() whose keys are substrings of length s
    and whose value for a given key is a set of epitopes that match
    the substring at the k'th position.
    eg: ldict[3][ABCD] = {xxxABCDxx, ... }
    or: ldict[4][ABCD] = {xxxxABCDx, ... }
    or: ldict[7][ABCD] = {CDxxxxxAB, ... } ## Note wrap-around
    '''

    s=int((E-1)/offby) ## length of substring

    ldict = [ dict() for k in range(E) ]
    for xepi in fullepiset:
        epiepi = xepi+xepi ## enables wrap-around
        for k in range(E):
            ldict[k].setdefault(epiepi[k:k+s], set() ).add(xepi)

    return s,ldict

def offby_struct_candidates(E,s,ldict,epi):
    '''
    Returns a set of epitopes that are candidates for 
    being within 'offby' characters of 'epi'

    Runtime to setup candidate list is O(E*N)
       where N is number of epitopes in subepiset
    Runtime to screen candidates is O(E*N*M)
       where O(E) is time to run hamming
       where M is average size of candidate set oblist[epi]
    (The reason we want substring length s as large as possible 
     is to reduce the size of this candidate pool)
    (The reason we want E separate sets for each s-mer 
     is also to reduce the size of the candidate pool; eg,
     we don't want xxxABCDxx to be a candidate match for xxxxABCDx)
    '''

    candi = set()
    epiepi = epi+epi 
    for k in range(E):
        try:
            candi.update(ldict[k][epiepi[k:k+s]])
        except KeyError:
            #print("epi=",epi,"may not be in the sequence set")
            pass

    return candi

def mkoblist(E,fullepiset,subepiset=None,offby=1):
    """Given a container (set, list, etc) of epitopes (epi's); produce a
    dict() called oblist such that oblist[epi] is a set of all the
    epitopes whose hamming distance from epi is less than or equal to
    offby.  Do NOT include epi in that list.

    if subepiset is supplied, it should be a subset of fullepiset, and
    the dict().keys will only include the subpeiset epi's, but the
    dict().values will be sets that collectively include all the epi's
    in fullepiset

    """
    ## Return oblist, a dict() whose keys are epitopes in subepiset
    ## oblist[epi] is set()
    ## of epitopes xepi such that 1 <= hamming(epi,xepi) <= offby
    ## with xepi taken from fullepiset

    oblist=dict()
    if subepiset is None:
        subepiset=fullepiset

    if offby==0:
        for epi in subepiset:
            oblist[epi]=set()
        return oblist ## return full dict() of emptysets

    s,ldict = offby_struct_prepare(E,fullepiset,offby)

    for epi in subepiset:
        oblist[epi]=offby_struct_candidates(E,s,ldict,epi)
        
        ## At this point oblist[epi] contains candidates xepi
        ## that contain a continguous substring match with epi
        ## This is a superset of xepi's that have hamming(epi,xepi)<=offby
        ## So we'll remove xepi's for which hamming > offby
        ## (Also, we remove epi from oblist, so hamming >= 1)

        #oblist[epi].remove(epi)
        #oblist[epi].discard(epi) ## for situations in which epi is not in oblist[epi]
        xepi_toofar = [ xepi for xepi in oblist[epi] 
                        if hamming(xepi,epi) > offby ]
        for xepi in xepi_toofar:
            oblist[epi].remove(xepi)

    return oblist

