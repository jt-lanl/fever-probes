'''
hamming: compute hamming distance between two strings
that is: count number of positions in which the two strings differ
eg: hamming("ABCDEF","AGXDEF") returns 1
'''
import functools
import operator

@functools.lru_cache(maxsize=2**22)
def hamming(s: str, t: str) -> int:
    ## http://code.activestate.com/recipes/499304-hamming-distance/
    return sum(map(operator.ne, s, t))

