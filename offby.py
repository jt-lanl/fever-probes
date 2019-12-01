from typing import List, Set, Dict, Tuple

from collections import Counter
from functools import lru_cache

from hamming import hamming
from probe import kmerify as probe_kmerify
from mkoblist import mkoblist

class kmeroffby():

    def __init__(self, K: int, offby: int) -> None:
        self.K = K
        self.offby = offby
        self.prepared = False

    def all_kmers_exact(self, seqs: List) -> Set[str]:
        kmers: Set[str] = set()
        for s in seqs:
            kmers.update( probe_kmerify(self.K, s.seq) )
        return kmers

    def all_kmers_wcount_exact(self, seqs: List, pool: Set = None) -> Dict[str,int]:
        ''' 
        return a Counter such that cnt[kmer] is number of seqs in which kmer appear;
        if pool is specified, then only include kmers in pool (saves space, but not time)
        '''
        cnt: Dict[str,int] = Counter()
        for s in seqs:
            kset = probe_kmerify(self.K, s.seq)
            if pool is not None:
                kset &= pool
            cnt.update(kset)            
        return cnt

    def prepare(self, fullepiset: Set[str], pool: Set[str] = None) -> None:
        ''' prepare for off-by calculations by setting up structure '''
        if self.prepared:
            raise RuntimeError("Already prepared!")
        if pool is None:
            pool = fullepiset
        self.oblist = mkoblist(self.K,fullepiset,subepiset=pool,offby=self.offby)
        #for kk in pool:
        #    self.oblist[kk].add(kk) ## include self
        self.prepared = True

    def prepare_fromseqs(self, seqs: List) -> None:
        if self.prepared:
            raise RuntimeError("Already prepared!")
        fullepiset = self.all_kmers_exact(seqs)
        self.prepare(fullepiset)

    def size(self) -> Tuple[int,int]:
        a = len(self.oblist)
        c = sum(len(self.oblist[epi]) for epi in self.oblist)
        return a,c

    def allkmers(self) -> Set[str]:
        ''' returns all the kmers in the database oblist '''
        ## these will be kmers from the "pool" in prepare();
        ## could be a subset of (or different from?) kmers in seqs
        assert( self.prepared )
        return set( self.oblist.keys() )

    @lru_cache(maxsize=None)
    def kmerify_withneighbors(self,s) -> Set[str]:
        ''' 
        for a single sequence s, produce a set of kmers that are in the sequence,
        plus all kmers that are neighbors (ie, that are within offby) of those kmers
        '''
        kmer_set = probe_kmerify(self.K, s.seq) 
        return set(kk for kk in self.oblist.keys()
                   if not self.oblist[kk].isdisjoint(kmer_set))
        
    def all_kmers_wcount(self, seqs: List, pool: Set[str] = None) -> Dict[str,int]:
        ''' for all kmers in the pool, return the coverage of those kmers in the 
        sequence set; ie, the number of sequences that have a kmer that is near one
        of the kmers in the pool.  if pool not specified, then do this for all the kmers
        in the database
        '''
        assert( self.prepared )

        obcnt: Dict[str,int] = Counter()

        if pool is None:
            for s in seqs:
                nbrs = self.kmerify_withneighbors(s)
                obcnt.update( nbrs )
        else:
            pool = set(pool) & set(self.oblist.keys())
            for s in seqs:
                nbrs = self.kmerify_withneighbors(s)
                obcnt.update( nbrs & pool )
            #kmer_set = probe_kmerify(self.K, s.seq) 
            ## Want: for each kk in kpool, are any kmers within offby?
            #obcnt.update(kk for kk in pool 
            #             if not self.oblist[kk].isdisjoint(kmer_set))

            # slow, straightforward approach
            #    for kk in pool:
            #        if any( hamming(kk,x) <= offby for x in kmer_set ):
            #            obcnt[kk] += 1

                
        return obcnt

    def coverage(self, seqs: List, kmers: List) -> int:
        ''' return the number of sequences that are covered (within offby) by at least one kmer in kmers '''
        assert( self.prepared )
        c = 0
        for s in seqs:
            kmer_set = probe_kmerify(self.K, s.seq) 
            for kmer in kmers:
                if not kmer_set.isdisjoint( self.oblist[kmer] ):
                    c += 1
                    break
        return c

    def covered_sequences(self, seqs: List, kmers: List[str], 
                          i_mean_uncovered: bool = False) -> List:
        assert( self.prepared )
        retseqs = []
        for s in seqs:
            kmer_set = probe_kmerify(self.K, s.seq) 
            s_covered=False
            for kmer in kmers:
                if not kmer_set.isdisjoint( self.oblist[kmer] ):
                    s_covered = True
                    break
            if s_covered and not i_mean_uncovered:
                retseqs.append(s)
            if not s_covered and i_mean_uncovered:
                retseqs.append(s)
        return retseqs

    def uncovered_sequences(self, seqs: List, kmers: List[str]) -> List:
        return self.covered_sequences(seqs,kmers,i_mean_uncovered=True)
                
if __name__ == "__main__":

    import readseq

    seqs: List = readseq.read_fasta("Data/pa.fasta",rmdash=True)

    K: int = 31
    offby: int = 2
    kob = kmeroffby(K,offby)
    print("preparing...",end="",flush=True)
    kob.prepare_fromseqs(seqs)
    print("ok")
    print("size of ob:",kob.size())

    kmers = [seqs[0].seq[:K]]
    print(f"Coverage: {kmers[0]} {kob.coverage(seqs,kmers)}/{len(seqs)}")


    

    

        


     



    
       
