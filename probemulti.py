import sys

import argparse
from collections import Counter
import random

from hamming import hamming
import probe
import offby
import readseq

import tictoc

def getargs():
    ap = argparse.ArgumentParser()
    paa = ap.add_argument
    paa("inseqfile",
        help="Input sequence file in fasta/mase/tbl format")
    paa("-n",type=int,default=1,
        help="Number of probes")
    paa("-K",type=int,required=True,
        help="K value for K-mers")
    paa("--gc","-g",type=int, default=0,
        help="minimum number of G/C in K-mer")
    paa("--gcfrac",type=float, default=0,
        help="minimum fraction of G/C in K-mer")
    paa("--stem","-s",type=int, default=0,
        help="maximum length of stem in stemloop structure")
    paa("--offby","-o",type=int, default=0,
        help="string matches can be off by this many characters")
    paa("--patt",
        help="only use sequences whose names match this pattern (eg, 'Zika')")
    paa("--rmdup",action="store_true",
        help="remove duplicate sequences")
    paa("--fix",nargs='+',
        help="include (fix == don't alter) these probes")
    paa("--tweak",action="store_true",
        help="tweak solution (sometimes achieves better coverage)")
    paa("--random","-r",action="store_true",
        help="use random initial conditions")
    paa("--uncovered",action="store_true",
        help="print out a list of the uncovered sequences")
    paa("--verbose","-v",action="count",
        help="verbosity")
    args = ap.parse_args()
    return args

def greedy_n_probes(kob, N, seqs, initkmers=[]):
    ''' returns a list of N K-mer probes, which collectively cover many sequences'''

    kmer_list = initkmers.copy()

    for i in range(N):
        while len(kmer_list) > N-1:
            kmer_list.pop(0)
        vvprint(i,kmer_list)
        tic(f"i={i}: uncovered sequences...")
        subseqs = kob.uncovered_sequences(seqs,kmer_list)
        toc("ok {:.2e} sec")
        vvprint("subseq:",len(subseqs))
        if len(subseqs)==0:
            break
        tic(f"all_kmers_wcount_frompool_offby [seq={len(subseqs)}]...")
        if 0:
            pool = kob.all_kmers_exact(subseqs)
            vprint(f"[pool={len(pool)}]...",end="")
        else:
            pool = None        
        covx = kob.all_kmers_wcount(subseqs,pool)
        toc("ok {:.4f} sec")

        def keyfcn(kmer):
            # sort criterion is coverage first, GC content second
            cv = covx[kmer]
            gc = probe.gc_content(kmer)
            return (cv,gc)
            
        knext = max(covx,key=keyfcn)
        vprint(i,covx[knext],knext)
        kmer_list.append(knext)

    return kmer_list


def printcoverage(seqs,kmerlist,offby=0):
    K = len(kmerlist[0])
    cab = "+".join(["%4d" % probe.coverage_offby(seqs,[k],offby) for k in kmerlist])
    kmerstr = " ".join(kmerlist)
    gcstr = " ".join(["%2d" % probe.gc_content(k) for k in kmerlist])
    slenstr = " ".join([probe.longestStem(k) for k in kmerlist])
    c = probe.coverage_offby(seqs,kmerlist,offby)
    print("%3d %.3f %4d [%s] %2d %s %s %s" %
          (K,c/len(seqs),c,cab,
           probe.max_distance(seqs,kmerlist),
           kmerstr,gcstr,slenstr))

def main(args):

    seqs = readseq.read_seqfile(args.inseqfile,rmdash=True)
    vvprint("Number of Sequences:",len(seqs))

    if args.patt:
        seqs = probe.filter_seq_names_by_pattern(seqs,args.patt)
        if len(seqs) == 0:
            raise RuntimeError(f"Pattern /{args.patt}/ not found in any names")
        vprint("Number of Sequences:",len(seqs))

    if args.rmdup: ## remove duplicate sequences
        seqstrset = set() ## set of sequence strings
        seqlist = list()  ## set of sequences (name,str)
        for s in seqs:
            if s.seq not in seqstrset:
                seqlist.append(s)
            seqstrset.add(s.seq)
        seqs = seqlist
        vprint("Number of Sequences:",len(seqs))

    GCmin = max(args.gc, int(0.99+args.K*args.gcfrac))
    vprint("min GC:",GCmin)

    tic("Prepare to prepare...")
    kob = offby.kmeroffby(args.K,args.offby)
    allkmers = kob.all_kmers_exact(seqs)
    filtkmers = probe.filter_kmerset(allkmers.copy(),gc=GCmin,stem=args.stem)
    if args.fix:
        filtkmers.update( args.fix )
    toc("ok {:.4f} sec")
    vprint("original kmers:",len(allkmers))
    vprint("filtered kmers:",len(filtkmers))


    tic("Prepare kmer-offby...")
    kob.prepare(allkmers,pool=filtkmers)
    toc("ok {:.4f} sec")
    vprint("kob.size: %d kmers, %d neigbors" % kob.size())

    fix_kmers=[]
    if args.fix:
        fix_kmers.extend( args.fix )
        N = args.n - len(fix_kmers)
        xseqs = kob.uncovered_sequences(seqs, fix_kmers)
        if len(xseqs)==0:
            print("Full coverage from fixed sequences!")
        vvprint("xseqs:",len(xseqs))
    else:
        N = args.n
        xseqs = seqs
    
    init_kmers=[]
    if args.random:
        ## Initialize with random kmers (not quite so greedy...)
        rnd_kmers = random.sample(filtkmers,N)
        init_kmers.extend( rnd_kmers )
    

    kmer_list = greedy_n_probes(kob,N,xseqs,initkmers=init_kmers)

    if args.tweak:
        ## OOPS, Tweak breaks it!
        print("Before tweak:")
        printcoverage(seqs,fix_kmers + kmer_list,offby=args.offby)

        kmer_list = greedy_n_probes(kob,N,xseqs,initkmers=kmer_list)

        print("After tweak:")
        printcoverage(seqs,fix_kmers + kmer_list,offby=args.offby)

    full_kmer_list = fix_kmers + kmer_list
    tic("Print coverage...\n")
    ## Print coverage for various offby's
    print(" n   -exact--   -offby-1   -offby-2  kmer                            GC LL SL stem")
    for i,kmer in enumerate(full_kmer_list):
        kmers = full_kmer_list[:i+1]

        cx = probe.coverage_offby(seqs,kmers,offby=0)
        c1 = probe.coverage_offby(seqs,kmers,offby=1)
        co = probe.coverage_offby(seqs,kmers,offby=2)
        stem = probe.longestStem(kmer)
        print("%2d %4d %.3f %4d %.3f %4d %.3f  %s %2d %2d %2d %s" % 
              (len(kmers),cx,cx/len(seqs),c1,c1/len(seqs),
               co,co/len(seqs),kmer,probe.gc_content(kmer),
               probe.lenStemLoop(kmer,stem),len(stem),stem))

        if args.verbose and args.uncovered:
            print("Covered by this probe:")
            for s in kob.covered_sequences(seqs,[kmer]):
                print("+",s.name)
            print("Uncovered sequences so far:")
            for s in kob.uncovered_sequences(seqs,kmers):
                print("-",s.name)

            


    ## Ok, which sequences were not covered?
    if args.uncovered:
        print("Uncovered sequences:")
        for s in kob.uncovered_sequences(seqs,full_kmer_list):
            print("",s.name)
        
    toc("coverage: {:.4f} sec")

    try:
        ## Diagnostics for cache'ing (aka memo-izing) functions
        vprint("hamming:",hamming.cache_info())
        vprint("kmerify:",probe.kmerify.cache_info())
        vprint("withneighbors",kob.kmerify_withneighbors.cache_info())
    except:
        pass



if __name__ == "__main__":

    args = getargs()

    def vprint(*parg,**kwarg):
        if args.verbose:
            print(*parg,file=sys.stderr,flush=True,**kwarg)

    def vvprint(*parg,**kwarg):
        if args.verbose and args.verbose>1:
            print(*parg,file=sys.stderr,flush=True,**kwarg)

    def vvprintcoverage(*parg,**kwarg):
        if args.verbose and args.verbose>1:
            printcoverage_pair(*parg,**kwarg)

    def donothing(*pargs,**kwargs): pass
    if args.verbose:
        tic = tictoc.tic 
        toc = tictoc.toc
    else:
        tic = toc = donothing



    main(args)
                


  
