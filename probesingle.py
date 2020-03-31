from __future__ import division, print_function

import sys
import argparse
from collections import Counter

import readseq
from hamming import hamming
import probe
import offby

import tictoc

def getargs():
    ap = argparse.ArgumentParser()
    paa = ap.add_argument
    paa("inseqfile",
        help="Input sequence file in fasta/mase/tbl format")
    paa("-K",type=int,required=True,
        help="K value for K-mers")
    paa("--gc","-g",type=int, default=0,
        help="minimum number of G/C in K-mer")
    paa("--gcfrac",type=float, default=0,
        help="minimum fraction of G/C in K-mer")
    paa("--stem","-s",type=int, default=0,
        help="maximum length of stem in stemloop structure")
    paa("--offby","-o",type=int,default=0,
        help="string matches can be off by this many characters")
    paa("--patt",
        help="only use sequences whose names match this pattern (eg, 'Zika')")
    paa("--rmdup",action="store_true",
        help="remove duplicate sequences")
    paa("--exact","-x",
        help="input sequence file from which probe must match exactly")
    paa("--eval","-e",
        help="evaluate specified K-mer")
    paa("--unconstrained",action="store_true",
        help="ignore GC content and STEM length constraints")
    paa("--verbose","-v",action="count",
        help="verbosity")
    args = ap.parse_args()
    return args
        
def printcoverage_header(K):
    print(" K  --exact--  --offby-1  --offby-2  kmer" + " "*(K-3) + "GC LL SL stem")
  
def printcoverage(seqs,ksingle):
    kmers=[ksingle]
    cx = probe.coverage_offby(seqs,kmers,offby=0)
    c1 = probe.coverage_offby(seqs,kmers,offby=1)
    co = probe.coverage_offby(seqs,kmers,offby=2)
    stem = probe.longestStem(ksingle)
    print("%2d %4d %.3f %4d %.3f %4d %.3f  %s %2d %2d %2d %s" % 
              (len(ksingle),cx,cx/len(seqs),c1,c1/len(seqs),
               co,co/len(seqs),ksingle,probe.gc_content(ksingle),
               probe.lenStemLoop(ksingle,stem),len(stem),stem))

def common_misfits(seqs,ksingle,offby=1):

    K = len(ksingle)

    cnt = Counter()
    for s in seqs:
        kmer_set = probe.kmerify(K,s.seq)
        snbrs=set()
        for x in kmer_set:
            if offby-1 < hamming(ksingle,x) <= offby:
                snbrs.add(x)
        cnt.update(snbrs)

    z = cnt.most_common(15)
    print("most commmon mismatches")
    print("     %s" % ksingle)
    for p,n in z:
        h = hamming(p,ksingle)
        print("%4d %s [%d]" % (n,p,h))


def main(args):

    seqs = readseq.read_seqfile(args.inseqfile,rmdash=True)
    vprint("Number of Sequences:",len(seqs))
    vprint("Number of Distinct Sequences:",len(set(s.seq for s in seqs)))

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

    seqpool = []
    if args.exact:
        seqpool = readseq.read_seqfile(args.exact,rmdash=True)
        vprint("seqpool:",len(seqpool),"sequences")

    if args.eval:
        printcoverage_header(args.K)
        printcoverage(seqs,args.eval)
        if args.exact:
            printcoverage(seqpool,args.eval)
            print(args.eval in seqpool[0].seq)
            vprint(seqpool[0].seq)
        common_misfits(seqs,args.eval,offby=args.offby)
        return

    tic("Prepare from seqs...")
    kob = offby.kmeroffby(args.K,args.offby)
    kob.prepare_fromseqs(seqs)
    toc("ok {:.4f} secs")

    tic("all_kmers_wcount...\n")
    if args.exact:
        # only include kmers drawn from seqpool
        covpool = kob.all_kmers_exact(seqpool)
        vprint(len(covpool),"pool kmers")
        cov = kob.all_kmers_wcount(seqs,pool=covpool)
        vprint(len(cov),"pool kmers in seqs, min count:",min(cov.values()))
    else:
        cov = kob.all_kmers_wcount(seqs)
    toc("...all_kmers_wcount: {:.4f} sec")

    vprint("cov:",len(cov))
    
    ## all kmers, sorted by coverage (and secondarily by GC content)
    kall = sorted(cov.keys(),key=lambda x: (cov[x],probe.gc_content(x)),reverse=True)
    vprint("kall:",len(kall))
    print("Unconstrained:")
    printcoverage_header(args.K)
    printcoverage(seqs,kall[0])

    if args.unconstrained:
        return ## only look at unconstrained soln

    already_printed = set()
    for stemlen in (args.K,4,3,2):
        kallst = [k for k in kall if not probe.stemfound(k,length=1+stemlen)]
        vprint(len(kall),"--st-->",len(kallst))
        printcoverage_header(args.K)
        for gc in range(GCmin,25):
            ## list gets shorter each time...
            kallst = [k for k in kallst if probe.gc_content(k) >= gc]
            vvprint("gc: ",gc,"-->",len(kallst))
            covcurrent = 0
            for ksingle in kallst:
                if cov[ksingle] < covcurrent:
                    break
                covcurrent = cov[ksingle]
                if ksingle not in already_printed:
                    #print(" %2d  %2d" % (stemlen,gc), end="")
                    #print("[%4d]" % (cov[ksingle],),end="")
                    printcoverage(seqs,ksingle)
                    already_printed.add(ksingle)

                


  
if __name__ == "__main__":

    args = getargs()

    def vprint(*pargs,**kwargs):
        if args.verbose:
            print(*pargs,file=sys.stderr,flush=True,**kwargs)

    def vvprint(*pargs,**kwargs):
        if args.verbose and args.verbose>1:
            print(*pargs,file=sys.stderr,flush=True,**kwargs)

    def donothing(*pargs,**kwargs): pass
    if args.verbose:
        tic = tictoc.tic 
        toc = tictoc.toc
    else:
        tic = toc = donothing

    main(args)

