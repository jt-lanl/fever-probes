import sys
import re
from collections import Counter
import argparse

import probe
import offby
import readseq

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
    paa("--offby","-o",type=int, default=0,
        help="string matches can be off by this many characters")
    paa("--patt",
        help="only use sequences whose names match this pattern (eg, 'Zika')")
    paa("--rmdup",action="store_true",
        help="remove duplicate sequences")
    paa("--verbose","-v",action="count",
        help="verbosity")
    args = ap.parse_args()
    return args

def best_pairs(kob,seqs,pool=None):
    ''' return a list of the best pairs of kmers '''
    
    ## Get coverage dictionary for kmers
    ## cov[kmer] = #seqs in which kmer appears
    cov = kob.all_kmers_wcount(seqs,pool=pool)

    ## could maybe use a more sophisticated sort key
    ## that also accounted for gc content
    kall = sorted(cov.keys(), key=cov.get, reverse=True)

    kx = None   ## kx=(ki,kj) is best pair so far
    kxlist = [] ## list of kx's
    covx = 0    ## coverage for best pair so far
    for i,ki in enumerate(kall):
        #if i>60:
        #    print("break early!: covi=",covi,cov[ki])
        #    break ### TEMPORARY!!
        covi = cov[ki]
        vprint("i=",i,"#pairs=",len(kxlist),"cov-pair:",covx,"cov-one:",covi)
        if covi*2 < covx:
            break

        ## Get the subset of sequences that are not covered by ki
        subseqs = kob.uncovered_sequences(seqs,[ki])
        vvprint(len(subseqs),"uncovered sequences")
        if len(subseqs)==0:
            print("Full coverage from single sequence:",ki)
            kxlist.append([ki,ki])
            return kxlist

        ## kall[i+1:] is the pool of remaining kmers
        covres = kob.all_kmers_wcount(subseqs,kall[i+1:])
        covj = max(covres.values())
        if covi+covj < covx:
            continue
        if covi+covj > covx:
            covx = covi+covj
            kxlist=[]
        kres = [k for k in covres if covres[k]==covj]
        for kj in kres:
            kx = [ki,kj]
            kxlist.append(kx)
            #vvprint("[%d+%d]" % (covi,covj),end=" ")
            vvprintcoverage(seqs,kx,cov,covres)

    return kxlist

def fastprintcoverage(seqs,kpair,cov,covres):
    ki,kj = kpair
    K = len(ki)
    ca = cov[ki]
    cb = cov[kj]
    c = cov[ki] + covres[kj]
    print("** %3d %.3f %4d [%4d+%4d] %s %s %2d %2d %s %s" %
          (K,c/len(seqs),c,ca,cb,
           #ck.max_distance(seqs,kpair),
           kpair[0],kpair[1],
           probe.gc_content(kpair[0]),
           probe.gc_content(kpair[1]),
           probe.longestStem(kpair[0]),
           probe.longestStem(kpair[1])))
    
def printcoverage_pair(seqs,kpair,offby=0):
    ''' print coverage for a pair of kmers:
    this is the slow version that doesn't take advantage of the offby accelerations;
    but it therefore serves as a check on that faster (but maybe less robust) code
    '''
    K = len(kpair[0])
    ca = probe.coverage_offby(seqs,kpair[:1],offby)
    cb = probe.coverage_offby(seqs,kpair[-1:],offby)
    c =  probe.coverage_offby(seqs,kpair,offby)
    s = "   %3d %.3f %4d [%4d+%4d] %2d %s %s %2d %2d %s %s" % (
           K,c/len(seqs),c,ca,cb,
           probe.max_distance(seqs,kpair),
           kpair[0],kpair[1],
           probe.gc_content(kpair[0]),
           probe.gc_content(kpair[1]),
           probe.longestStem(kpair[0]),
           probe.longestStem(kpair[1]))
    return s

def main(args):

    seqs = readseq.read_seqfile(args.inseqfile,rmdash=True)
    vvprint("Number of Sequences:",len(seqs))

    if args.patt:
        patt = re.compile(args.patt)
        seqs = [s for s in seqs if patt.search(s.name)]
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
    toc("ok {:.4f} sec")
    vprint("original kmers:",len(allkmers))
    vprint("filtered kmers:",len(filtkmers))


    tic("Prepare kmer-offby...")
    kob.prepare(allkmers,pool=filtkmers)
    toc("ok {:.4f} sec")
    vprint("kob.size: %d kmers, %d neigbors" % kob.size())

    tic("Begin timing best_pairs...\n")
    kpairlist = best_pairs(kob,seqs,pool=filtkmers) ## include filtkmers?
    toc("Best pairs: {:.4f} sec")

    ## sort best pairs by GC content
    def pair_mingc_key(kpair):
        ki,kj = kpair
        return -1*sum([probe.gc_content(ki),
                       probe.gc_content(kj)])
        
    kpairlist = sorted(kpairlist,key=pair_mingc_key)

    
    npair = min([len(kpairlist),10])
    print("Found",len(kpairlist),"equivalent pairs; printing",npair,"of them")
    tic()
    for kpair in kpairlist[:npair]:
        s = printcoverage_pair(seqs,kpair,args.offby)
        print(s)
    toc("Print coverage: {:.4f} sec")

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
            fastprintcoverage(*parg,**kwarg)

    tic,toc = tictoc.tictoc_ifverbose(True) #args.verbose)

    main(args)
                
           
                


  
