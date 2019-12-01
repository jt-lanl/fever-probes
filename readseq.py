#!/usr/bin/env python
## Library of routines for reading sequence files (mase, fasta, tbl, seq)
## Type "seq" is raw sequences (no names)

## read_seqfile() encapsulates all of these functions

import re

FILETYPES = ["mase", "fasta", "tbl", "seq"]

class SequenceSample:
    def __init__(self,name,seq):
        self.name = name
        self.seq = seq

def auto_type(filename):
    suffix_is = lambda s: filename.lower().endswith("."+s)
    for t in FILETYPES:
        if filename.lower().endswith("."+t):
            return t
    return ""

def read_seqfile(filename,type="auto",rmdash=False):
    atype = auto_type(filename)

    if type != "auto" and type not in FILETYPES:
        raise RuntimeError("type="+type+" not supported")
    if type != "auto" and atype and atype != type:
        raise RuntimeError("type="+type+" but file appears of type "+atype)

    if type == "auto":
        type = atype
    if type == "mase":
        seqlist = read_mase(filename)
    elif type == "fasta":
        seqlist = read_fasta(filename)
    elif type == "tbl":
        seqlist = read_tbl(filename)
    elif type == "seq":
        seqlist = read_rawseq(filename)
    else:
        raise RuntimeError("Unknown type ["+type+
                           "] of sequence file ["+filename+"]")

    if rmdash:
        re_dash = re.compile('-')
        for s in seqlist:
            #print("s=",s)
            s.seq = re_dash.sub('',s.seq)

    return seqlist

def read_rawseq(filename):
    seqlist = []
    with open(filename,'r') as f:
        for line in f:
            line = line.strip()
            seqlist.append( SequenceSample('',line) )
    return seqlist

def read_mase(filename):
    re_comment = re.compile('^;.*')
    re_white = re.compile('\s+')

    seq_samples=[]
    cur_name=''
    cur_seq=''
    with open(filename,'r') as mase:
        for line in mase: #.readlines():
            line = line.strip()
            ## This should not be necessary!
            if re_white.match(line):
                raise RuntimeError("Whitespace in MASE file: %s" % (filename,))
            if cur_name:
                if re_comment.match(line):
                    sample = SequenceSample(cur_name,cur_seq)
                    seq_samples.append(sample)
                    cur_name=''
                    cur_seq=''
                else:
                    cur_seq += line
            else:
                line = re_comment.sub('',line)
                if not line:
                    continue
                cur_name=line
    if cur_name:
        ## append the last sequence in the file
        sample = SequenceSample(cur_name,cur_seq)
        seq_samples.append(sample)
    return seq_samples

def read_fasta(filename,rmdash=False):
    re_seqname = re.compile('^>')
    re_white = re.compile('\s+')
    re_dash = re.compile('-')

    seq_samples=[]
    cur_name=''
    cur_seq=''
    with open(filename,'r') as fasta:
        for line in fasta:
            line = line.strip()
            ## This should not be necessary!
            if re_white.match(line):
                raise RuntimeError("Whitespace in FASTA file: %s" % (filename,))
            if cur_name:
                if re_seqname.match(line):
                    sample = SequenceSample(cur_name,cur_seq)
                    seq_samples.append(sample)
                    cur_name = re_seqname.sub('',line)
                    cur_seq=''
                else:
                    if rmdash:
                        line = re_dash.sub('',line)
                    ## convert to uppercase
                    line = line.upper()
                    cur_seq += line
            else:
                if re_seqname.match(line):
                    cur_name = re_seqname.sub('',line)
                    cur_seq = ''
                else:
                    raise RuntimeError("Seq name not yet defined")

    if cur_name:
        ## append the last sequence in the file
        sample = SequenceSample(cur_name,cur_seq)
        seq_samples.append(sample)

    return seq_samples

def read_tbl(filename,rmdash=False):
    '''
    read_tbl: read .tbl file, return list of sequences
    rmdash=False: if True, remove "-" from string; useful for unaligned
    '''
    re_comment    = re.compile('^\#.*')
    re_leadwhite  = re.compile('^\s*')
    re_trailwhite = re.compile('\s*$')
    re_badchar    = re.compile('[\#\$\*X]')
    re_dash       = re.compile('-')

    seq_samples=[]
    if not filename:
        return seq_samples

    with open(filename,'r') as tbl:
        for line in tbl.readlines():
            line = re_comment.sub('',line)
            line = re_leadwhite.sub('',line)
            line = re_trailwhite.sub('',line)
            if not line:
                continue
            tokens = line.split()
            if len(tokens) != 2:
                print("Invalid line in tbl file:")
                print("[",line,"]")
                continue
            (name,seq) = tokens[:2]
            if seq:
                seq = re_badchar.sub('x',seq)
                if rmdash:
                    seq = re_dash.sub('',seq)
                sample = SequenceSample(name,seq)
                seq_samples.append(sample)

    return seq_samples

if __name__ == "__main__":
    
    import argparse
    argparser = argparse.ArgumentParser()
    paa = argparser.add_argument
    paa("filename",
        help="Name of input file")
    paa("--type","-t",default="auto",
        help="Type of input file (mase, fasta, tbl, etc)")
    paa("--rmdash",action="store_true",
        help="Remove dashes from sequences")
    args = argparser.parse_args()
    
    seqlist = read_seqfile(args.filename,type=args.type,rmdash=args.rmdash)
    n = min(5,len(seqlist))
    for s in seqlist[:n]:
        print(s.name,s.seq[:10],"...")
        
