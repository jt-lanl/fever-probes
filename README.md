# FEVER Probe Design

FEVER = Fast Evaluation of Viral Emerging Risks

Three main programs are included in this package: `probesingle.py`, `probepair.py`, and `probemulti.py`.  As the names imply, the first seeks a single probe that optimizes coverage, the second seeks a pair of probes that collectively optimize coverage, and the third seeks *n* probes, where *n* is specified by the user.  In terms of operational usage, `probemulti.py` can use *n=1* to obtain a single probe and *n=2* to provide a pair of probes, and this is "all you really need."  But `probesingle.py` and `probepair.py` are useful for research purposes as they provide further analyses.

> `probesingle.py` will provide multiple solutions to the single probe problem, over a range of GC content and stem length constraints, better enabling the user to identify the trade-off between coverage and these quantities.

> `probepair.py` is almost entirely for research purposes.  It seeks the optimal pair of probes by an exhaustive search of all pairs; this is in contrast to the `probemulti.py` code, which uses a more straightforward iterative process.

### Coverage

The probe is a K-mer (string of K bases) which is said to *cover* a sequence if it appears in that sequence.  A cocktail of probes is said to *cover* a sequence if at least one of the probes in the cocktail covers (ie, appears in) the sequence.  Coverage is defined as the fraction of sequences in a database that are covered by the K-mer or cocktail of K-mers.

Note that the database that is used either to optimize coverage or to evaluate coverage is a list of sequences, and that those sequences do *not* need to be aligned.


## PYTHON FILES

All the code in this repository is written in python.
        
        hamming.py
        mkoblist.py
        offby.py
        probe.py
        probemulti.py
        probepair.py
        probesingle.py
        readseq.py
        tictoc.py

## SUPPORT LIBRARIES

readseq: reads fasta/mase/tbl and raw sequence files

offby: class encapsulates mkoblist
    imports mkoblist, probe, hamming

mkoblist: creates data structure to speed up off-by kmer matching
    imports hamming

probe: various utilities for kmers, 
    imports hamming

tictoc: this is imported by some code just to help benchmark performance; it is not logically necessary.

## MAIN EXECUTABLES

probesingle:
    imports readseq, offby, probe, hamming

probepair:
    imports readseq, offby, probe

probemulti: 
    imports readseq, offby, probe, hamming

## USAGE

The main executables use the `argparse` library for providing parameters on the command line. One feature of this package is that you can get a brief list of all the options by using `-h` or `--help` on the command line:

    python probesingle.py -h

which yields:

    usage: probesingle.py [-h] -K K [--gc GC] [--gcfrac GCFRAC] [--stem STEM]
                          [--offby OFFBY] [--patt PATT] [--rmdup] [--exact EXACT]
                          [--eval EVAL] [--unconstrained] [--verbose]
                          inseqfile
    
    positional arguments:
      inseqfile             Input sequence file in fasta/mase/tbl format
    
    optional arguments:
      -h, --help            show this help message and exit
      -K K                  K value for K-mers
      --gc GC, -g GC        minimum number of G/C in K-mer
      --gcfrac GCFRAC       minimum fraction of G/C in K-mer
      --stem STEM, -s STEM  maximum length of stem in stemloop structure
      --offby OFFBY, -o OFFBY
                            string matches can be off by this many characters
      --patt PATT           only use sequences whose names match this pattern (eg,
                            'Zika')
      --rmdup               remove duplicate sequences
      --exact EXACT, -x EXACT
                            input sequence file from which probe must match
                            exactly
      --eval EVAL, -e EVAL  evaluate specified K-mer
      --unconstrained       ignore GC content and STEM length constraints
      --verbose, -v         verbosity
    
    

## SIMPLE EXAMPLE:

    python probesingle.py -K31 Data/pa.fasta

produces many candidate probes (with varying GC and hairpin stem length values) that optimize exact coverage.  There are many options; eg:

    python probesingle.py -K31 Dat/pa.fasta --exact Data/HongKong68.fasta

restricts candidate probes to those that appear (exactly) in the HongKong sequence


### MORE DETAILS ON OPTIONS:

#### I: `probesingle.py`, `probepair.py`, and `probemulti.py` all take the following command-line arguments; note that the first two are required, and the rest are optional

  `inseqfile`             Input sequence file in fasta/mase/tbl format

  The user must always specify an input sequence file.

  `-K K`                  K value for K-mers

  Although our research almost always uses K=31, this is not the default; in fact there is no default -- the user is required to specify a K value.

  `--gc GC, -g GC`        minimum number of G/C in K-mer
  
  `--gcfrac GCFRAC`       minimum fraction of G/C in K-mer

  Default values are zero, which means that by default there are no G/C constraints.  GC is an integer (with K=31, we typically used GC=15 as a minimum), and GCFRAC is a fraction (which is less commonly specified, but is useful when doing a range of K values).  If both are specified, then the maximum of GC and K*GCFRAC is used.

  `--stem STEM, -s STEM`  maximum length of stem in stemloop structure

  There is some confusing terminology here.  The overall probe is a stem loop, and it is constructed by appending stems to the basic probe design that this code produces.  However, the STEM here refers only to the basic probe, and corresponds to the length of a short subsequence that appears along with its reverse complement in the probe.  If the subsequence is long, then it has a propensity to bind with its reverse complement and produce a hairpin in the probe.  This is undesirable, so we seek probe designs with small STEM values.
  Note the specifying STEM=0 implies that this hairpin propensity condition is not considered at all.  Default value is STEM=0, but we typically use STEM=3.

  `--offby OFFBY, -o OFFBY`
                        string matches can be off by this many characters

  A probe is said to "cover" a sequence if that probe appears as a subsequence of that sequence.  This is strictly the case when OFFBY=0, which is its default value.  But if OFFBY=1, for instance, then if a sequence includes a subsequence that nearly matches the probe (that is off by one character), then that sequence will be covered by that probe.
  Note that exact comparison is a much faster operation than approximate comparison, and so setting OFFBY>0 leads to slower runtimes.

  `--patt PATT`           only use sequences whose names match this pattern (eg, 'Zika')

  This option is useful for selecting a subset of the sequences in the inseqfile. It considers only sequences whose names match the specified pattern.  This pattern can be a python-style regular expression, and the pattern can appear anywhere in the name.

  `--rmdup`               remove duplicate sequences

  If two sequences in the `inseqfile`, presumably with different names, have identical sequences, then remove one of them from consideration.  

  `--verbose, -v`         verbosity

  When specified, the program provides more intermediate information.  You can specify levels of verbosity either by multiple invocations (eg, "-v -v" is verbose=2).

#### II: Options available on `probesingle.py`:

  `--exact EXACT, -x EXACT`
                        input sequence file from which probe must match
                        exactly
			
  This option allows the user to provide a second inseqfile, in this case expected to be a file with only one sequences, and demands that the identified probe covers this sequence exactly. For example, using "--exact Data/HongKong68.fasta" demands that the probe be drawn from a pool that includes K-mers that appear in the original Hong Kong flu sequence.					

  `--eval EVAL, -e EVAL`  evaluate specified K-mer

  This is a handy utility that enables the user to evaluate a candidate probe by just including the probe on the command line (eg, "-e CTGAGGCTGAGAAGCAACTCCAACAATATGC")

  `--unconstrained`       ignore GC content and STEM length constraints

  In its usual operation, probesingle.py tries to supply a range of solutions that satisfy a range of GC and STEM constraints; with this option set, it only provides an unconstrained solution.

#### III: Options available for `probemulti.py`:

It is worth remarking that although probesingle.py is guaranteed to find the optimal probe, in the sense of maximizing coverage subject to GC content and STEM length constraints, no such guarantee is available for the probemulti.py code.  The default algorithm is straightforward and greedy: the first probe is chosen to optimize coverage over the full set of sequences; the second probe optimizes coverage over the sequences not covered by the first probe; the third probe optimizes coverage over sequences not covered by the first two probes; and so on.  Two deviations (improvements?) are provided with the '--tweak' and '--random' options.  These are more for research than for practical use, however, since we have observed with most of our influenza data that these deviation rarely provide actual improvement.  Still, they might, particularly on other datasets.

  `-n N`                  Number of probes

  User specifies how many probes are desired.  These are probes that collectively try to optimize coverage.  In the multi-probe scenario, a sequence is said to be covered if any one of the probes covers it.

  `--fix FIX [FIX ...]`   include (fix == don't alter) these probes

  User can specify specific probes to be included.  For instance if n=4 and two probes are specified with the FIX option, then the code will try to optimize coverage by finding two further probes such that the four of them collectively cover as many sequences as possible.

  `--tweak`               tweak solution (sometimes achieves better coverage)

  The default algorithm greedily picks N probes, one at a time, each one optimizing the coverage of the sequences not covered by the previous probes.  When this option is invoked, a second iteration is made through the probes.

  `--random, -r`          use random initial conditions

  This provides another deviation from the default algorithm.  Here, N probes are assigned at random, and then a tweak-style iteration is made through the data.  Since this is random, it can be run multiple times, and can get different answers each time.  This enables a hopeful user to make many runs, and then keep the best result.  (In practice, this best result is often no better than the result obtained by the default variant of the algorithm.)

  `--uncovered`           print out a list of the uncovered sequences

  It can be handy, on output, to list which sequences have not yet been covered by the probes.  This does not affect the running of the algorithm, just what is reported on output.

#### IV: Options available only with `probepair.py`:

Even though `probepair.py` is ostensibly the most experimental of the three provided codes, it does not have any options that are not also available on `probesingle.py`.

# COPYRIGHT

(c) 2021. Triad National Security, LLC. All rights reserved.

This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

# LICENSE

This program is open source under the BSD-3 License.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# LANL C-number

Notice of Open Source Copyright Assertion for C20108 FEVER Probe Design Algorithm


