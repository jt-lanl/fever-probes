## Usage
## from toctoc import tic,toc
## In text of code
##    tic()
##    do stuff
##    print(toc(),"seconds so far")
##    do more stuff
##    print(toc(),"seconds total")
## Or:
##    tic("Begin stuff...")
##    do stuff
##    toc("ok: {} sec")
## Or: 
##    tic("Begin...")
##    do stuff
##    print("ok",toc(),"seconds",file=sys.stderr)
##
## Useful alternative, to set up in __main__:
##    import tictoc
##    tic,toc = tictoc.tictoc_ifverbose(args.verbose)
##


import sys
import time

class tictoc:
    def __init__(self):
        self.start_time=None
    def tic(self):
        self.start_time = time.time()
        return self.start_time
    def toc(self):
        return time.time() - self.start_time

T = tictoc() ## a single global clock

def tic(msg=None):
    if msg is not None:
        print(msg,file=sys.stderr,flush=True,end="")
    T.tic()

def toc(msg=None):
    ttoc = T.toc()
    if msg is not None:
        msg = msg.format(ttoc)
        print(msg,file=sys.stderr,flush=True)
    else:
        print("ok {} sec".format(ttoc))    
    return ttoc

def tictoc_ifverbose(verbose):
    '''
    returns functions that can be used as tic
    define returned functions to be tic,toc if verbose is True
    otherwise, return return functions tic,toc if verbose is True
    '''
    if verbose:
        tick = tic 
        tock = toc
    else:
        def donothing(*pargs,**kwargs): pass
        tick = tock = donothing
    return tick,tock
            


if __name__ == "__main__":

    tic("Begin...")
    for i in range(10000):
        xi = sum([k*k for k in range(i)])
    toc()

    theclock = tictoc()
    theclock.tic()
    print("A differnt clock begins...",flush=True)
    for i in range(10000):
        xi = sum([k*k for k in range(i)])
    tic("Default clock...")
    for i in range(10000):
        xi = sum([k*k for k in range(i)])
    toc()
    for i in range(10000):
        xi = sum([k*k for k in range(i)])
    ttoc = theclock.toc()
    print("Different clock elapsed:",ttoc,"sec")




