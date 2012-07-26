#! /usr/bin/env python
#coding=utf-8
"""Module Description
Copyright (c) 2012 Jinyan HUANG <jhuang@cephb.fr>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License.
@meqa
@status:  experimental
@version: $Revision$
@author:  Jinyan HUANG
@contact: jhuang@cephb.fr

"""
import sys, string, re, math, os
from optparse import OptionParser
import logging
import time
import subprocess

# ------------------------------------
# constants
# ------------------------------------
logging.basicConfig(level=20,
                    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w"
                    )

# ------------------------------------
# Misc functions
# ------------------------------------
error   = logging.critical
warn    = logging.warning
debug   = logging.debug
info    = logging.info

def prepare_optparser ():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = "usage: %prog <-f bamfilename  -d fileOut>[options]"
    description = "-f is bamfilename, -d is for result output."
    optparser = OptionParser(version="%prog",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-f","--bamfilename",dest="bamfilename",type="string",
                         help="wt file.")
    optparser.add_option("-d","--outdir",dest="outdir",type="string",
                         help="mut file.")                     
    return optparser
def opt_validate (optparser):
    """Validate options from a OptParser object.
    Ret: Validated options object.
    """
    (options,args) = optparser.parse_args() 
    if not options.bamfilename:
        optparser.print_help()
        sys.exit(1)
    if not options.outdir:
        optparser.print_help()
        sys.exit(1)
    return options
def samstat(samstatpath,bamfilename,outdir):
    if os.path.exists(bamfilename):
        cmd="%ssamstat %s" % (samstatpath,bamfilename)
        info(cmd)
        p=subprocess.Popen(['/bin/bash', '-c', cmd])
        sts = os.waitpid(p.pid, 0)
        info ('#Congratulations! See %s.html for the mapping quality control!' %(bamfilename))
    else:
        info("file %s not exist." % (bamfilename))
def main():
    # read the options and validate them
    tstart=time.time()
    options=opt_validate(prepare_optparser())
    info ("Run python script.")
    bamfilename = options.bamfilename
    outdir = options.outdir
    samstat(bamfilename,outdir)
    tend =time.time()
    elapsed = (tend-tstart)/60.0
    info("Total run time is:%s mins",elapsed)
if __name__=="__main__":
    try:
        main()
        info ("Successful run!!!")
    except KeyboardInterrupt:
        warn("User interrupts me! ;-) See you!")
        sys.exit(0)
