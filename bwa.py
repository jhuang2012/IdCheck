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
    usage = "usage: %prog <-f wt file -s species -d fileOut>[options]"
    description = "-f is fastq file name, -s is species, -d is for result output."
    optparser = OptionParser(version="%prog",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-f","--fqfilename",dest="fqfilename",type="string",
                         help="wt file.")
    optparser.add_option("-s","--species",dest="species",type="string",
                         help="species.")
    optparser.add_option("-d","--dirOut",dest="dirOut",type="string",
                             help="Result output directory name.")
    return optparser
def opt_validate (optparser):
    """Validate options from a OptParser object.
    Ret: Validated options object.
    """
    (options,args) = optparser.parse_args()
    if not options.fqfilename:
        optparser.print_help()
        sys.exit(1)
    if not options.species:
        optparser.print_help()
        sys.exit(1)
    if not options.dirOut:
            optparser.print_help()
            sys.exit(1)
    return options
def bwa_aln(bwapath,fqfilename,bwaindex,bwaaln,outdir,name):
    if os.path.exists(fqfilename):
        cmd="%sbwa aln %s %s %s > %s/%s.sai" % (bwapath,bwaaln,bwaindex,fqfilename,outdir,name)
        info(cmd)
        p=subprocess.Popen(['/bin/bash', '-c', cmd])
        sts = os.waitpid(p.pid, 0)
        info ('#... cong! Mapping %s is done!' %(name))
    else:
        info(fqfilename)
        return False
def bwa_sampe(bwapath,bwasampe,bwaindex,fqfilename1,fqfilename2,outdir,name):
    if os.path.exists(fqfilename1) and os.path.exists(fqfilename2):
        cmd="%sbwa sampe %s %s %s/%s.sai %s/%s.sai %s %s > %s/%s.sam" % (bwapath,bwasampe,bwaindex,outdir,name+"1",outdir,name+"2",fqfilename1,fqfilename2,outdir,name)   
        info(cmd)
        p=subprocess.Popen(['/bin/bash', '-c', cmd])
        sts = os.waitpid(p.pid, 0)
        info ('#... cong! Mapping bwa_sampe %s is done!' %(name))
    else:
        return False
def bwa_samse(bwapath,bwasamse,bwaindex,fqfilename1,outdir,name):
    if os.path.exists(fqfilename1):
        cmd="%sbwa samse %s %s %s/%s.sai %s > %s/%s.sam" % (bwapath,bwasamse,bwaindex,outdir,name+"1",fqfilename1,outdir,name)   
        info(cmd)
        p=subprocess.Popen(['/bin/bash', '-c', cmd])
        sts = os.waitpid(p.pid, 0)
        info ('#Congratulations! Mapping bwa_sampe %s is done!' %(name))
    else:
        return False
def main():
    # read the options and validate them
    tstart=time.time()
    options=opt_validate(prepare_optparser())
    info ("Run python script.")
    fqfilename = options.fqfilename
    species = options.species
    dirOut = options.dirOut
    bwa_aln(fqfilename,species,outdir)
    bwa_sampe(fqfilename,species,outdir)
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
