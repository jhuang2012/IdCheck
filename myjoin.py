#! /usr/bin/python
#coding=utf-8
"""Module Description
Copyright (c) 2009 Jinyan HUANG <jhuang@cephb.fr>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License.

@status:  experimental
@version: $Revision$
@author:  Jinyan HUANG
@contact: jhuang@cephb.fr

"""
import sys, string, re, math, os
from optparse import OptionParser
import logging
import time

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
error   = logging.critical		# function alias
warn    = logging.warning
debug   = logging.debug
info    = logging.info

def prepare_optparser ():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = "usage: %prog <-a wt file -b mut file -f fileOut>[options]"
    description = "-a is wt, -b is mut, -f is for result output."
    optparser = OptionParser(version="%prog",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-a","--fileA",dest="fileA",type="string",
                         help="smaller file.")
    optparser.add_option("-b","--fileB",dest="fileB",type="string",
                         help="bigger file.")
    optparser.add_option("-F","--fileda",dest="fileda",type="string",
                             help="filed for a. Separate by ,")
    optparser.add_option("-f","--filedb",dest="filedb",type="string",
                             help="filed for b, Separate by ,")
    optparser.add_option("-m","--match",dest="match",type="string",
                             help="output only match? y or n, default n")
                        
    return optparser
def opt_validate (optparser):
    """Validate options from a OptParser object.
    Ret: Validated options object.
    """
    (options,args) = optparser.parse_args()
    # if gdb not given, print help, either BED or WIG must be given 
    if not options.fileB:
        optparser.print_help()
        sys.exit(1)
    if not options.fileA:
        optparser.print_help()
        sys.exit(1)
    return options
def reada(fn,filed):
    dict={}
    fd=[]
    if not filed:
        fd.append(1)
    else:
        sp=filed.split(",")
        for i in sp:
            fd.append(int(i))
    for line in open(fn,"r"):
        sp=line.split("\t")
        key=""
        for i in fd:
            key+=sp[i-1]
        key=key.strip()
        dict[key]=line.strip()
    return dict
def output(dict,fn,filed,match):
    fd=[]
    if not filed:
        fd.append(1)
    else:
        sp=filed.split(",")
        for i in sp:
            fd.append(int(i))
    for line in open(fn,"r"):
        line=line.strip()
        sp=line.split("\t")
        key=""
        for i in fd:
            key+=sp[i-1]
        if match=="y":
            if dict.has_key(key):
                print('%s\t%s' % (dict[key],line))
        else:
            if dict.has_key(key):
                print('%s\t%s' % (dict[key],line))
                del dict[key]
    if match!="y":
        for key in dict:
            print('%s' % (dict[key]))
def main():
    # read the options and validate them
    tstart=time.time()
    options=opt_validate(prepare_optparser())
    #info ("Run python script.")
    fileA = options.fileA
    fileB = options.fileB
    fileda=options.fileda
    filedb=options.filedb
    match=options.match
    dict=reada(fileA,fileda)
    output(dict,fileB,filedb,match)
    tend =time.time()
    elapsed = (tend-tstart)/60.0
    #info("Total run time is:%s mins",elapsed)
if __name__=="__main__":
    try:
        main()
        #info ("Successful run!!!")
    except KeyboardInterrupt:
        warn("User interrupts me! ;-) See you!")
        sys.exit(0)
