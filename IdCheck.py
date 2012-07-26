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
import glob
import ConfigParser
import re
import operator
import libIdCheck
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
    usage = "usage: %prog <-f fastq file -o fileOut>[options]"
    description = "-f is fastq file name, -o is directory for result output."
    optparser = OptionParser(version="%prog",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-m","--map",dest="map",type="string",
                         help="map file.")
    optparser.add_option("-p","--ped",dest="ped",type="string",
                         help="ped file.")
    optparser.add_option("-b","--bamfile",dest="bamfile",type="string",
                         help="bam file. If your reads are not mapped, please using -f and -r options.")
    optparser.add_option("-f","--read1",dest="read 1 file.(Forward read.)",type="string",
                         help="read 1 file. If it is single read, please put the single read as read1.")
    optparser.add_option("-r","--read2",dest="read 2 file.(Reverse read.)",type="string",
                         help="read 2 file. If it is single read, please do not use this option.")
    optparser.add_option("-s","--sampleID",dest="sampleID",type="string",
                         help="sample info ID, e.g., demo")
    optparser.add_option("-c","--configfile",dest="configfile",type="string",
                         help="set the configuration file for IdCheck. In this configuration file, you can finer control IdCheck. The configuration file is end with cfg. e.g., IdCheck.cfg")
    optparser.add_option("-o","--output",dest="output",type="string",
                         help="result file output path.")
    return optparser
def opt_validate (optparser):
    """Validate options from a OptParser object.
    Ret: Validated options object.
    """
    (options,args) = optparser.parse_args()
    if not options.map:
        logging.info("Please set the map file for genotype.")
        optparser.print_help()
        sys.exit(1)
    if not options.ped:
        logging.info("Please set the ped file for genotype.")
        optparser.print_help()
        sys.exit(1)
    if not options.sampleID:
        logging.info("Please set an ID for the sample.")
        optparser.print_help()
        sys.exit(1)
    if not options.output:
        logging.info("Please set the output path for the result.")
        optparser.print_help()
        sys.exit(1)
    if not options.bamfile and not options.read1:
        logging.info("At least the bam file or the fastq file(read1) should be set.")
        optparser.print_help()
        sys.exit(1)
    if not options.configfile:
        logging.info("Please set the configuration file.")
        optparser.print_help()
        sys.exit(1)
    return options
def main():
    tstart=time.time()
    options=opt_validate(prepare_optparser())
    info ("We are now running IdCheck tools.")
    idcheck=idcheck(options)
    idcheck.putIdCheckInfo()
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
