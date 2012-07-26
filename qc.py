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
import itertools, warnings
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
    usage = "usage: %prog <-f fastq file -o fileOut>[options]"
    description = "-f is fastq file name, -o is directory for result output."
    optparser = OptionParser(version="%prog",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-f","--fa",dest="fafilename",type="string",
                         help="fastq file name.")
    optparser.add_option("-o","--DirOut",dest="DirOut",type="string",
                             help="Directory for result output.")
                        
    return optparser
def opt_validate (optparser):
    """Validate options from a OptParser object.
    Ret: Validated options object.
    """
    (options,args) = optparser.parse_args()
    if not options.fafilename:
        optparser.print_help()
        sys.exit(1)
    if not options.DirOut:
        optparser.print_help()
        sys.exit(1)
    return options
#########################
## Utils
#########################


class FileOrSequence( object ):
   """ The construcutor takes one argument, which may either be a string,
   which is interpreted as a file name (possibly with path), or a
   connection, by which we mean a text file opened for reading, or 
   any other object that can provide an iterator over strings 
   (lines of the file).
   
   The advantage of passing a file name instead of an already opened file
   is that if an iterator is requested several times, the file will be
   re-opened each time. If the file is already open, its lines can be read
   only once, and then, the iterator stays exhausted.      

   From HTSeq
   """
   
   def __init__( self, filename_or_sequence ):      
      self.fos = filename_or_sequence
      self.line_no = None
      
   def __iter__( self ):
      self.line_no = 1
      if isinstance( self.fos, str ):
         if self.fos.lower().endswith( ( ".gz" , ".gzip" ) ):
            lines = gzip.open( self.fos )
         else:
            lines = open( self.fos )
      else:
         lines = self.fos
      for line in lines:
         yield line
         self.line_no += 1
      if isinstance( self.fos, str ):
         lines.close()
      self.line_no = None
         
   def __repr__( self ):
      if isinstance( self.fos, str ):
         return "<%s object, connected to file name '%s'>" % (
            self.__class__.__name__, self.fos )
      else:
         return "<%s object, connected to %s >" % (
            self.__class__.__name__, repr( self.fos ) )   
            
   def get_line_number_string( self ):
      if self.line_no is None:
         if isinstance( self.fos, str ):
            return "file %s closed" % self.fos
         else:
            return "file closed"
      if isinstance( self.fos, str ):
         return "line %d of file %s" % ( self.line_no, self.fos )
      else:
         return "line %d" % self.line_no

class FastqReader( FileOrSequence ):
   """A Fastq object is associated with a FASTQ self.file. When an iterator
   is requested from the object, the FASTQ file is read.
   
   qual_scale is one of "phred", "solexa", "solexa-old".
   
   From HTSeq
   """

   def __init__( self, file_, qual_scale = "phred" ):
      FileOrSequence.__init__( self, file_ )
      self.qual_scale = qual_scale
      if qual_scale not in ( "phred", "solexa", "solexa-old" ):
         raise ValueError, "Illegal quality scale."
      
   def __iter__( self ):
      fin = FileOrSequence.__iter__( self )
      while True:
         id1  = fin.next()
         seq  = fin.next()
         id2  = fin.next()
         qual = fin.next()
         if qual == "":
            if id1 != "":
               warnings.warn( "Number of lines in FASTQ file is not "
                  "a multiple of 4. Discarding the last, "
                  "incomplete record" )
            break
            
         if not qual.endswith( "\n" ):
            qual += "\n"
         if not id1.startswith( "@" ):
            raise ValueError( "Primary ID line in FASTQ file does"
               "not start with '@'. Either this is not FASTQ data or the parser got out of sync." )
         if not id2.startswith( "+" ):
            raise ValueError( "Secondary ID line in FASTQ file does"
               "not start with '+'. Maybe got out of sync." )
         if len( id2 ) > 2 and id1[1:] != id2[1:]:
            raise ValueError( "Primary and secondary ID line in FASTQ"
               "disagree." ) 
         yield [ seq[:-1], qual[:-1] ]
         
def qulity(fafilename,outdir):
    if os.path.exists(fafilename):
        readfile = FastqReader( fafilename )
        readlen = 0
        reads = readfile
        qmax=0
        for r in itertools.islice( reads, 2000 ):
            if r[0]>readlen:
                readlen=len(r[0])
            for i in r[1]:
                if ord(i)>qmax:
                    qmax=ord(i)
        max_qual = 40
        offset=qmax-max_qual
        baseA=[0 for x in range(readlen)]
        baseC=[0 for x in range(readlen)]
        baseG=[0 for x in range(readlen)]
        baseT=[0 for x in range(readlen)]
        baseN=[0 for x in range(readlen)]
        qual=[[] for x in range(readlen)]
        rn=0
        fqual=open(outdir+"/qual.txt","w")
        for r in reads:
            rn+=1
            m=0
            for i in r[0]:
                if i=="A":
                    baseA[m]+=1
                if i=="C":
                    baseC[m]+=1
                if i=="G":
                    baseG[m]+=1
                if i=="T":
                    baseT[m]+=1
                if i=="N":
                    baseN[m]+=1
                m+=1
            m=0
            for i in r[1]:
                #qual[m].append(ord(i)-64)
                s=ord(i)-offset
                fqual.write("%s\t" % s)
                m+=1
            fqual.write("\n")
            if rn>1000000:
                break
        fqual.close()
        tmpfname=os.path.basename(fafilename).split('.')[0]
        fw=open(outdir+"/"+tmpfname+"_qa.R","w")
        fw.write("pdf(\"%s\")\n" % (outdir+"/"+tmpfname+"_qa.pdf"))
        fw.write("a=c%s\n" % str(tuple(baseA)))
        fw.write("c=c%s\n" % str(tuple(baseC)))
        fw.write("g=c%s\n" % str(tuple(baseG)))
        fw.write("t=c%s\n" % str(tuple(baseT)))
        fw.write("n=c%s\n" % str(tuple(baseN)))
        fw.write("col=c(\"red\",\"blue\",\"orange\",\"black\",\"grey\")\n")
        fw.write("plot(seq(1:%s),100*a/%s,ylim=c(0,100),xlim=c(1,%s),ann=FALSE,axes=FALSE,type=\"l\",col=col[1])\n" % (readlen,rn,readlen))
        fw.write("par(new=TRUE)\n")
        fw.write("plot(seq(1:%s),100*c/%s,ylim=c(0,100),xlim=c(1,%s),ann=FALSE,axes=FALSE,type=\"l\",col=col[2])\n" % (readlen,rn,readlen))
        fw.write("par(new=TRUE)\n")
        fw.write("plot(seq(1:%s),100*g/%s,ylim=c(0,100),xlim=c(1,%s),ann=FALSE,axes=FALSE,type=\"l\",col=col[3])\n" % (readlen,rn,readlen))
        fw.write("par(new=TRUE)\n")
        fw.write("plot(seq(1:%s),100*t/%s,ylim=c(0,100),xlim=c(1,%s),ann=FALSE,axes=FALSE,type=\"l\",col=col[4])\n" % (readlen,rn,readlen))
        fw.write("par(new=TRUE)\n")
        fw.write("plot(seq(1:%s),100*n/%s,ylim=c(0,100),xlim=c(1,%s),xlab=\"Position in read(bp)\",ylab=\"Base percentage\",type=\"l\",col=col[5],main=\"Per base sequence content\")\n" % (readlen,rn,readlen))
        fw.write("legend(\"topright\",c(\"A\",\"C\",\"G\",\"T\",\"N\"),cex=1,text.col=col,bty=\"n\",horiz=\"TRUE\")\n")
        
        fw.write("phred<-matrix(scan(\"%s/qual.txt\",what=0),ncol=%s,byrow=TRUE)\n" % (outdir,readlen))
        
        fw.write("m<-apply(phred, 2,median)\n")
        fw.write("names(phred)<-seq(1:%s)\n" % readlen)
        fw.write("col<-rep(\"blue\",%s)\n" % readlen)
        fw.write("col[m<20]<-\"red\"\n")
        fw.write("col[m<10]<-\"grey\"\n")
        fw.write("boxplot(phred,outline=FALSE,xlim=c(0,%s),ylim=c(1,40),col=col,main=\"Per base sequence quality\",xlab=\"Position in read(bp)\",ylab=\"Quality Score\")\n" % readlen)
        
        fw.write("plot(density(rowMeans(phred)),xlab=\"Quality Score\",ylab=\"Density\",main=\"Sequence quality scores density\",col=\"blue\")\n")
                
        fw.write("plot(seq(1:%s),100*(g+c)/%s,ylim=c(0,100),xlim=c(1,%s),xlab=\"Position in read(bp)\",ylab=\"Base(G+C) percentage\",type=\"l\",col=col[2],main=\"Per base GC content\")\n" % (readlen,rn,readlen))
                        
        fw.write("plot(density(100*(g+c)/%s),xlim=c(35,55),xlab=\"GC Percentage\",ylab=\"Density\",type=\"l\",col=col[2],main=\"Per sequence GC density\")\n" % (rn))
        
        fw.write("plot(seq(1:%s),100*n/%s,ylim=c(0,100),xlim=c(1,%s),xlab=\"Position in read(bp)\",ylab=\"Base N percentage\",type=\"l\",col=col[2],main=\"Per base N content\")\n" % (readlen,rn,readlen))
        
        fw.write("dev.off()")
        fw.close()
        
        cmd="R CMD BATCH %s/%s" % (outdir,tmpfname+"_qa.R")
        info(cmd)
        p=subprocess.Popen(['/bin/bash', '-c', cmd])
        sts = os.waitpid(p.pid, 0)
        p=subprocess.Popen(['/bin/bash', '-c', "rm %s/%s" % (outdir,"qual.txt")])
        info ('#Congratulations! See %s for the sequence quality control!' %("qa.pdf"))
    else:
        sys.exit(1)

def main():
    # read the options and validate them
    tstart=time.time()
    options=opt_validate(prepare_optparser())
    info ("Run python script.")
    fafilename = options.fafilename
    outdir = options.DirOut
    qulity(fafilename,outdir)
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
