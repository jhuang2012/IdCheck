import sys, string, re, math, os
from optparse import OptionParser
import logging
import time
import subprocess
import glob
import ConfigParser
import re
import operator
from qc import qulity
from bwa import bwa_aln,bwa_sampe,bwa_samse
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
class idcheck:
    def __init__(self,options):
        #parameters from command line
        self.mapfile=options.map
        self.bamfile=options.bamfile
        self.fout=options.output.strip()
        if self.fout[-1]!="/":
            self.fout=self.fout+"/"
        self.ped=options.ped
        self.name=options.sampleID
        self.rnaseqSNP=self.fout+self.name+"-rnaseq.map"
        try:
            self.fseqR1=options.f 
        except:
            self.fseqR1=""
        try:
            self.fseqR2=options.r 
        except:
            self.fseqR2=""
        config = ConfigParser.ConfigParser()
        config.read(options.configfile)
        #parameters from configuration file.
        self.genome=config.get('Genome', 'refgenome')
        self.SequenceQA=config.get('pipeline', 'SequenceQA')
        self.Align=config.get('pipeline', 'Align')
        self.MapQA=config.get('pipeline', 'MapQA')
        self.IdCheck=config.get('pipeline', 'IdCheck')
        #parameters for BWA.
        self.bwapath=config.get('toolsPath', 'bwa')
        self.bwaindex=config.get('BWA', 'bwaindex')
        self.bwaaln=config.get('BWA', 'bwaindex')
        self.bwaindex=config.get('BWA', 'bwaindex')
        self.bwasampe=config.get('BWA', 'bwasampe')
        
        self.samstatpath=config.get('toolsPath', 'samstat')
        
        #IdCheck path
        strfilepath = os.path.realpath(__file__) 
        IdCheckPath = "%s/" % (os.path.dirname(strfilepath),)


    def putIdCheckInfo(self):
        if self.SequenceQA.lower()=='y':
            info("Sequencing quality access:")
            info(" Yes.")
            info("Sequence quality control for %s",self.fseqR1)
            qulity(self.fseqR1,self.fout)
            if self.fseqR2!="":
                info("Sequence quality control for %s",self.fseqR2)
                qulity(self.fseqR2,self.fout)
        else:
            info("Sequencing quality access:")
            info(" No.")
        if self.Align.lower()=='y':
            info("Mapping the read to reference genome:")
            info(" Yes.")
            bwa_aln(self.bwapath,self.fseqR1,self.bwaindex,self.bwaaln,self.fout,self.name+"1")
            if self.fseqR2!="":
                bwa_aln(self.bwapath,self.fseqR2,self.bwaindex,self.bwaaln,self.fout,self.name+"2")
                bwa_sampe(self.bwapath,self.bwasampe,self.bwaindex,self.fseqR1,self.fseqR2,self.fout,self.name)
            else:
                bwa_samse(self.bwapath,bwasamse,self.bwaindex,self.fseqR1,self.fout,self.name)
        else:
            info("Mapping the read to reference genome:")
            info(" No.")
        if self.MapQA.lower()=='y':
            info("Doing quality accessing of sequence mapping:")
            info(" Yes.")
            bamfile=self.name+".bam"
            samstat(self.samstatpath,self.fout+bamfile,self.fout)
        else:
            info("Doing quality accessing of sequence mapping:")
            info(" No.")
        if self.IdCheck.lower()=='y':
            info("Sample Identity checking:")
            info(" Yes.")
            self.runIdCheck()
        else:
            info("Sample Identity checking:")
            info(" No.")
    
    def complement(self,s): 
        """Return the complementary sequence string.""" 
        basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
        letters = list(s)
        letters = [basecomplement[base] for base in letters] 
        letters.sort()
        return ''.join(letters)
    def snpCompare(self,af1,af2):
        snp1=[af1[0],af1[1]]
        snp1.sort()
        snp2=[af2[0],af2[1]]
        snp2.sort()
        allele1="".join(snp1)
        allele2="".join(snp2)
        if allele1!=allele2 and allele1!=self.complement(allele2):
            return False
        else:
            return True
    def cleanline(self,line):
        line=line.strip()
        pattern1="\+[0-9]+[ACGTNacgtn]+"
        pattern2="\-[0-9]+[ACGTNacgtn]+"
        pattern3="\^."
        rexp = re.compile(pattern1)
        line1=re.sub(rexp,"",line)
        rexp = re.compile(pattern2)
        line2=re.sub(rexp,"",line1)
        rexp = re.compile(pattern3)
        line3=re.sub(rexp,"",line2)
        return line3
    def findmatch(self,s,mlist):
        cnt={}
        for c in s:
            if c in mlist:
                cnt[c] = cnt.get(c,0) + 1
        if len(cnt)>0:
            return max(cnt.iteritems(), key=operator.itemgetter(1))[0]
        else:
            return False
    def mpileup(self,mapfile,fout,genome,rnaseqSNP,bamfile,name):
        maptmp=fout+name+"maptmp"
        mpileuptmp=fout+name+"mpileup"
        cmd1="cut -f 1,4 %s > %s" % (mapfile,maptmp)
        cmd2="samtools mpileup -C50 -Q25 -q1 -l %s -f %s %s |awk '{if($4>99)print}'> %s" % (maptmp,genome,bamfile,mpileuptmp)
        cmd3="rm %s" % (maptmp)
        for cmd in [cmd1,cmd2,cmd3]:   
            info(cmd)
            p=subprocess.Popen(['/bin/bash', '-c', cmd])
            sts = os.waitpid(p.pid, 0)
        fw=open(rnaseqSNP,"w")
        for line in open(mpileuptmp):
            sp=line.split("\t")
            sp[4]=self.cleanline(sp[4])
            sp[4]=sp[4].upper()
            allele=""
            line="\t".join(sp)
            s=""
            if self.findmatch(sp[4],["A","C","G","T"]):
                s= self.findmatch(sp[4],["A","C","G","T"])
                allele+=s
            if self.findmatch(sp[4],[",",".","$"]):
                if len(allele)==0:
                    allele+=sp[2]+sp[2]+"\t0"
                else:
                    r=1.0*sp[4].count(s)/int(sp[3])
                    cutoff=0.05
                    if r>cutoff and r<(1-cutoff):
                        if r>0.5:
                            r=1-r
                        allele+=sp[2]
                    if r>=(1-cutoff):
                        allele=s+s #Find, but little
                    if r<=cutoff:
                        allele=sp[2]+sp[2]
                    if r>=0.5:
                        r=1-r
                        allele+="\t"+str(r)
                    else:
                        allele+="\t"+str(r)
            else:
                allele+=s+"\t0"
            outstr="%s\t%s\t%s\n" % (sp[0],sp[1],allele.upper())
            fw.write(outstr)
        fw.close()
    def extractSNP(self,mapfile,fout,ped,rnaseqSNP,name):
        tmpsnp=fout+name+"tmpsnp"
        fname=fout+name+"inv.txt"
        fw=open(fname,"w")
        fw.write("\t".join(name.split("_")))
        fw.close()
        cmd1="python %smyjoin.py -F 1,2 -f 1,4 -m y -a %s -b %s |cut -f 6 > %s"\
            % (IdCheckPath,rnaseqSNP,mapfile,tmpsnp)
        compare=fout+name+"-compare"
        
        cmd2="%splink --ped %s --map %s --keep %s --extract %s --tab --recode --out %s" % (IdCheckPath,ped,mapfile,fname,tmpsnp,compare)
        cmd3="rm %s" % (fname)
        cmd4="rm %s" % (tmpsnp)
        for cmd in [cmd1,cmd2,cmd3,cmd4]:   
            info(cmd)
            p=subprocess.Popen(['/bin/bash', '-c', cmd])
            sts = os.waitpid(p.pid, 0)
        if not os.path.exists(compare+".ped"):
            info("The genotype file not exist!")
            sys.exit()
    def compare(self,rnaseq,fout,name,result):
        map=fout+name+"-compare.map"
        ped=fout+name+"-compare.ped"
        tmp=fout+name+"tmp"
        tmp1=fout+name+"tmp1"
        cmd1="cut -f 7- %s | awk -f /env/ceph/home/jhuang/bin/rnaseq/pipeline/tools/transpose1.awk |paste %s - > %s"\
            % (ped,map,tmp)
        cmd2="python %smyjoin.py -F 2 -f 4 -m y -a %s -b %s |sed 's/ //g'> %s"\
            % (IdCheckPath,rnaseq,tmp,tmp1)
        for cmd in [cmd1,cmd2]:   
            info(cmd)
            p=subprocess.Popen(['/bin/bash', '-c', cmd])
            sts = os.waitpid(p.pid, 0)
        matchnum=0
        totalnum=0
        fw=open(result,"w")
        for line in open(tmp1):
            sp=line.split()
            if sp[8].strip()!="00":
                totalnum+=1
                if self.snpCompare(sp[2].strip(),sp[8].strip()):
                    matchnum+=1
                    outstr="%s\t%s\t%s\t%s\t%s\tOK\n" % (sp[0],sp[1],sp[2],sp[8],sp[3])
                else:
                    outstr="%s\t%s\t%s\t%s\t%s\tProblem\n" % (sp[0],sp[1],sp[2],sp[8],sp[3])
                fw.write(outstr)
        outstr="Summary:\n"
        outstr+="total:\t%s\n" % (totalnum)
        outstr+="match:\t%s\n" % (matchnum)
        outstr+="Percentage:\t%s\n" % (1.0*matchnum/totalnum)
        fw.write(outstr)
        fw.close()
        cmd="rm %s %s %s %s" % (map,ped,tmp,tmp1)
        p=subprocess.Popen(['/bin/bash', '-c', cmd])
        sts = os.waitpid(p.pid, 0)
    def runIdCheck(self):
        self.mpileup(self.mapfile,self.fout,self.genome,self.rnaseqSNP,self.bamfile,self.name)
        self.extractSNP(self.mapfile,self.fout,self.ped,self.rnaseqSNP,self.name)
        result=self.fout+self.name+"-res.txt"
        self.compare(self.rnaseqSNP,self.fout,self.name,result)