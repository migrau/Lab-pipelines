#############################################################################
#                                                                           #
# Description                                                               #
#                                                                           #    
#############################################################################

#INPUT: 
#(1) Fasta files downloaded from GISAID (to be downsampled) and (2) Fasta file downloaded from GISAID (to don't be downsampled) (Australia):
#-IDs MUST be like `>A/Mali/31/2019__|_EPI_ISL_377868_|_2019-01-31` (default)
#-file name MUST contain the region. E.g. "h1n1pm_Europe.fa"
#OUTPUT: GTR and treetime trees. Metadata ready to be loaded in Figtree.

#STEPS:
#1. Add the region to the ids in the fasta.
#2. Select X samples per month per fasta.
#3. Merge and align.
#4. Create metadata
#5. RAxML tree
#6. Treetime tree.
#*Manually, Figtree/Illustrator edition.

import subprocess, sys, os, glob 
from os.path import join
from os.path import basename
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
import random
from datetime import datetime as dt
import time

# Input parameters  ------------------------------------------------------------------------
#

#folder with the gisaid fasta files (to be downsamples).
IFA = config["ifa"]
#Output folder
workspace= config["out"]
#samples NOT downsampled (e.g. Australia)
ausfasta= config["ausfa"]
#Downsampling number. E.g. n=5 -> Only 5 samples per month per fasta.
nsamplesPerMonth=config["nsamples"]

#threads for raxml and mafft (alignment)
threads=4
#RAxML bootstrap
bs=20

#mafft, raxml, treetime should exist in $PATH

## Functions -------------------------------------------------------------------
#

def multi2singleFasta(content):
    block=[]
    res=""
    dataList= content.split("\n")
    for line in dataList:
        if line.startswith('>'):
            if block:
                res+=''.join(block) + '\n'
                block = []
            res+=line+ '\n'
        else:
            block.append(line.strip())

    if block:
        res+=''.join(block) + '\n'
    return res

#date to fraction number
def toYearFraction(date):
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction

#Given a fasta, it selects X samples per each month
def samplesXmonth(fafile,n):
    ids_per_month={}

    #Another way could be: create a dict. (https://www.biostars.org/p/314630/)
    #record_dict = SeqIO.to_dict(SeqIO.parse('seqs.fa', 'fasta'))
    #Create the auxiliary ids_per_month dict but use record_dict for the output

    #group ids per month
    for index,record in enumerate(SeqIO.parse(fafile, 'fasta')):
        linesp=record.id.split("_|_")
        if linesp[-1].count("-")==2:
            #take month
            if linesp[-1].split("-")[0]+"_"+linesp[-1].split("-")[1] not in ids_per_month.keys():
                ids_per_month[linesp[-1].split("-")[0]+"_"+linesp[-1].split("-")[1]]=[]
            ids_per_month[linesp[-1].split("-")[0]+"_"+linesp[-1].split("-")[1]].append(record)
        else:
            print("WARNING. Record discarded because incomplete date: "+ record.id)

    res=""
    #take n random per month
    for month in ids_per_month.keys():
        print("Samples in month: "+month+": "+str(len(ids_per_month[month])))
        if len(ids_per_month[month])<n:
            sampling=ids_per_month[month]
        else:
            sampling = random.sample(ids_per_month[month],int(n))
        print("Selected "+ str(len(sampling))+"samples:")
        for ke in sampling:
            #no duplicates
            if ke.id not in res:
                print(ke.id)
                res+=">"+ke.id+"\n"+ke.seq+"\n"
        print("_-_-_-_-_-_-_")
    return res

#Replace all the elements from a list, present in a string
def replace_all(text, mlist):
    for i in mlist:
        text = text.replace(i, "_")
    return text

#Add location (extracted from the file-name) and remove special characters and duplicated entries
def addRegionNoDup(fafile):
    res=""
    myrec=[]
    region=fafile.split(".")[-2].split("_")[-1]
    mlist=["\t", " ", ":", ",", ")", "(", ";", "]", "[", "'"]
    for index,record in enumerate(SeqIO.parse(fafile, 'fasta')):
        line=replace_all(record.id, mlist)
        line=line+"|"+region+"\n"
        if line not in myrec:
            myrec.append(line)
            res+=">"+line+record.seq+"\n"
    return res

#Extract metada, ready to load in Treetime and Figtree.
def extractMeta(fafile):
    datescsv="name, date\n"
    metacsv=""
    for index,record in enumerate(SeqIO.parse(fafile, 'fasta')):
        if record.id.split("_|_")[2].count("-")==2:
            ldate=record.id.split("_|_")[2].split("|")[0].split("-")
            fract_date=toYearFraction(dt(int(ldate[0]), int(ldate[1]), int(ldate[2]), 0, 0, 0))
            datescsv+=record.id+", "+str(fract_date)+"\n"
            metacsv+=record.id+"\t"+str(fract_date)+"\t"+record.id.split("|")[3][:-1]+"\n"
        elif "Month_and_day_unknown" in record.id:
            datescsv+=record.id+", "+record.id.split("_|_")[2].split("__")[0]+"\n"
            metacsv+=record.id+"\t"+record.id.split("_|_")[2].split("(")[0][:-1]+"\t"+record.id.split("|")[3][:-1]+"\n"
        else:
            print ("incomplete date: " + record.id)
    return datescsv, metacsv


#Check file-names
SAMPLES=[]
for file in glob.glob(IFA+"*.fasta"):
    SAMPLES.append('.'.join(file.split("/")[-1].split(".")[:-1]))
        
# Rules ------------------------------------------------------------------------
# 

rule all:
    input:
        expand(workspace+'reads/{sample}_loc.fasta', sample=SAMPLES),
        expand(workspace+'reads'+str(nsamplesPerMonth)+'/{sample}_loc.fasta', sample=SAMPLES),
        workspace+'reads'+str(nsamplesPerMonth)+'/all.fasta',
        workspace+'reads'+str(nsamplesPerMonth)+'/dates.csv',
        workspace+'reads'+str(nsamplesPerMonth)+'/meta.csv',
        workspace+'reads'+str(nsamplesPerMonth)+'/RAxML_bestTree.tree.localRaxml',
        workspace+'reads'+str(nsamplesPerMonth)+'/treetime/divergence_tree.nexus'

#SINGLE READS
FASTA_DIR = IFA
PATTERN = '{sample}.fasta'

#AddRegion. Remove duplicates. Remove special characters
rule addRegion:
   input:
        fa=IFA+"{sample}.fasta",
   output:
        locfa=workspace+'reads/{sample}_loc.fasta',
   run:
        with open(output.locfa,'w') as fo:
            fo.write(str(addRegionNoDup(input.fa)))

#AddRegion. Remove duplicates. Remove special characters.
rule addRegionNoM:
   input:
        aus=ausfasta
   output:
        locfaaus=workspace+'reads/noM_loc.fasta',
   run:
        with open(output.locfaaus,'w') as fo:
            fo.write(str(addRegionNoDup(input.aus)))

#Select X samples per month
rule samplesPerMonth:
    input:
        locfa=workspace+'reads/{sample}_loc.fasta',
    output:
        locfaM=workspace+'reads'+str(nsamplesPerMonth)+'/{sample}_loc.fasta'
    params:
        n=nsamplesPerMonth
    run:
        with open(output.locfaM,'w') as fo:
            fo.write(str(samplesXmonth(input.locfa, params.n)))

#Merge previous results and align
rule mergeAlign:
    input:
        samples=expand(workspace+'reads'+str(nsamplesPerMonth)+'/{sample}_loc.fasta',  sample=SAMPLES),
        locfaaus=workspace+'reads/noM_loc.fasta'
    output:
        merge=workspace+'reads'+str(nsamplesPerMonth)+'/all.fasta',
        alignment=workspace+'reads'+str(nsamplesPerMonth)+'/all.alg.fasta'
    params:
        threads=threads
    shell:"""
         cat {input.samples} {input.locfaaus} >  {output.merge}
         mafft --thread {params.threads} {output.merge} > {output.alignment}
     """

#CHECK AND FIX ALIGNMENT MANUALLY BEFORE MAKE THE TREES (Or annotate and trim automatically)

#Extract metadata
rule extractMeta:
    input:
        alignment=workspace+'reads'+str(nsamplesPerMonth)+'/all.alg.fasta'
    output:
        dates=workspace+'reads'+str(nsamplesPerMonth)+'/dates.csv',
        meta=workspace+'reads'+str(nsamplesPerMonth)+'/meta.csv'
    run:
        dates, meta = extractMeta(input.alignment)
        with open(output.dates,'w') as fo:
            fo.write(dates)
        with open(output.meta,'w') as fo:
            fo.write(meta)

#Phylogenies
rule treeRAxML:
    input:
        fa=workspace+'reads'+str(nsamplesPerMonth)+'/all.alg.fasta'
    output:
        raxml=workspace+'reads'+str(nsamplesPerMonth)+'/RAxML_bestTree.tree.localRaxml'
    params:
        ofolder=workspace+'reads'+str(nsamplesPerMonth)+"/",
        ofile="tree.localRaxml",
        bs=bs,
        threads=threads
    shell:"""
        raxmlHPC-PTHREADS -T {params.threads} -m GTRGAMMA -s {input.fa} -w {params.ofolder} -n {params.ofile} -p 12345 -x 12345 -f a -# {params.bs}
    """

#Phylogenies
rule treetime:
    input:
        fa=workspace+'reads'+str(nsamplesPerMonth)+'/all.alg.fasta',
        dates=workspace+'reads'+str(nsamplesPerMonth)+'/dates.csv',
        raxml=workspace+'reads'+str(nsamplesPerMonth)+'/RAxML_bestTree.tree.localRaxml'
    output:
        treetime=workspace+'reads'+str(nsamplesPerMonth)+'/treetime/divergence_tree.nexus'
    params:
        ofolder=workspace+'reads'+str(nsamplesPerMonth)+"/"
    shell:"""
        treetime --aln {input.fa} --tree {input.raxml} --dates {input.dates} --outdir {params.ofolder}treetime
    """

