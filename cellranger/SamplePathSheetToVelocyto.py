#Searches through subdirectories of Arg1 path to find files matching Arg3
#Output Frame path is arg 2. Writes file of containing dir and file
#usage python3 FileTreeFinder.py PathIn PathOut Regex
import sys
import os
import pathlib
import pandas
import re
from itertools import compress

#Usage
sys.argv.append("~/Desktop/FQfiles.txt")
sys.argv.append("~/Desktop/MacaqueBrainVelocytoNoMaskRun.sh")
sys.argv.append("/ye/yelabstore2/mtschmitz/refdata-celranger-mmul8-toplevel")
sys.argv.append("$SCRATCH_JOB/")
#sys.argv.append("/ye/yelabstore2/mtschmitz/seq/AlignedMacaqueBrainSeq/Exonic/")

df =  pandas.read_csv(os.path.expanduser(sys.argv[1]),sep="\t")
csvfile = open( os.path.expanduser(sys.argv[2]), 'w')

csvfile.write('#!/bin/bash\n\
#$ -o ~/log\n\
#$ -e ~/log\n\
#$ -cwd\n\
#$ -j y\n\
#$ -pe smp 1\n\
#$ -l mem_free=32G\n\
#$ -l h_rt=300:00:00\n\n\n')

#csvfile.write('SCRATCH_JOB=/scratch/$USER/jobs/$JOB_ID\n\

csvfile.write('SCRATCH_JOB=/wynton/scratch/mschmitz1/Exonic\n\
export PATH="/ye/yelabstore2/mtschmitz/utils/samtools-1.8:$PATH"\n\
mkdir -p $SCRATCH_JOB\n\
cd $SCRATCH_JOB/\n\n')


pathdone=[]
namedone=[]
for index, row in df.iterrows():
    if row['path'] not in pathdone:
        outputName=row['name']+"_Out"
        pathsplit=row['path'].split(os.sep)
        if row['name'] in namedone:
            matchingIndices = [x == row['name'] for x in namedone]
            matchPaths=list(compress(pathdone, matchingIndices))
            nestList=[x.split(os.sep) for x in matchPaths]
            flatList = [item for subList in nestList for item in subList]
            pathsplitBool=[x not in flatList for x in pathsplit]
            uniquePartsOfPath= "-".join(list(compress(pathsplit, pathsplitBool)))
            outputName = outputName+"_"+uniquePartsOfPath
        csvfile.write("##"+row['name']+"\n")
        csvfile.write("/ye/yelabstore2/mtschmitz/utils/python3/bin/velocyto run10x ")
        csvfile.write(sys.argv[4]+outputName+" ")
        csvfile.write(sys.argv[3]+"/genes/genes.gtf")
        csvfile.write("\n\n")
        namedone.append(row['name'])
        pathdone.append(row['path'])

csvfile.close()
