#Searches through subdirectories of Arg1 path to find files matching Arg3
#Output Frame path is arg 2. Writes file of containing dir and file
#usage python3 FileTreeFinder.py PathIn PathOut Regex
import sys
import os
import pathlib
import pandas
import re
from itertools import compress

#Usage python3 SamplePathSheetToCellRanger.py PATHtoCellRanger ~/MacaqueE65Files.txt ~/MacaqueBrainE65JobRun.sh refdata-celranger-mmul8-toplevel /ye/yelabstore2/mtschmitz/seq/AlignedMacaqueBrainSeq/Exonic/
#Usage python3 SamplePathSheetToCellRanger.py ~/MosFiles.txt ~/MacaqueBrainE65JobRun.sh refdata-celranger-mmul8-toplevel /ye/yelabstore2/mtschmitz/seq/AlignedMacaqueBrainSeq/Exonic/

sys.argv.append("~/Desktop/FQfiles.txt")
sys.argv.append("~/Desktop/MacaqueBrainJobRun.sh")
sys.argv.append("refdata-celranger-mmul8-toplevel")
sys.argv.append("/ye/yelabstore2/mtschmitz/seq/AlignedMacaqueBrainSeq/Exonic/")
sys.argv.append("/home/mschmitz1/utils/cellranger-3.0.2/cellranger")
sys.argv.append("/home/mschmitz1/matthew/tmp/cellrangerrun")
sys.argv.append("/home/mschmitz1/matthew/refdata")


if False:
    header='#!/bin/bash\n\
#$ -o ~/log\n\
#$ -e ~/log\n\
#$ -cwd\n\
#$ -j y\n\
#$ -pe smp 2\n\
#$ -l mem_free=64G\n\
#$ -l h_rt=300:00:00\n\n\n'
else:
    header='#!/bin/bash\n\
#PBS -o /home/mschmitz1/log\n\
#PBS -e /home/mschmitz1/log\n\
#PBS -l nodes=1:ppn=4 -l vmem=148gb\n\n\n'

df =  pandas.read_csv(os.path.expanduser(sys.argv[1]),sep="\t")
csvfile = open( os.path.expanduser(sys.argv[2]), 'w')
print(df)
csvfile.write(header)
csvfile.write('SCRATCH_JOB='+sys.argv[6]+'\n\
export PATH="'+sys.argv[5] +':$PATH"\n\
mkdir -p $SCRATCH_JOB\n\
## 1. Copy input files from global disk to local scratch\n\
cp -r ' +sys.argv[7]+sys.argv[3]+ ' $SCRATCH_JOB/\n\
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
        csvfile.write("cp -r "+row['path']+" $SCRATCH_JOB/\n")
        csvfile.write("gunzip $SCRATCH_JOB/"+pathsplit[len(pathsplit)-1]+"/*.gz\n")
        csvfile.write(sys.argv[5] + " count ")
        csvfile.write("--id="+outputName+" ")
        csvfile.write("--fastqs=$SCRATCH_JOB/"+ pathsplit[len(pathsplit)-1] +" ")
        csvfile.write("--sample="+ row['name']+ " ")
        csvfile.write("--transcriptome=$SCRATCH_JOB/"+sys.argv[3]+" ")
        csvfile.write("--jobmode=sge\n")
        csvfile.write("\n")
        csvfile.write("mv $SCRATCH_JOB/"+outputName+" "+sys.argv[4]+"\n")
        csvfile.write("rm -r $SCRATCH_JOB/"+pathsplit[len(pathsplit)-1])
        csvfile.write("\n\n")
        namedone.append(row['name'])
        pathdone.append(row['path'])

csvfile.write('cd $SCRATCH_JOB/\ncd ..\nrm -rf $SCRATCH_JOB/')

csvfile.close()
