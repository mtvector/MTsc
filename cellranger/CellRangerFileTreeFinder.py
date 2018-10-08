#Searches through subdirectories of Arg1 path to find files matching Arg3
#Output Frame path is arg 2. Writes file of containing dir and file
#usage python3 FileTreeFinder.py PathIn PathOut Regex
import sys
import os
import pathlib
import pandas
import re

#sys.argv.append("~/code/data")
#sys.argv.append("~/Desktop/UseOutputs.txt")
#sys.argv.append("filtered_gene_bc_matrices_h5.h5")

#Compile the regex/search string for the hunt
searcher=re.compile(sys.argv[3])

myFile = pathlib.Path(os.path.expanduser(sys.argv[1]))
print(myFile.is_dir())

df = pandas.DataFrame(columns=['path','name','file'])

# traverse root directory, and list directories as dirs and files as files
for root, dirs, files in os.walk(os.path.expanduser(sys.argv[1])):
    path = root.split(os.sep)
    for file in files:
        if searcher.search(file):
            sampleName=path[len(path)-2]
            splitted=re.split("_",sampleName)
            df=df.append({'path':str(os.path.join(os.sep,*path)),
                'sample_name':splitted[0],
                'sample_full_run_name':'_'.join(splitted[0:]),
                'file':file},
                ignore_index=True)
df.to_csv(path_or_buf=os.path.expanduser(sys.argv[2]),index=False,sep="\t")
