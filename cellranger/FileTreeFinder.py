#Searches through subdirectories of Arg1 path to find files matching Arg3
#Output Frame path is arg 2. Writes file of containing dir and file
#usage python3 FileTreeFinder.py PathIn PathOut Regex
import sys
import os
import pathlib
import pandas
import re

#sys.argv.append("~/Downloads")
#sys.argv.append("~/Desktop/UseFiles.txt")
#sys.argv.append(".txt")

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
            splitted=re.split("_",file)
            df=df.append({'path':str(os.path.join(os.sep,*path)),
                'name':'_'.join(splitted[:len(splitted)-4]),
                'file':str(file)},
                ignore_index=True)
df.to_csv(path_or_buf=os.path.expanduser(sys.argv[2]),index=False,sep="\t")
