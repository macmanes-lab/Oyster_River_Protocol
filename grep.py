import re
import sys

# python wanted.file input.file

with open(sys.argv[1],'r') as f:
    queries = [l.strip() for l in f]

with open(sys.argv[2],'r') as f:
    for line in f:
        for query in queries:
            if query in line:
                print line


a=set(open(sys.argv[1],'r'))
for line in a:
    for query in queries:
        if query in line:
            print query 
