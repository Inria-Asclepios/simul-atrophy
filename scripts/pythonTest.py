#!/usr/bin/env python
import subprocess
#import sys

lists = subprocess.check_output("cat ../src/lists",shell = True) 
names = lists.split()

#print lists
#print names
for name in names:
    subprocess.call("echo " + name, shell = True)

