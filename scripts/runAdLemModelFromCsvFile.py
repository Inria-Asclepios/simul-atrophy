#!/usr/bin/env python
# This script reads a configuration file that specifies all
# atrophy maps, masks, lambdaType options(scalar and/or tensor) and lambdaFiles for which
# the AdLemModel is to be run.
# The first word of each line is the heading followed by the filenames/options which are inputs
# corresponding to that heading. Each element in the row should be separated by a single whitespace.
# The headings (first column) are: atrophy mask lambdatype lambdaFile
# This script will create all the combinations possible for the provided options and files.
# And calls a script that runs the program AdLemModel for each of these combinations of options and files.

# Import csv
import csv #module to read the contents of file and put it into a dictionary
import sys 
import itertools #to get cartesian product of different lists.
import subprocess #for shell executable command
import optparse #to parse the command line input to this script.

#parse the command line inputs
parser = optparse.OptionParser()
parser.add_option('-p', '--patientID', dest='patient', help='patient ID name which is also a valid directory containing the atrophy map')
parser.add_option('-c', '--configFile', dest='configFile', help='configFile to be read')
parser.add_option('-j', '--jobName', dest='jobName', help='pbs jobscript filename')

(options, args) = parser.parse_args()
if options.patient is None:
    options.patient = raw_input('Enter a valid existing patient directory name')
if options.configFile is None:
    options.configFile = raw_input('Enter a valid existing configuration file')
if options.jobName is None:
    options.jobName = raw_input('Enter a valid existing pbs jobscript name')


inFile = open(options.configFile,'r') #open the csv file

# Start empty dictionary
inConfig = {}

# Loop over each line
for line in inFile.readlines():
    inConfig[line.split()[0]] = line.split()[1:] # get a list of words with split, then create dictionary
    # with first word as key and the remaining words as the list.

# Get all the combinations possible except using lambdaTypes.
tensorComb = list(itertools.product(inConfig['atrophy'],inConfig['mask'],inConfig['lambdaFile']))
scalarComb = list(itertools.product(inConfig['atrophy'],inConfig['mask']))

commandBase = './runAdLemModelCluster.py -p ' + options.patient + ' -j ' + options.jobName

# Iterate the relevant combinations for two possible lambdaTypes: scalar and/or tensor.
for tensorOrScalar in inConfig['lambdaType']:
    if tensorOrScalar == 'scalar':
        for choices in scalarComb:
            command = commandBase + ' -a ' + choices[0] + ' -m ' + choices[1] + ' -u false -l nofile -f ' + choices[0].split('.')[0]+'_'+choices[1].split('.')[0]+'_scalar_'
            #print command
            subprocess.call(command, shell=True)
    elif tensorOrScalar == 'tensor':
        for choices in tensorComb:
            command = commandBase + ' -a ' + choices[0] + ' -m ' + choices[1] + ' -u true -l ' + choices[2] + ' -f ' + choices[0].split('.')[0]+'_'+choices[1].split('.')[0]+'_'+choices[2].split('.')[0]+'_'
            #print command
            subprocess.call(command, shell=True)

inFile.close()