#!/usr/bin/env python3
import sys, os, argparse
from operator import itemgetter
import random

parser=argparse.ArgumentParser(usage='''\
FiterContig.py
''')
parser.add_argument('-k', type=int, default=2, help='minimum length of contigs')
parser.add_argument('-n', type=int, default=20, help='the number of selected features')
#parser.add_argument('-P', type=str, default="N", help='prefix of first group')
#parser.add_argument('-T', type=str, default="T", help='q-value threshold')
parser.add_argument('inMat', help='input matrix file')
parser.add_argument('outMat', help='output matrix file')
args=parser.parse_args()
inMat = open(args.inMat, 'r')
#outMat = open(args.outMat, 'w')
#upper = args.U
#lower = args.L

IDs = []
for line in inMat.readlines():
	linn = line.rstrip().split("\t")
	IDs.append(linn[0])
half = int(len(IDs)/2)

group1= []
group2= []

random.shuffle(IDs)
#print(newIDs)
group1 =set( IDs[0:args.k])
group2 = set(IDs[args.k:args.k*2])
group1_sm = set(IDs[args.k*2:args.k*3])
group2_sm = set(IDs[args.k*3:args.k*4])
	
filelist =[]
groupA = group1
groupB = group2
groupA_sm = group1_sm
groupB_sm = group2_sm

#print(group1)
#print(group2)

print(group1)
print(group2)

for i in range(args.n):
	outfilename = args.outMat+"_"+str(i)+".tsv"
	outfile = open(outfilename, "w")
	filelist.append(outfilename)

	if i == args.n/2:
		groupA = group2
		groupB = group1
		groupA_sm = group2_sm
		groupB_sm = group1_sm

	for j in range(len(IDs)):
		if IDs[j] in groupA:
			abundance = random.uniform(15.0, 20.0)
		elif IDs[j] in groupB:
			abundance = random.uniform(0.0, 5.0)
		elif IDs[j] in groupA_sm:
			abundance = random.uniform(10.0, 15.0)
		elif IDs[j] in groupB_sm:
			abundance = random.uniform(5.0, 10.0)
		else:
			abundance = random.uniform(4.0, 7.0)
		#index = random.sample(range(len(header)-1), threshold)
	


		outfile.write(IDs[j]+"\t"+str(abundance)+"\n")
	outfile.close()

print(filelist)








































