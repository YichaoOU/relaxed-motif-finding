from Bio import SeqIO as io
from StringIO import StringIO as si
from Bio import pairwise2 as pa
import sys
import os
import re
motif={}
#motif[motif]=array
#array=[length,mismatch,gap,occurrence]

def sub(m):
	motif=""
	mis=0
	gap=0
	match=0
	for l in m:
		if l == "-":
			gap = 1
		elif l == ".":
			mis += 1
		elif l == "|":
			match += 1
		else:
			motif += l
	#print "motif:",motif
	return [motif,match,mis,gap]		


#re.sub(r'[^\w]', ' ', s)
def count_motif(list):
	global motif
	for m in list:
		[key,match,mis,gap] = sub(m)
		if match <= 2:
			break
		if motif.has_key(key):
			motif[key][2] += 1
		else:
			motif[key]=[0,0,0,0]
			motif[key][0]=match
			motif[key][1]=mis
			motif[key][2]=gap
			motif[key][3]=1

			
# local alignment by water and parse file
def water(seq1,seq2):
	command = "/bin/bash -c \"water -asequence <( echo \""+ seq1 +"\" ) -bsequence <( echo \""+seq2+"\" ) -gapopen 10.0 -gapextend 0.5 -outfile temp\""
	os.system(command)
	local = open("temp").readlines()				
	final = ""				
	for i in range(31,len(local),4):				
		#print "local ",i," length: ",len(local[i])			
		b = local[i+1][21:71]			
		a = local[i][21:71]								
		if len(local[i]) <= 78:			
			a = a[:len(b)]
			#print "i:",i
			#print "a:",a
			#print "b:",b
			for j in range(len(b)):		
				if b[j] == " ":	
					final += "-"
				else:	
					final += a[j]+b[j]
			break		
		for j in range(len(b)):			
			if b[j] == " ":		
				final += "-"	
			else:		
				final += a[j]+b[j]	
	#print "final: ",final
	count_motif(re.split(r'--+',final[:-1]))

	
# read in
handle = open(sys.argv[1])
pos = list(io.parse(handle,"fasta"))
pos_size = len(pos)
for i in range(pos_size):
	for j in range(i+1,pos_size):
		water(str(pos[i].seq),str(pos[j].seq))
	
#reverse_complement
#seq1 = str(pos[1].seq)
#seq2 = str(pos[2].seq)
'''
# water -asequence <( echo "AAAA" ) -bsequence <( echo "TTT" ) -gapopen 10.0 -gapextend 0.5 -outfile s
command = "/bin/bash -c \"water -asequence <( echo \""+ seq1 +"\" ) -bsequence <( echo \""+seq2+"\" ) -gapopen 10.0 -gapextend 0.5 -outfile temp\""
os.system(command)
local = open("temp").readlines()
final = ""
for i in range(31,len(local),4):
	
	b = local[i+1][21:71]
	a = local[i][21:71]
	#print a
	#print b
	
	
	if len(local[i]) <= 77:
		a = a[:len(b)]
		for j in range(len(b)):
			if b[j] == " ":
				final += "-"
			else:
				final += a[j]+b[j]
		break
	for j in range(len(b)):
		if b[j] == " ":
			final += "-"
		else:
			final += a[j]+b[j]
'''
#print final



for key in motif.keys():
	print key,motif[key][0],motif[key][1],motif[key][2],motif[key][3]
	
