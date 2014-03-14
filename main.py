# A relaxed motif finding problem
'''
Given a postive (or negtive) dataset, find the top k motifs in m length


'''
pos_word={}
neg_word={}
import sys
m=int(sys.argv[3])
def count_pos(seq):
	global m
	global pos_word
	for i in range(len(seq)-m+1):
		key=seq[i:i+m]
		if pos_word.has_key(key):
			pos_word[key]+=1
		else:
			pos_word[key]=1

def count_neg(seq):
	global m
	global neg_word
	for i in range(len(seq)-m+1):
		key=seq[i:i+m]
		if neg_word.has_key(key):
			neg_word[key]+=1
		else:
			neg_word[key]=1


from Bio import SeqIO as io

# read in
handle = open(sys.argv[1])
pos = io.parse(handle,"fasta")
handle = open(sys.argv[2])
neg = io.parse(handle,"fasta")
# .id .seq
#pos_size = 0
#neg_size = 0
#baseline = 0.2*pos_size*()
for seq in pos:
	#print seq.seq
	#exit()
	count_pos(str(seq.seq))
	#pos_size+=1
for seq in neg:
	count_neg(str(seq.seq))
	#neg_size+=1
#print pos_word
#print neg_word
diff_word={}
#enrich_coef = 0.3
#diff_coef = 0.7	
for key in pos_word.keys():
	if neg_word.has_key(key):
		diff_word[key] = pos_word[key]*(pos_word[key]-neg_word[key])
	else:
		diff_word[key] = pos_word[key]*pos_word[key]

	
for w in sorted(diff_word, key=diff_word.get, reverse=True):
	print w, diff_word[w]