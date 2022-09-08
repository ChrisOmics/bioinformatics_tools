#!/usr/bin/python
################################
#Created by: Chris Robles
#Date: 06/11/16
#
#Haplotype Aligner
#Input: reads file (command line argument)
#Output: Haplotypes with baseline method and markov chain method
###############################


import sys
import re
import numpy as np
import os


infile = open(sys.argv[1], 'r')

haplofn = (sys.argv[1]).replace("reads", "haplotypes")
#try:
	##print "True Haplotypes:"
	#os.system("cat %s"%haplofn)
#except:
	#print("cant find haplotype file")

def read_store(infile):
	#make an 2 d array of each read and position
	read_id = 0
	rl=6
	posReadDict = dict()
	posReadlist = list()
	for line in infile:
		haplo_length = len(line.strip())
		r = re.search(r'\d+', line.strip())
		read = r.group()
		pos = r.start()

		#posReadDict[read_id] = [''.join(read),pos]
		posReadlist.append([read,pos])

	#print posReadlist

	return(haplo_length, posReadlist)


def haplo_assembler_baseline(readlist):
	#put read into one bin at each position, and then add to a different bin if it is a new read
	#initiialize variables
	haplotype1=""
	haplotype2=""
	prev_read = readlist[0][0]
	prev_pos = 0
	firstrd = ""
	secread = ""
	for read in readlist:

		pos = int(read[1])
		rd = str(read[0])
		if rd != prev_read and pos != prev_pos:
			haplotype1 = haplotype1[:pos] +  rd + haplotype1[pos+len(rd):]
			firstrd = rd
		if rd != prev_read and pos == prev_pos and rd != firstrd:
			haplotype2 = haplotype2[:pos] +  rd + haplotype2[pos+len(rd):]
			secread = rd
		elif rd != prev_read and pos == prev_pos and rd == firstrd:
			continue
			#print rd
			#print haplotype2

		elif rd != prev_read and pos == prev_pos and rd == secread:
			continue
		elif rd == prev_read and pos == prev_pos:
			continue
		elif rd == prev_read and pos != prev_pos and rd == firstrd:
			haplotype1 = haplotype1[:pos] +  rd + haplotype1[pos+len(rd):]
			firstrd = rd
		elif rd == prev_read and pos != prev_pos and rd == secread:
			haplotype2 = haplotype2[:pos] + rd + haplotype2[pos+len(rd):]
			secread = rd
		prev_read = rd
		prev_pos = pos
	print haplotype1
	print haplotype2

def tran_matrix(hlen, readlist):
	##try to figure out how to make list of list of As, Bs, etc
	#A is 0 to 0, B is 0 to 1, C is 1 to 0, D is 1 to 1

	#list of dictionaries with transition probabilities


	#then go throguh read and add at each position
	tm = np.zeros((hlen,4))
	pos = 0
	read = []
	for read in readlist:
		pos = read[1]
	#	print pos
		rd = str(read[0])
		for i, base in enumerate(rd):
		#i keeps track of where in read
			try:
				#next read is saved
				next_rd = rd[i+1]
				#calculate transition probability
				if int(rd[i]) == 0 and int(next_rd) ==0:
					tm[pos+i][0]+=1
				elif int(rd[i]) == 0 and int(next_rd) ==1:
					tm[pos+i][1]+=1
				elif int(rd[i]) == 1 and int(next_rd) ==0:
					tm[pos+i][2]+=1
				elif int(rd[i]) == 1 and int(next_rd) ==1:
					tm[pos+i][3]+=1
	#
	#
	#
					continue
			except IndexError:
				continue
	#print trans_matrix
	return tm


def haplotype_builder(transmatrix, firsthap):
	tm = transmatrix
	base = firsthap
	haplotype=[str(firsthap)]
	for pos in tm:
		prob_zz = pos[0]/sum(pos)
		prob_zo = pos[1]/sum(pos)
		prob_oz = pos[2]/sum(pos)
		prob_oo = pos[3]/sum(pos)
		if base == 0 and prob_zz > prob_zo:
			newbase = 0
			haplotype.append(str(newbase))
			base = newbase
		elif base == 0 and prob_zz < prob_zo:
			newbase = 1
			haplotype.append(str(newbase))
			base = newbase
		elif base ==1 and prob_oz > prob_oo:
			newbase = 0
			haplotype.append(str(newbase))
			base = newbase
		elif base ==1 and prob_oz < prob_oo:
			newbase = 1
			haplotype.append(str(newbase))
			base = newbase


	print "".join(haplotype)


def main():
	print "Baseline:"
	hlen, rlist = read_store(infile)
	haplo_assembler_baseline(rlist)
	print "Method:"
	tm = tran_matrix(hlen, rlist)
	haplotype_builder(tm, 0)
	haplotype_builder(tm,1)
if __name__ == '__main__':
	main()
#then store the string in a dictionary as position:read
#but if the position is the same, then have to rename it maybe

#if pos and read exsists, add a .1 to name. else continue
#then check all positions and see if 1 levenschtein distance away.
#also add a counter for how many reads ther are
#if 1 levenstein distance way then remove read


#another way is to go through each line, if read is 1 lev away from last read, then add to haplotype.
#else add to haplotype2
