# -*- coding: cp1257 -*-
#Last update: 04.04.18
#Calculate Henikoff Weights 
#Steven Henikoff and Jorja G. Henikoff (1994) "Position-based Sequence Weights" 

import sys, getopt

def getDictItems(dictionary):
	#If python 2
	if sys.version_info[0] < 3:
		arr = dictionary.iteritems()
	else:
		arr = dictionary.items()

	return arr
		
def writeListToFile(filename,list):
	outputfile = open(filename, "w")
	for i in range (0,len(list)):
		outputfile.write(str(list[i][0])+"\t"+str(list[i][1])+"\n")
	outputfile.close()

def readFastaToArray(filePath):
	try:
		file = open(filePath, "r")
	except:
		print ("ERROR - Reading " + str(filePath))

	list = []
	entry = []
	while True :
		
		line = file.readline()
		if line == "" :
			#last seq
			entry.append(seq)
			list.append(entry)
			break
		
		lineList = line.split("\n")
		line = lineList[0]
		if(">" in line):
			#previous seq
			if(len(entry)>0):
				entry.append(seq)
				list.append(entry)
			
			seq = ""
			name = line.split(">")[1]
			entry = []
			entry.append(name)
		else:
			seq = seq + line
		
	file.close()
	return list

#Weights based on codons
def calculateHenikoffWeights(sequences):
	weights = []
	sequences_count = len(sequences) 
	msa_length = len(sequences[0][1]) #length of the first sequence in a MSA
	for i in range (0, sequences_count):
		weights.append([sequences[i][0],0])

	for column in range(0, msa_length, 3):
		codons = {}
		for i in range (0, sequences_count):
			if (not len(sequences[i][1]) == msa_length):
				print("ERROR - Sequences are with different length")
				sys.exit()

			sequence = sequences[i][1]
			codon = sequence[column]+sequence[column+1]+sequence[column+2]
			codon = codon.lower()
			if (codons.get(codon) is not None):
				codons[codon] +=1
			else:
				codons[codon] = 1

		for i in range (0, sequences_count):
			sequence = sequences[i][1]
			codon = sequence[column]+sequence[column+1]+sequence[column+2]
			codon = codon.lower()
			r = len(codons) #r is the number of different residues in the position a
			s = codons.get(codon) #s is the number of times the particular residue appears in the position.
			weight = 1/(r*s)
			weights[i][1] += weight
	
	return weights
		
# Weights based on single nucleotide	
def calculateHenikoffWeightsSingleNucleotide(sequences):

	weights = []
	sequences_count = len(sequences) 
	msa_length = len(sequences[0][1])
	for i in range (0, sequences_count):
		weights.append([sequences[i][0],0])

	for column in range(0, msa_length):
		nucleotides = {}
		for i in range (0, sequences_count):
			if (not len(sequences[i][1]) == msa_length):
				print("ERROR - Sequences are with different length")
				sys.exit()
			sequence = sequences[i][1]
			nucleotide = sequence[column]
			nucleotide = nucleotide.lower()
			if (nucleotides.get(nucleotide) is not None):
				nucleotides[nucleotide] +=1
			else:
				nucleotides[nucleotide] = 1
		
		for i in range (0, sequences_count):
			sequence = sequences[i][1]
			nucleotide = sequence[column]
			nucleotide = nucleotide.lower()
			r = len(nucleotides) #r is the number of different residues in the position a
			s = nucleotides.get(nucleotide) #s is the number of times the particular residue appears in the position.
			weight = 1/(r*s)
			weights[i][1] += weight
	return weights

def normalizeHenikoffWeights(weights):
	sum = 0
	for i in range (0, len(weights)):
		sum += weights[i][1]
	for i in range (0, len(weights)):
		weights[i][1] = weights[i][1] / sum
	return weights
	

def usage():
	print('heinikoff_weights.py -i <codon alignment in FASTA format>')
	print('-i, --input \t a codon alignment in FASTA format')
	print('-o, --output \t output file, default is weights.txt')
	print('-v, --verbose \t print process steps')
	
def main():
	try:
		opts, args = getopt.getopt(sys.argv[1:], "h:i:o:v", ["help", "input=","output="])
	except getopt.GetoptError as err:
		# print help information and exit:
		print(err) 
		usage()
		sys.exit(2)

	for opt, arg in opts:
		if opt == "-v":
			global verbose
			verbose = True
		elif opt in ("-h", "--help"):
			usage()
			sys.exit()
		elif opt in ("-i", "--input"):
			global input
			input = arg
		elif opt in ("-o", "--output"):
			global output
			output = arg
		else:
			assert False, "unhandled option"


if __name__ == '__main__':

	#Global variables
	input = None
	output = None
	verbose = False #A verbose mode is an option available - provides additional details as to what the computer is doing
	
	main()
	if (input):
		if(verbose):
			print('reading sequences to memory') 
		
		sequences = readFastaToArray(input)
		if(verbose):
			print('calculate weights') 
		rawResult = calculateHenikoffWeights(sequences)
		if(verbose):
			print('normalize weights') 
		result = normalizeHenikoffWeights(rawResult)
			
		if(output != None):
			writeListToFile(output,result)
		else:
			for i in range (0, len(result)):
				print(result[i][0] + "\t" + str(result[i][1]))
		
		if(verbose):
			print('completed') 

	else:
		print("ERROR - No input. See --help for more details")


	
