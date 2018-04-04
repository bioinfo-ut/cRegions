# -*- coding: cp1257 -*-
#Last update: 08.03.18
#Calculate Henikoff Weights  -> Steven Henikoff and Jorja G. Henikoff (1994) "Position-based Sequence Weights" 

import sys, getopt

def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False
		
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
			
			#new seq
			seq = ""
			name = line.split(">")[1]
			entry = []
			entry.append(name)
		else:
			seq = seq + line
		
	file.close()
	return list
	
def calculateHenikoffWeights(sequences):

	weights = []
	sequences_count = len(sequences) 
	msa_length = len(sequences[0][1])
	for i in range (0, sequences_count):
		weights.append([sequences[i][0],0])

	for column in range(0, msa_length):
		letters = [0,0,0,0,0] #A ,C ,G , T/U ,gap
		for i in range (0, sequences_count):
			sequence = sequences[i][1]
			letter = sequence[column]
			letter = letter.lower()
			if (letter == "a"):
				letters[0] += 1
			elif (letter == "c"):
				letters[1] += 1
			elif (letter == "g"):
				letters[2] += 1
			elif (letter == "t" or letter == "u"):
				letters[3] += 1
			elif (letter == "-"):
				letters[4] += 1

		r = sum(1 for x in letters if x > 0) # "r" is the number of different nucleotides
		weight = 1.0/r 
		
		for i in range (0, sequences_count):
			sequence = sequences[i][1]
			letter = sequence[column]
			letter = letter.lower()
			if (letter == "a"):
				weights[i][1] += weight/letters[0] # letters[0] aka "s" is the number of times the particular nucleotide appears in the position.
			elif (letter == "c"):
				weights[i][1] += weight/letters[1]
			elif (letter == "g"):
				weights[i][1] += weight/letters[2]
			elif (letter == "t" or letter == "u"):
				weights[i][1] += weight/letters[3]
			elif (letter == "-"):
				weights[i][1] += weight/letters[4]
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
	print('-i, --input \t input file, a codon alignment')
	print('-o, --output \t output file, default is weights.txt')
	print('-v, --verbose \t print process steps')
	
def main():
	try:
		opts, args = getopt.getopt(sys.argv[1:], "h:i:o:v", ["help", "input=","output="])
	except getopt.GetoptError as err:
		# print help information and exit:
		print(err)  # will print something like "option -a not recognized"
		usage()
		sys.exit(2)

	for opt, arg in opts:
		if opt == "-v":
			global verbose
			verbose = True
		elif opt in ("-h", "--help"):
			usage()
			sys.exit()
		elif opt in ("-i", "-in", "--input"):
			global input
			input = arg
		elif opt in ("-o", "-out", "--output"):
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
	if(verbose):
		print('reading sequences to memory') 
	sequences = readFastaToArray(input)
	if(verbose):
		print('calculate Henikoff Weights') 
	rawresult = calculateHenikoffWeights(sequences)
	if(verbose):
		print('normalize Henikoff Weights') 
	result = normalizeHenikoffWeights(rawresult)
		
	if(output != None):
		writeListToFile(output,result)
	else:
		for i in range (0, len(result)):
			print(result[i][0] + "\t" + str(result[i][1]))
	
	if(verbose):
		print('completed') 


	
