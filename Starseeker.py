from Bio import SeqIO
import collections
import RNA
import timeit
import datetime
import sys
import operator
import argparse

def string_to_pairs(bracket):
	open_pos = []
	result = {}
	for i in range(0, len(bracket)):
		if bracket[i] == '(':
			open_pos.append(i+1)
		elif bracket[i] == ')':
			last_open = open_pos.pop()
			result[last_open] = i+1
		elif bracket[i] == ".":
			result[i+1] = "*"
	reverse = {value: key for key, value in result.items()}
	final = result.copy()
	final.update(reverse)
	del final['*']
	return final

def get_star(precursor, mature, pairs):
	start = precursor.find(mature) + 1
	if start == 0:
		return("NOT FOUND")
	finish = start + len(mature) - 1
	extend_start = 0
	extend_finish = 0
	while pairs[start] == "*":
		start = start+1
		extend_start = extend_start + 1
		if extend_start == len(mature):
			return "error during the process"
	while pairs[finish] == "*":
		finish = finish - 1
		extend_finish = extend_finish + 1
		if extend_finish == len(mature):
			return "error during the process"
	star = precursor[pairs[finish] - extend_finish + 1:pairs[start] + extend_start + 2]
	if star == "":
		return "error during the process"
	elif len(star) > 40: 
		return "error during the process"
	else:
		return star

if __name__ == '__main__':

	usage = "Program that finds star..."
	parser = argparse.ArgumentParser(description=usage)

	parser.add_argument('-m', '--m', '-mir', dest='mir',
	                   help='file w/ mature miRNA sequences (FASTA format)')
	parser.add_argument('-p', '--p', '-pre', dest='pre',
	                   help='file w/ pre-mirRNA sequences (FASTA format)')
	parser.add_argument('-o', '--o', '-out', dest='output', default="output.txt",
	                   help='output file')
	parser.add_argument('-l', '--l', '-log', action='store_true', dest="log",
		               help="If this option selected the log file is created")

	args = parser.parse_args()

	pre_file = args.pre
	mat_file = args.mir

	start = timeit.default_timer()
	output = open(args.output, "w")
	
	handle1 = open(pre_file, "rU")
	precursors = list(SeqIO.parse(handle1, "fasta"))
	handle1.close

	handle2 = open(mat_file, "rU")
	matures = list(SeqIO.parse(handle2, "fasta"))
	handle2.close

	preseqs = []
	matseqs = []
	headers = []

	for record in precursors:
		preseqs.append([str(record.id), str(record.seq)])

	for record in matures:	
		matseqs.append(str(record.seq))
		headers.append(str(record.id))

	howmany = {key: 0 for key in matseqs}

	entries = []

	for mat in matseqs:
		for i in range(0, len(preseqs)):
			if mat in preseqs[i][1]:
				entries.append([preseqs[i][0], preseqs[i][1], mat, " "])
				i = i + 1
				howmany[mat] = howmany[mat] + 1
			else:
				continue
	
	helplist = []
	for record in matures:
		helplist.append([str(record.id), str(record.seq)])
	final = {key: None for key in headers}
	
	for record in matseqs:
		for i in helplist:
			if record == i[1]:
				final[i[0]] = howmany[i[1]]
	
	final_sorted = sorted(final.items(), key=operator.itemgetter(1), reverse = True)
	
	for entry in entries:
		result = RNA.fold(entry[1])
		entry[3] = result[0]
	
	e = []
	for entry in entries:
		if entry not in e:
			e.append(entry)
		else:
			continue	
	
	for entry in e:
		output.write(">" + entry[0] + " *star*")
		output.write('\n')
		output.write(get_star(entry[1], entry[2], string_to_pairs(entry[3])))
		output.write('\n')
	
	output.close
	stop = timeit.default_timer()
	if args.log == True:
		log = open("log.txt", "w")
		log.write(str(len(preseqs)) + " sequences existed in precursor file.")
		log.write('\n')
		log.write(str(len(matseqs)) + " sequences existed in matures file.")
		log.write('\n')
		log.write(str(len(entries)) + " pairs were created.")
		log.write('\n')
		log.write(str(len(e)) + " stars were created.")
		log.write('\n')
		log.write("The runtime of the program was ")
		runtime = stop - start
		precision = 6
		runtime = round(runtime * 10**precision) / 10**precision
		log.write(str(runtime) + " seconds.")
		log.write('\n')
		now = datetime.datetime.now()
		log.write("The analysis was performed on " + str(now) + ".")
		log.write('\n')
		log.write('\n')
		log.write("miRNA - times paired")
		log.write('\n')
		for record in final_sorted:
			log.write(str(record[0]) + " " + str(record[1]))
			log.write('\n')
		log.close()		
