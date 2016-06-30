# This script calculates the Gini index and GC percent over 100bp windows,
# starting from the beginning of the coding region, and taking 100bp windows
# until the end of the coding region.

# A filtering process is used to discard invalid regions, which must satisfy the
# following criteria:
# 1. The minimum coverage on any base is 15
# 2. The average mismatch coverage on AC bases is at least 15

# Average mismatch coverage = Average NUMBER of mutations on AC positions

# The Gini index is calculated by the Gini over the AC mutation percentages, on 
# valid regions which pass filter.


# Import necessary packages:

import heapq
import matplotlib 
matplotlib.use('agg') 
import numpy as np
import scipy
from scipy import stats
import cPickle as pickle
import csv
import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages
import pysam as py


def main():
	
	print 'Calculating gini...'
	coding_gc, coding_gini, utr_gc, utr_gini, filter_discard, bad_genes, min_coverage, good_genes = get_gini('MegZ37_38.p')

	# Now, we plot the correlation between gini index and gc_percent in coding regions
	# only, in utr regions only and then combined
	print 'Regions not passing filters: ', filter_discard
	data_pts  = 0 
	with open('Good_genes_and_regions.txt','w') as fout:
		fout.write('@ Number of Discarded Genes: {} \n @Minimum Coverage: {}'.format(bad_genes, min_coverage))
		for gene in good_genes:
			good_regions = good_genes[gene]
			for tup in good_regions:
				data_pts += 1
				fout.write(">{}\t{}\t{}\n \n".format(gene, tup[0], tup[1]))
	print 'Number of data points: ', data_pts

	comb_gc = coding_gc + utr_gc
	comb_gini = coding_gini + utr_gini

	comb_gini=np.array(comb_gini)
	comb_gc = np.array(comb_gc)

	coding_gc = np.array(coding_gc)
	coding_gini = np.array(coding_gini)
	utr_gc = np.array(utr_gc)
	utr_gini = np.array(utr_gini)

	assert len(comb_gc) == len(comb_gini)
	pearson, p = stats.pearsonr((comb_gc), comb_gini)
	print 'overall p', p
	pearson_utr, p = stats.pearsonr(utr_gini, utr_gc)
	print 'utr p', p
	pearson_cod, p = stats.pearsonr(coding_gini, coding_gc)
	print 'coding p', p
	print pearson, 'overall', pearson_utr, 'utr', pearson_cod, 'coding'
	plt.scatter(coding_gini, (coding_gc), color = 'red')
	plt.scatter(utr_gini, (utr_gc), color = 'blue')
	print 'num points', len(coding_gini) + len(utr_gini)
	plt.grid(True)
	plt.ylabel('GC Percent')
	plt.xlabel('Gini Index')
	plt.title('Correlation between GC Percent and Gini Index in Coding Regions \n Pearson Correlation Coefficient = {} \n UTR r = {}, Coding r = {}'.format(pearson, pearson_utr, pearson_cod ))
	fig = plt.gcf
	plt.savefig('GiniCombined_withp.eps', format='eps', dpi=1000)
	plt.clf()

	return None

def get_gini(pickle_file):
	# Define variables
	inc = 50 				#Size of sliding window 
	mut_filter = 20			#Minimum mismatch average
	count_filter = 100		#Minimum coverage filter
	min_coverage = 150		#Minimum average coverage for a gene
	num_genes = 100 			#Consider only these many genes, with the best coverage
	# Define options 
	use_coding_annotation = True
	normalize_untreated = False

	# Define input files 
	#pickle_file = 'MegZ37_TruSeq1_S79_L007_R1_001.strip.human.canon.mk.sam_Mutations_positive.p'									#Pickled dictionary of mutations 
	untreated_pickle_file = 'fxr2_untreated_pickle.p'				#Pickled dictionary of mutations on untreated control
	ref_file = './human_transcriptome/human_trans.fasta.txt'		#Contains reference sequences
	start_file = 'human.canonical.cds.txt'							#Contains start and stop position of coding regions for each gene

	# Load dictionary of mutations on the positive strand
	mut_dict = pickle.load( open(pickle_file, "rb" ) )

	if normalize_untreated:
		untreated_mut_dict = pickle.load( open(untreated_pickle_file, 'rb'))

	# If already calculated, open reference and coding region dictionary
	try: 
		ref_dict = pickle.load(open('ChrPos_to_Base.p', 'rb'))
		pos_dict = pickle.load(open('Gene_to_StartStop.p', 'rb'))
	except getattr(__builtins__,'FileNotFoundError', IOError):
		ref_dict = get_reference_sequence(ref_file)
		pos_dict = get_start_and_end(start_file)

	# Initialize variables

	bad_genes = 0			#Counts how many genes have been discarded
	good_genes = {}			#Dictionary containing genes and the [(start,stop),...(start,stop)] positions of valid regions
	utr_gini = []			#List of gini indexes on UTR windows
	utr_gc = []				#List of GC percent on UTR windows
	coding_gini = []		#List of gini indexes on coding windows
	coding_gc = []			#List of GC percent on coding windows
	best_genes = []
	filter_discard = 0 
	
	# Obtain the best covered genes 
	genes = []
	gene_cvg = []
	for gene in mut_dict:
		genes.append(gene)
		gene_cvg.append(np.sum(np.array(mut_dict[gene]['mutations_counts'][:-2])))
	genes = np.array(genes)
	gene_cvg = np.array(gene_cvg)
	idx = np.argpartition(gene_cvg, -num_genes)
	best_genes = genes[idx[-num_genes:]]

	for gene in best_genes:	
		mut_values = mut_dict[gene]['mutations_counts'][:-2]
		total_counts = mut_dict[gene]['total_counts'][:-2]

		if normalize_untreated:
			ctl_mut_values = untreated_mut_dict[gene]['mutations_counts'][:-2]
			ctl_total_counts = untreated_mut_dict[gene]['total_counts'][:-2]

		# Discard genes with less than the minimum average coverage
		if np.mean(np.array(total_counts)) < min_coverage:
			bad_genes += 1
		else:
			if use_coding_annotation:
				# Define the beginning and end of the coding region
				try:
					start, end = pos_dict[gene]
					start = int(start)
					end = int(end)
				except (KeyError, ValueError):
					print 'Key Error or ValueError accessing start and end of coding region'
					start, end = (0,len(mut_values))
			else:
				start = 0
				end = len(mut_values)	


			# Keep a list of which positions are G's or C's
			gc_count = [0]*(len(mut_values))

			for pos in range(len(mut_values)):
				#Determine the base at that position
				try:
					base = ref_dict[gene][pos]
				except IndexError:
					print 'IndexError accessing base position {} on gene {} in reference dictionary'.format(pos, gene)

				# Only consider mutations occurring at A or C. Add mutation value and count at that position
				if base.upper() not in set(['A','C']):
					mut_values[pos] = -1.0
				if base in set(['C','G']):
					gc_count[pos] = 1

			# Now, calculate the Gini index and GC percent along sliding windows on this gene

			# These indices are the start of 100bp windows 
			coding_indices = [(i, i+inc, True) for i in range(start,end,inc)]
			if len(coding_indices) > 0:
				print 'Coding gene!: ', gene, len(coding_indices)
			utr_indices=[(i, i + inc, False) for i in range(0, start, inc)] + [(i, i+inc, False) for i in range(end, len(mut_values), inc)]

			for (start_pos, end_pos, coding) in coding_indices + utr_indices:
				valid_region = True
				# Get the windows of exactly <inc>bp
				if start_pos < 0:
					valid_region = False
				try:
					counts = total_counts[start_pos: end_pos]
					mutations = mut_values[start_pos: end_pos]
					gc = np.array(gc_count[start_pos:end_pos])
				except KeyError:
					valid_region = False

				# Now test to see if these regions pass all the filters

				if valid_region:
					if counts == []:
						valid_region = False
					else:
						observed = np.array([mutations[i]*1.0/(max(1,counts[i])) for i in range(len(counts))])
						if normalize_untreated:
							for pos in range(len(observed)):
								try:
									observed[pos] *= (1.0*ctl_total_counts[pos])/ctl_mut_values[pos]
								except ZeroDivisionError:
									pass
						idx = np.where((observed < 0.15) & (observed >= 0))
						observed = observed[observed >= 0]
						observed = observed[observed < 0.15]
						observed = np.sort(observed)
						observed = observed.tolist()

						if observed == [] or min(counts) < count_filter or np.mean(np.array(mutations)[idx]) < mut_filter or min(heapq.nlargest(1,[mutations[i]*100.0/max(1,total_counts[i]) for i in range(len(mutations))]))>15:
							#print 'Coverage', min(counts) < count_filter
							#print 'Mutation', np.mean(observed) < mut_filter
							valid_region = False
							filter_discard += 1

				if valid_region:
					gini = find_gini(observed)
					gc_percent = np.sum(gc)*1.0/len(gc)
					if (gc_percent - gini) <= -0.2 :
						# pass
						print 'Weird', gene, start_pos, end_pos, observed, gini 
					if not gini:
						print 'Invalid Gini value at on reference {}. mutations: {}. counts: {}'.format(ref_dict[gene][start_pos:end_pos], mutations, counts)
					else:
						if coding:	
							coding_gini.append(gini)
							coding_gc.append(gc_percent)
						else:
							utr_gini.append(gini)
							utr_gc.append(gc_percent)
						try:
							good_genes[gene] += [(start_pos, end_pos)]
						except KeyError:
							good_genes[gene] = [(start_pos, end_pos)]

	return (coding_gc, coding_gini, utr_gc, utr_gini, filter_discard, bad_genes, min_coverage, good_genes)



def find_gini(T):
	#Where T is a sorted list of signal values! i.e mutation percentage
	sumy=0
	n = len(T)
	for count in T: 
		sumy = sumy + (n+1.0- (T.index(count)+1.0))*count
	sumAll = float(sumy)/sum(T)
	b1 = n+ 1.0 -2.0*sumAll
	giniv = b1/n
	return giniv 

def get_start_and_end(start_file):
	'''dictionary from gene to coding start and end position
		positions[chrName]: (start, stop)
	'''
	f = open(start_file)
	positions = {}
	reader=csv.reader(f,delimiter='\t')
	for line, start, stop in reader:
		#Get reference name
		chrNameini = line.strip().strip(">") # storing the chr name in a variable
		try:
			index= chrNameini.index("H")
		except ValueError:
			pass

		chrName = chrNameini[0:index-1]

		if "PREDICTED" in chrName:
			retract=chrName.index('PRED')
			chrName=chrName[0:retract-1]

		#Populate dictionary with start and stop positions
		try:
			positions[chrName] = (start, stop)
		except KeyError:
			print line 
	pickle.dump(positions, open('Gene_to_StartStop.p','wb'))
	return positions


def get_reference_sequence(ref_file):
	"""This is a more elegant version that uses the full file

	returs d[chromo,position] =base"""

	seq_for_chrName = ""

	chrName = ""

	mydict = {}

	f= open(ref_file)

	for line in f:
		if line.find(">") >= 0:
		# put the previous chr data, if there is any, in the dict.

			if (chrName != ""):
				mydict[chrName] = seq_for_chrName
			# Time for new chr and seq.
			chrNameini = line.strip().strip(">") # storing the chr name in a variable
			try:
				index= chrNameini.index("H")
			except ValueError:
				pass
			chrName = chrNameini[0:index-1]

			if "PREDICTED" in chrName:
				retract=chrName.index('PRED')
				chrName=chrName[0:retract-1]
			seq_for_chrName = "" # clearing the variable here.

		elif line.find("<") == -1 :
		# concatenate all strings that don't have > into a string.
			seq_for_chrName = seq_for_chrName + line.strip()#"\n")
	pickle.dump(mydict, open('ChrPos_to_Base.p', 'wb'))
	return mydict 

if __name__ == '__main__':                      
	main()
