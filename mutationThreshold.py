# Constructs a plot of the fraction of mRNAs in the dataset which 
# have an average mismatch coverage over a certain value, for 
# different amount of reads -50, 100, 200 million etc until we achieve >25% 
# fraction of mRNAs

# General Procedure
# 1. Sample some number of reads, N
# 2. Count the average number of mismatches per gene 
# 3. Determine the fraction of reads which exceed a certain value
# 4. Repeat for higher values of N

# Import relevant packages
import matplotlib 
import pysam as py
matplotlib.use('agg') 
import numpy as np
import scipy
from scipy import stats
import cPickle as pickle
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import re

def main():
	# Define variables
	
	inc = 100		#Size of rolling windows over genes
	avg_mismatch_cvg = [i for i in xrange(0,40,1)]+[i for i in xrange(40,155,5)]
	gene_fractions = []
	use_dict = True	#Use pre-existing list of mutations
	plot_simulated = True 	# Plot 1 billion simulated reads


	# Sam alignment files
	sam_filename = 'MegZ37_38_sorted.bam'
	samfile = py.AlignmentFile(sam_filename)
	num_reads = int(samfile.mapped/1e06)
	print 'Number of reads:', num_reads
	N_vals = [50, 100, 198]
	#N_vals = [i for i in xrange(50, num_reads, 25)]
	if use_dict:		
		try:
			gene_fractions = pickle.load(open( 'gene_fracs.p', "rb" ))
			print 'Got fractions...', gene_fractions
		except getattr(__builtins__,'FileNotFoundError', IOError):
			print 'Getting fractions...'
			for N in N_vals:
				frac_genes = get_fracs(samfile, N, avg_mismatch_cvg, inc)
				print 'frac_genes', frac_genes
				print 'Got a fraction for N = {} ...'.format(N)
				gene_fractions.append(frac_genes)
			pickle.dump( gene_fractions, open( 'gene_fracs.p', "wb" ))
			print 'Got fractions...', gene_fractions
	else:
		print 'Getting fractions...'
		for N in N_vals:
			frac_genes = get_fracs(samfile, N, avg_mismatch_cvg, inc)
			print 'frac_genes', frac_genes
			print 'Got a fraction for N = {} ...'.format(N)
			gene_fractions.append(frac_genes)
		pickle.dump( gene_fractions, open( 'gene_fracs.p', "wb" ))
		print 'Got fractions...', gene_fractions
	print 'Plotting graph'
	color = ['r', 'g', 'b', 'y']
	yis15 = np.array([15]*len(gene_fractions[0]))
	gene_fractions = np.array(gene_fractions)
	avg_mismatch_cvg = np.array(avg_mismatch_cvg)
	#gene_fractions = [[17409,1123,442,319,247,200,173,145,127,110,101,93,86,83,78,74,67,66,63,59,58,54,49,47,45,43,41,39,39,39,38],[17409,1123,442,319,247,200,173,145,127,110,101,93,86,83,78,74,67,66,63,59,58,54,49,47,45,43,41,39,39,39,38]]
	lines = []
	for i in range(len(N_vals)):
		lines.append(plt.plot(avg_mismatch_cvg, gene_fractions[i], color = color[i], label=N_vals[i]))
		plt.annotate('{}'.format(str(gene_fractions[i][20])[0:5]), xy=(20,gene_fractions[i][20]+0.001), fontsize = 4)	
	
	if plot_simulated:
		N_vals.append(1000)
		print 'Simulating reads...'
		#Plot 1 billion simulated reads
		try:
			sim_fractions = pickle.load(open( 'gene_fracs_5_simulated.p', "rb" ))
		except getattr(__builtins__,'FileNotFoundError', IOError):
			sim_fractions = simulate_reads(samfile, 198, avg_mismatch_cvg, inc)
			pickle.dump( sim_fractions, open( 'gene_fracs_5_simulated.p', "wb" ))
		print '...Done simulation.'
		print sim_fractions
		lines.append(plt.plot(avg_mismatch_cvg, sim_fractions, color = 'teal', label=1000))
		plt.annotate('{}'.format(str(sim_fractions[20])[0:5]), xy=(20,sim_fractions[20]+0.001), fontsize = 4)		

	plt.axvline(x=10, color='khaki')	#Plot a line at coverage 10
	plt.axvline(x=15, color='gold')	#Plot a line at coverage 15
	plt.axvline(x=20, color='goldenrod')	#Plot a line at coverage 20
	
	plt.grid(False)
	plt.ylabel('Fraction of Genes')
	plt.xlabel('Average Mismatch Coverage')
	plt.title('Fraction of Reads Exceeding Average Mismatch Coverage')
	plt.legend(['N = {} million'.format(N_vals[i]) for i in range(len(N_vals))])
	plt.savefig('CoveragePlotEPS.eps', format='eps', dpi=1000)

	return None 

def simulate_reads(samfile, N, avg_mismatch_cvg, inc):
	'''
	Given a samfile considers reads 10 times each until N reads are obtained
	'''
	return get_fracs(samfile, N, avg_mismatch_cvg, inc, 5)

def get_fracs(samfile, N, avg_mismatch_cvg, inc, wt = 1.0):
	'''
	Given a sam file and a number of reads to consider,
	returns the fraction of genes whose average mismatch 
	coverage exceed given values
	'''
	

	max_mismatch = []	#Max mismatch average over genes 
	fraction_genes = []

	#For N reads, calculate the number of mutations per position on each gene
	pos_mutations, neg_mutations, skipped_reads, skipped_multiple = count_mutations(samfile, N, wt)


	#For each gene, calculate the best window average
	for k,gene in enumerate(samfile.references):
		pos_d = pos_mutations[gene]['mutations_counts']
		length = samfile.lengths[k]
		windows = [(i, i+inc) for i in range(0, length ,inc)]

		mut_avg = 0
		for (start,end) in windows:
			if end > length:
				end = length
			avg = np.mean(np.array(pos_d[start:end]))
			if avg > mut_avg:
				mut_avg = avg 
		max_mismatch.append(mut_avg)

	#Determine the fraction of genes whose max window average exceeds a certain threshold
	max_mismatch = np.array(max_mismatch)
	num_genes = len(max_mismatch)
	for mismatch_cvg in avg_mismatch_cvg:
		#print 'hey', max_mismatch, mismatch_cvg
		good_regions = np.where(max_mismatch >= mismatch_cvg)[0]
		#print 'they', mismatch_cvg, len(good_regions)*1.0/num_genes
		fraction_genes.append(len(good_regions)*1.0/num_genes)

	return fraction_genes

def count_mutations(samfile, N, wt):
	'''
	Count mutations per base along the genes in the file sam_filename
	for up to N reads.
	'''

	N = N*int(1e6) 

	unmpapped = 0 

	# Create the dictionary containing the mutation counts for each reference sequence
	mutations_positive = {}
	mutations_negative = {}
	mds = {}        #A dictionary containing the MD tag for each reference read name.

	for i,ref in enumerate(samfile.references):    #(i, reference[i])
		mutations_positive[ref] = {'mutations_counts':[0]*(samfile.lengths[i]+2)}
		mutations_negative[ref] = {'mutations_counts':[0]*(samfile.lengths[i]+2)}

		mds[ref] = []

	# Iterate through the reads and detect the mutations
	# Keep track of discarded reads, due to multiple alignments or insertions/deletions etc
	skipped_reads = 0
	skipped_multiple = 0
	num_check = 0
	print 'Starting at i = 0, with N = ', N
	for read in samfile:
		if read.is_unmapped:
			unmpapped += 1
			continue 
		if num_check >= N:
			print 'Got to {} reads, cutting off'.format(num_check)
			break
		else:
			num_check += 1
			# Only consider reads without insertion ('I'), deletion ('D'), skipped region from the ref('N') or padding ('P')
			# Only consider reads which are single mapped
			multiple_mapped = False
			try:
				if read.get_tag('AS') == read.get_tag('XS'):
					multiple_mapped = True
					skipped_multiple += 1
			except KeyError:
				pass 

			if len(set(read.cigarstring).intersection(set('IDNP'))) != 0 or multiple_mapped: 
				skipped_reads +=1      #Keep count of lost reads
			else:
				#The read contains no insertions, deletions, unknowns or P, and does not have multiple similar alignments

				md = read.get_tag('MD')             #read the MD tag from the file
				md_matches = re.findall(r'\d+',md)  #this means get all the numbers
				md_muts = re.findall(r'\D+',md)     #this means get everything but the numbers
				mut_start = read.reference_start    #here we are putting it already in terms of ref start not read start 

				intoRead = 0                        #keeps track oh how many bases into the read our mutation is
		 
				for j,mut in enumerate(md_muts):    #i, mut[i]
					
					mut_start += int(md_matches[j]) #mutation position relative to reference start 
					mut = re.split('\^',mut)
					mut = ''.join(mut)              #Get rid of spurious ^ characters. e.g. T^CC becomes TCC 
					mut_end = mut_start + len(mut)

					#Add mutations to dictionary of positions vs count, depending on whether the read is the negative strand or the positive strand
					intoRead+= int(md_matches[j])
					for i in xrange(mut_start,mut_end):
						
						if intoRead >=6:            #only update the mutation counts if the initial match was more than 6 bases 
							
							if read.is_reverse:
								mutations_negative[read.reference_name]['mutations_counts'][i] += wt
								
							else:
								mutations_positive[read.reference_name]['mutations_counts'][i] += wt

					intoRead += len(mut)
					mut_start = mut_end

				#print read.reference_name, 'ref end +1', read.reference_end + 1
				mds[read.reference_name].append(md)
	print 'Considered {} reads'.format(num_check)	 
	print 'unmpapped: ', unmpapped
	samfile.reset()
	return(mutations_positive, mutations_negative, skipped_reads, skipped_multiple)


if __name__ == '__main__':                      
	main()
