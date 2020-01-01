#!/usr/bin/python3
__author__ = "Houcemeddine Othman"
__credits__ = "Paractical, reproducible research PHINDaccess 2020"
__maintainer__ = "Houcemeddine Othman"
__email__ = "houcemoo@gmail.com"

import matplotlib.pylab as plt 
import gzip
import sys

# check if you have python 3
if sys.version_info[0] < 3:
    raise Exception("You are not using python3")

class read_input:
	"""docstring for prepare"""
	def __init__(self, fastafile, bedfile):
		with gzip.open( fastafile, 'rt') as inputfasta:
			self.fasta =  inputfasta.readlines()

		with open(bedfile, 'r' ) as inputbed: 
			self.bed = inputbed.readlines()

	def prepare_seq(self):	
		"""Will join all the lines of the sequence in one string"""		
		self.sequence = "".join(line.strip() for line in self.fasta[1:])

class compute(read_input): 
	""" The class for the computing stage """
	def __init__(self, fastafile, bedfile):
		super().__init__(fastafile, bedfile )

	def computeGC(self):
		GC_content_pergene = []
		not_gene_coor = []		
		for line in self.bed: 
			start = int(line.split()[1]) -1
			end = int(line.split()[2]) - 1
			slice_gene = (self.sequence[start:end])
			geneLen = end-start
			GC =  occurenceGC(slice_gene) 
			try :
				GC_content_pergene.append(GC/geneLen)
			except:
				pass
		return GC_content_pergene

def occurenceGC(string):
	""" Calculate the number of GCs in a string """
	G_number = string.count('G')
	C_number = string.count('C')
	return G_number+C_number

def plot_histogram(list_of_GCs, nbins = 50, path2output="./" ):
	""" Plots the figure and generate PNG and SVG formats"""
	plt.hist(list_of_GCs, bins=nbins, histtype="barstacked", alpha = 0.5, edgecolor = "k")
	plt.ylabel('Probability')
	plt.xlabel('GC %')
	plt.savefig(path2output+"GC.png")
	plt.savefig(path2output+"GC.svg")

if __name__ == "__main__":
	myGC = compute("Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz", "data_HX/genes_x.bed")
	myGC.prepare_seq()
	print("Calculating GC content...")
	GCs = myGC.computeGC()
	print("FINISH")
	print("plotting to GC.png and GC.svg ...")
	plot_histogram(GCs)
	print("FINISH")


