#!/usr/bin/env python
# Script by Jason Kwong
# Checks for ctrA gene in Neisseria meningitidis

# Use modern print function from python 3.x
from __future__ import print_function

# Modules
import argparse
from argparse import RawTextHelpFormatter
import sys
import os
import subprocess
from subprocess import Popen
from io import StringIO
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from pkg_resources import resource_string, resource_filename

# Standard functions
# Log a message to stderr
def msg(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

# Log an error to stderr and quit with non-zero error code
def err(*args, **kwargs):
	msg(*args, **kwargs)
	sys.exit(1);

# ctrA isPCR
def ctrA_PCR(f, p, dbpath):
	ctrAPRIMERS = os.path.join( dbpath, 'ctrAPRIMERS' )
	ctrADB = os.path.join( dbpath, 'blast', 'ctrA' )
	resultBLAST = '-'
	resultPCR = '-'
	proc = subprocess.Popen(['isPcr', f, ctrAPRIMERS, 'stdout', '-minPerfect=10'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	PCRout = proc.communicate()[0].decode('UTF-8')
	if not PCRout:
		ctrABLAST = NcbiblastnCommandline(query=f, db=ctrADB, task='blastn', perc_identity=90, evalue='1e-20', outfmt='"6 sseqid pident length"', culling_limit='1')
		
		stdout, stderr = ctrABLAST()
		msg(stdout)
		if stdout:
			lenMATCH = 0
			line = stdout.split('\n')[0]
			amp = line.split('\t')
			resultBLAST = amp[1]		# Currently only takes top/first BLAST hit
	else:
		alleleSEQ = StringIO()
		alleleSEQ.write(PCRout)
		alleleSEQ.seek(0)
		for amplicon in SeqIO.parse(alleleSEQ, "fasta"):
			
			product = amplicon.description.split()
			msg(product)
			ampID = product[1]
			ampLEN = int(product[2][:-2])
			if ampLEN > 100 and ampLEN < 120:
				resultPCR = 'ctrA'
				resultBLAST = 'ctrA'
	result = '\t'.join([f, resultPCR, resultBLAST])
	if p:
		print(result)
	else:
		return result

def main():
	# Usage
	parser = argparse.ArgumentParser(
		formatter_class=RawTextHelpFormatter,
		description='Checks for ctrA gene in Neisseria meningitidis',
		usage='\n  %(prog)s FASTA-1 FASTA-2 ... FASTA-N')
	parser.add_argument('fasta', metavar='FASTA', nargs='+', help='FASTA file to search (required)')
	parser.add_argument('--db', metavar='DB', help='specify custom directory containing allele databases')
	parser.add_argument('--version', action='version', version=
		'=====================================\n'
		'%(prog)s v0.1\n'
		'Updated 22-Feb-2017 by Jason Kwong\n'
		'Dependencies: Python 2.x, BioPython, BLAST\n'
		'=====================================')
	args = parser.parse_args()
	
	DBpath = resource_filename(__name__, 'db')
	
	# Main
	print('\t'.join(['SAMPLE_ID', 'PCRresult', 'BLASTresult']))
	for f in args.fasta:
		ctrA_PCR(f, True, DBpath)
	
if __name__ == "__main__":
	main()