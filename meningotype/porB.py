#!/usr/bin/env python
# Script by Jason Kwong
# Extracts porB sequence from Neisseria meningitidis

# Use modern print function from python 3.x
from __future__ import print_function

# Modules
import argparse
from argparse import RawTextHelpFormatter
import sys
import os
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from pkg_resources import resource_string, resource_filename

# Import local modules
from . import nmen

# Standard functions
# Log a message to stderr
def msg(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

# Log an error to stderr and quit with non-zero error code
def err(*args, **kwargs):
	msg(*args, **kwargs)
	sys.exit(1);

# BLAST
def porBBLAST(f, blastdb):
	porB = 'new'
	blast_qseqid = '-'
	blast_pident = '-'
	blast_cov = '<99'
	porBRECR = None
	fBLAST = NcbiblastnCommandline(query=f, db=blastdb, outfmt="'6 qseqid sseqid pident length sstrand qstart qend sstart send slen'", dust='no', culling_limit=1)
	stdout, stderr = fBLAST()
	blastOUT = stdout.split('\t')
	if len(blastOUT) == 10:
		blast_qseqid = blastOUT[0]
		blast_sseqid = blastOUT[1]
		blast_pident = float(blastOUT[2])
		blast_length = int(blastOUT[3])
		blast_sstrand = blastOUT[4]
		blast_qstart = int(blastOUT[5])
		blast_qend = int(blastOUT[6])
		blast_sstart = int(blastOUT[7])
		blast_send = int(blastOUT[8])
		blast_slen = int(blastOUT[9])
		blast_cov = float(blast_length)/float(blast_slen)*100
		if blast_cov > 99:
			for s in SeqIO.parse(f, 'fasta'):
				if s.id == blast_qseqid:
					blastCONTIG = s.seq
			if blast_sstrand == 'plus':
				start = blast_qstart - blast_sstart
				end = start + blast_slen
				porBSEQ = blastCONTIG[start:end]
			else:
				end = blast_qend + blast_send - 1
				start = end - blast_slen
				porBSEQ = blastCONTIG[start:end].reverse_complement()
			porBRECR = SeqRecord(porBSEQ, id=f, description='PorB')
			if blast_cov == 100 and blast_pident == 100:
				porB = blast_sseqid
			elif blast_cov > 99:
				porB = ''.join([blast_sseqid, '-like'])
	result = [f, blast_qseqid, porB, str(blast_pident), str(blast_cov), porBRECR]
	return result

def main():
	# Usage
	parser = argparse.ArgumentParser(
		formatter_class=RawTextHelpFormatter,
		description='PorB typing of Neisseria meningitidis',
		usage='\n  %(prog)s FASTA-1 FASTA-2 ... FASTA-N')
	parser.add_argument('fasta', metavar='FASTA', nargs='+', help='FASTA file to search (required)')
	parser.add_argument('--db', metavar='DB', help='specify custom directory containing allele databases')
	parser.add_argument('--printseq', action='store_true', help='save porB allele sequences to file (default=off)')
	parser.add_argument('--version', action='version', version=
		'=====================================\n'
		'%(prog)s v0.1\n'
		'Updated 22-Feb-2017 by Jason Kwong\n'
		'Dependencies: Python 2.x, BioPython, BLAST\n'
		'=====================================')
	args = parser.parse_args()
	
	if args.db:
		DBpath = str(args.db).rstrip('/')
	else:
		DBpath = resource_filename(__name__, 'db')

	porBDB = os.path.join( DBpath, 'blast', 'porB' )

	# Main
	porBSEQS = []
	print('\t'.join(['SAMPLE_ID', 'CONTIG', 'PorB', '%ID', 'COV']))
	for f in args.fasta:
		result = porBBLAST(f, porBDB)
		print('\t'.join(result[:-1]))
		porBSEQS.append(result[5])

	# Print allele sequences to file
	if args.printseq:
		if porBSEQS:
			with open('printseq/porB_seqs.fasta', 'w') as output:
				SeqIO.write(porBSEQS, output, 'fasta')

if __name__ == "__main__":
	main()