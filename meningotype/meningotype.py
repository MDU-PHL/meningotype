#!/usr/bin/env python
# Script by Jason Kwong
# In silico typing for Neisseria meningitidis

# Use modern print function from python 3.x
from __future__ import print_function

# Usage
import argparse
from argparse import RawTextHelpFormatter
import sys
import os
import os.path
import StringIO
import urllib
import subprocess
from subprocess import Popen
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastxCommandline
from pkg_resources import resource_string, resource_filename

###### Script globals ##########################################################

# URLs to update database files
porA1URL = 'http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_neisseria_seqdef&page=downloadAlleles&locus=PorA_VR1'
porA2URL = 'http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_neisseria_seqdef&page=downloadAlleles&locus=PorA_VR2'
fetAURL = 'http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_neisseria_seqdef&page=downloadAlleles&locus=FetA_VR'

# allele sizes and serotype dictionary
alleleSIZE = {'A':92, 'B':169, 'C':74, 'W':129, 'X':65, 'Y':146}
seroDICT = {'sacB':'A', 'synD':'B', 'synE':'C', 'synG':'W', 'xcbB':'X', 'synF':'Y'}

porASEQS = []
fetASEQS = []
sero = None
porA = None
fet = None

# field separator
sep = '\t'

# Standard functions
# Log a message to stderr
def msg(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

# Log an error to stderr and quit with non-zero error code
def err(*args, **kwargs):
	msg(*args, **kwargs)
	sys.exit(1);

# Format database
def format(file):
	sed_all(file, '^.*<textarea name="concatenation".*?>', '')
	sed_all(file, '<\/textarea>.*$', '')
	sed_inplace(file, '&gt;', '>')

############### Database functions #############################################

# Update database files
def update_db(db_file, db_url):
	if os.path.isfile(db_file):
		os.rename(db_file, db_file+'.old')
	urllib.urlretrieve(db_url, db_file)

# Check database files are present
def check_primer_files(f):
	if not os.path.isfile(f):
		err('ERROR: Cannot locate file: "{}"'.format(f))

def check_db_files(f, db_url):
	if not os.path.isfile(f):
		update_db(f, db_url)

# Set up BLAST databases if not present
def makeblastDB(db, infile, dbtype):
	if dbtype == 'nucl':
		DBindex = db + '.nin'
	elif dbtype == 'prot':
		DBindex = db + '.pin'
	if not os.path.isfile(DBindex):
		proc = subprocess.Popen(['makeblastdb', '-in', infile, '-out', db, '-dbtype', dbtype], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

############### Main serotyping functions ######################################
def seroTYPE(f, seroprimers, allelesdb):
	seroCOUNT = []				# Setup list in case there are mixed/multiple hits
	proc = subprocess.Popen(['isPcr', f, seroprimers, 'stdout', '-minPerfect=10'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	PCRout = proc.communicate()[0]
	if not PCRout:
		sero = None
		seroBLAST = NcbiblastnCommandline(query=f, db=allelesdb, task='blastn', perc_identity=90, evalue='1e-20', outfmt='"6 sseqid pident length"')
		stdout, stderr = seroBLAST()
		lenMATCH = 0
		line = stdout.split('\n')[0]
		amp = line.split('\t')
#		ampID = amp[0]
#		ampPERC = amp[1]
#		ampLEN = int(amp[2])
		sero = amp[0]		# Currently only takes top/first BLAST hit
		if not sero:
			seroCOUNT.append('-')
		else:
			seroCOUNT.append(sero)
	else:
		alleleSEQ = StringIO.StringIO()
		alleleSEQ.write(PCRout)
		alleleSEQ.seek(0)
		for amplicon in SeqIO.parse(alleleSEQ, "fasta"):
			product = amplicon.description.split()
			ampID = product[1]
			ampLEN = int(product[2][:-2])
			sero = seroDICT[ampID]
			expLEN = int(alleleSIZE[sero])
			if ampLEN > (expLEN-6) and ampLEN < (expLEN+6):
				seroCOUNT.append(sero)
		alleleSEQ.close()
	return seroCOUNT

def finetypeBLAST(s, db):
	ft = None
	allele = None
	ftBLAST = NcbiblastxCommandline(query='-', db=db, outfmt='"6 qseqid sseqid pident length gaps"', seg='no')
	stdout, stderr = ftBLAST(stdin=s.format('fasta'))
	if stdout:
		BLASTout = stdout.split('\n')
		lenMATCH = 0
		for line in BLASTout:
			if line.strip():
				BLASTline = line.split('\t')
				if float(BLASTline[2]) == 100 and BLASTline[4] == '0':
					if int(BLASTline[3]) > int(lenMATCH):
						lenMATCH = BLASTline[3]
						ftRESULT = BLASTline[1]
						ft = ftRESULT.split('_')[2]
				elif not ft:
					ft = 'new'
	return str(ft)

def fineTYPE(f, finetypeprimers, pora1db, pora2db, fetdb):
	porACOUNT = []				# Setup list in case there are mixed/multiple hits
	fetACOUNT = []
	global porASEQS
	global fetASEQS
	proc = subprocess.Popen(['isPcr', f, finetypeprimers, 'stdout', '-maxSize=800', '-tileSize=10', '-minPerfect=8', '-stepSize=3'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	PCRout = proc.communicate()[0]
	alleleSEQ = StringIO.StringIO()
	alleleSEQ.write(PCRout)
	alleleSEQ.seek(0)
	for amplicon in SeqIO.parse(alleleSEQ, "fasta"):
		ampFILE = amplicon.description.split()
		ampID = ampFILE[1]
		ampLEN = ampFILE[2]					# Need to check amplicon length to exclude double hits?
		if ampID == 'porA':
			A1 = finetypeBLAST(amplicon, pora1db)
			A2 = finetypeBLAST(amplicon, pora2db)
			porA = A1 + ',' + A2
			porACOUNT.append(porA)
			porASEQ = amplicon.seq.upper()
			porARECR = SeqRecord(porASEQ, id=f, description='PorA')
			porASEQS.append(porARECR)
		if ampID == 'fetA':
			fet = finetypeBLAST(amplicon, fetdb)
			fetACOUNT.append(fet)
			fetASEQ = amplicon.seq.upper()
			fetARECR = SeqRecord(fetASEQ, id=f, description='FetA')
			fetASEQS.append(fetARECR)
	if len(porACOUNT) == 0:
		porACOUNT.append('-')
	if len(fetACOUNT) == 0:
		fetACOUNT.append('-')
	return porACOUNT, fetACOUNT

def results(f, seroCOUNT, porCOUNT, fetCOUNT):
	seroTYPE = '/'.join(seroCOUNT)
	if not porCOUNT:
		print(f + sep + seroTYPE)
	else:
		porTYPE = '/'.join(porCOUNT)
		fetTYPE = '/'.join(fetCOUNT)
		print(f + sep + seroTYPE + sep + porTYPE + sep + fetTYPE)

########## Meningotype main ####################################################

def main():
	parser = argparse.ArgumentParser(
		formatter_class=RawTextHelpFormatter,
		description='In silico typing for Neisseria meningitidis\n'
			'\nPCR Serotyping Ref: Mothershed et al, J Clin Microbiol 2004; 42(1): 320-328\n'
			'\nporA and fetA typing Ref: Jolley et al, FEMS Microbiol Rev 2007; 31: 89-96\n'
			'See also http://www.neisseria.org/nm/typing/',
		usage='\n  %(prog)s [OPTIONS] <fasta1> <fasta2> <fasta3> ... <fastaN>')
	parser.add_argument('fasta', metavar='FASTA', nargs='*', help='input FASTA files eg. fasta1, fasta2, fasta3 ... fastaN')
	# CSV option excluded due to syntax of porA finetype VR1,VR2
	#parser.add_argument('--csv', action='store_true', default=False, help='output comma-separated format (CSV) rather than tab-separated')
	parser.add_argument('--finetype', action='store_true', help='Perform porA and fetA fine typing (default=off)')
	parser.add_argument('--db', metavar='DB', help='Specify custom directory containing allele databases for porA/fetA typing\n'
		'Directory must contain database files "FetA_VR.fas", "PorA_VR1.fas", and "PorA_VR2.fas"')
	parser.add_argument('--printseq', action='store_true', help='Save porA, porB and fetA allele sequences to file (default=off)')
	parser.add_argument('--updatedb', action='store_true', default=False, help='update allele database from <www.ng-mast.net>')
	parser.add_argument('--test', action='store_true', default=False, help='run test example')
	parser.add_argument('--version', action='version', version=
		'=====================================\n'
		'%(prog)s v0.1-beta\n'
		'Updated 4-Jul-2016 by Jason Kwong\n'
		'Dependencies: isPcr, BLAST, BioPython\n'
		'=====================================')
	args = parser.parse_args()

	if args.db:
		DBpath = str(args.db).rstrip('/')
	else:
		DBpath = resource_filename(__name__, 'db')
	# Path to database files
	porA1alleles = os.path.join( DBpath, 'PorA_VR1.fas' )
	porA2alleles = os.path.join( DBpath, 'PorA_VR2.fas' )
	fetalleles = os.path.join( DBpath, 'FetA_VR.fas' )
	seroALLELES = os.path.join( DBpath, 'seroALLELES.fa' )
	allelesDB = os.path.join( DBpath, 'blast', 'seroALLELES' )
	seroPRIMERS = os.path.join( DBpath, 'seroPRIMERS' )
	finetypePRIMERS = os.path.join( DBpath, 'finetypePRIMERS' )
	porA1DB = os.path.join( DBpath, 'blast', 'porA1' )
	porA2DB = os.path.join( DBpath, 'blast', 'porA2' )
	fetDB = os.path.join( DBpath, 'blast', 'fet' )

	if args.updatedb:
		msg('Updating "{}" ... '.format(porA1alleles))
		update_db(porA1alleles, porA1URL)
		msg('Updating "{}" ... '.format(porA2alleles))
		update_db(porA2alleles, porA2URL)
		msg('Updating "{}" ... '.format(fetalleles))
		update_db(fetalleles, fetAURL)
		msg('Done.')
		sys.exit(0)


	if not os.path.exists(DBpath):
		err('ERROR: Cannot locate "db" directory at "{}"'.format(DBpath))
	check_primer_files(seroPRIMERS)
	check_primer_files(seroALLELES)
	check_primer_files(finetypePRIMERS)
	check_db_files(porA1alleles, porA1URL)
	check_db_files(porA2alleles, porA2URL)
	check_db_files(fetalleles, fetAURL)


	makeblastDB(allelesDB, seroALLELES, 'nucl')
	makeblastDB(porA1DB, porA1alleles, 'prot')
	makeblastDB(porA2DB, porA2alleles, 'prot')
	makeblastDB(fetDB, fetalleles, 'prot')

	# Test example to check meningotype works
	if args.test:
		TESTpath = resource_filename(__name__, 'test')
		testSEQS = [os.path.join( TESTpath, f ) for f in ['A.fna', 'B.fna', 'C.fna', 'W.fna', 'X.fna', 'Y.fna'] ]
		msg('\033[94mRunning meningotype.py on test examples ... \033[0m')
		msg('$ meningotype.py A.fna B.fna C.fna W.fna X.fna Y.fna')
		args.fasta = testSEQS

	# Create folder for output sequences if specified
	if args.printseq:
		try:
			if not os.path.exists('printseq'):
			    os.makedirs('printseq')
			else:
				err('"printseq" folder already exists in this directory.')
		except:
			err('Unable to create "printseq" folder in this directory.')

	# Run meningotype
	if len(args.fasta) == 0:
		message = "\033[91mEither use --test or specify at least one FASTA file.\033[0m"
		sys.stderr.write('error: {}\n'.format( message ) )
		parser.print_help()
		parser.exit(1)
	print('SAMPLE_ID' + sep + 'SEROTYPE' + sep + 'PorA_TYPE' + sep + 'FetA_TYPE')
	for f in args.fasta:
		seroCOUNT = seroTYPE(f, seroPRIMERS, allelesDB)
		if not args.finetype:
			results(f, seroCOUNT, None, None)
		else:
			ftRESULTS = fineTYPE(f, finetypePRIMERS, porA1DB, porA2DB, fetDB)
			porACOUNT = ftRESULTS[0]
			fetACOUNT = ftRESULTS[1]
			results(f, seroCOUNT, porACOUNT, fetACOUNT)

	# Print allele sequences to file
	if args.printseq:
		with open('printseq/porA_seqs.fasta', 'w') as output:
			SeqIO.write(porASEQS, output, 'fasta')
		with open('printseq/fetA_seqs.fasta', 'w') as output:
			SeqIO.write(fetASEQS, output, 'fasta')

	sys.exit(0)

if __name__ == "__main__":
	main()
