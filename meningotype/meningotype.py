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

# Import local modules
import scripts.nmen as nmen
import scripts.menwy as menwy
import scripts.ctrA as ctrA

###### Script globals ##########################################################

# URLs to update database files (OLD)
# porA1URL = 'http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_neisseria_seqdef&page=downloadAlleles&locus=PorA_VR1'
# porA2URL = 'http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_neisseria_seqdef&page=downloadAlleles&locus=PorA_VR2'
# fetAURL = 'http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_neisseria_seqdef&page=downloadAlleles&locus=FetA_VR'
# fHbpURL = 'http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_neisseria_seqdef&page=downloadAlleles&locus=fHbp_peptide'
# NHBAURL = 'http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_neisseria_seqdef&page=downloadAlleles&locus=NHBA_peptide'
# NadAURL = 'http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_neisseria_seqdef&page=downloadAlleles&locus=NadA_peptide'
# BASTURL = 'http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_neisseria_seqdef'

# URLs to update database files (REST API)
porA1URL = 'http://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/PorA_VR1/alleles_fasta'
porA2URL = 'http://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/PorA_VR2/alleles_fasta'
fetAURL = 'http://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/FetA_VR/alleles_fasta'
fHbpURL = 'http://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/fHbp_peptide/alleles_fasta'
NHBAURL = 'http://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/NHBA_peptide/alleles_fasta'
NadAURL = 'http://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/NadA_peptide/alleles_fasta'
BASTURL = 'http://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/53/profiles_csv'

# allele sizes and serotype dictionary
alleleSIZE = {'A':92, 'B':169, 'C':74, 'W':129, 'X':65, 'Y':146}
seroDICT = {'sacB':'A', 'synD':'B', 'synE':'C', 'synG':'W', 'xcbB':'X', 'synF':'Y'}

porASEQS = []
fetASEQS = []
fHbpSEQS = []
NHBASEQS = []
NadASEQS = []
sero = None
porA = None
fet = None
fHbp = None
NHBA = None
NadA = None

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
		seroBLAST = NcbiblastnCommandline(query=f, db=allelesdb, task='blastn', perc_identity=90, evalue='1e-20', outfmt='"6 sseqid pident length"', culling_limit='1')
		stdout, stderr = seroBLAST()
		lenMATCH = 0
		line = stdout.split('\n')[0]
		amp = line.split('\t')
		sero = amp[0]		# Currently only takes top/first BLAST hit
		if not sero:
			seroCOUNT.append('-')
		else:
			if sero == 'W' or sero == 'Y':
				sero = seroWY(f, sero)
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
				if sero == 'W' or sero == 'Y':
					sero = seroWY(f, sero)
				seroCOUNT.append(sero)
		alleleSEQ.close()
	return seroCOUNT

def seroWY(f, sero):
	wyTYPE = menwy.menwy(f, False)
	wy = wyTYPE.split('\t')[1]
	if wy == '-':
		return sero
	else:
		return wy

def finetypeBLAST(s, db):
	ft = None
	allele = None
	ftBLAST = NcbiblastxCommandline(query='-', db=db, outfmt='"6 qseqid sseqid pident length slen gaps evalue"', seg='no', query_gencode='11')
	stdout, stderr = ftBLAST(stdin=s.format('fasta'))
	if stdout:
		BLASTout = stdout.split('\n')
		lenMATCH = 0
		for line in BLASTout:
			if line.strip():
				BLASTline = line.split('\t')
				if float(BLASTline[2]) == 100 and BLASTline[5] == '0':		# check 100% match with no gaps
					if int(BLASTline[3]) == int(BLASTline[4]):				# check length of db subject = length of match
						if int(BLASTline[3]) > int(lenMATCH):				# select longest match
							lenMATCH = BLASTline[3]
							ftRESULT = BLASTline[1]
							ft = ftRESULT.split('_')[2]
#				elif not ft:
#					ft = 'new'
		if not ft:															# if amplicon detected, but no match in db, assign as new
			ft = 'new'
	return str(ft)

def bxtypeBLAST(s, db):
	bx = None
	allele = None
	bxBLAST = NcbiblastxCommandline(query='-', db=db, outfmt='"6 qseqid sseqid pident length slen gaps evalue"', seg='no', culling_limit='1', evalue='1e-100', query_gencode='11')
	stdout, stderr = bxBLAST(stdin=s.format('fasta'))
	if stdout:
		BLASTout = stdout.split('\n')
		lenMATCH = 0
		for line in BLASTout:
			if line.strip():
				BLASTline = line.split('\t')
				if float(BLASTline[2]) == 100 and BLASTline[5] == '0':		# check 100% match with no gaps
					if int(BLASTline[3]) == int(BLASTline[4]):				# check length of db subject = length of match
						if int(BLASTline[3]) > int(lenMATCH):				# select longest match
							lenMATCH = BLASTline[3]
							bxRESULT = BLASTline[1]
							bx = bxRESULT.split('_')[2]
#				elif not bx:
#					bx = 'new'
		if not bx:															# if amplicon detected, but no match in db, assign as new
			bx = 'new'
	else:																	# if no amplicon detected, assign as 0
		bx = '0'
	return str(bx)

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

def bxTYPE(f, bxPRIMERS, fHbpDB, NHBADB, NadADB):
	fHbpCOUNT = []				# Setup list in case there are mixed/multiple hits
	NHBACOUNT = []
	NadACOUNT = []
	global fHbpSEQS
	global NHBASEQS
	global NadASEQS
	proc = subprocess.Popen(['isPcr', f, bxPRIMERS, 'stdout', '-maxSize=3000', '-tileSize=7', '-minPerfect=8', '-stepSize=3'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	PCRout = proc.communicate()[0]
	alleleSEQ = StringIO.StringIO()
	alleleSEQ.write(PCRout)
	alleleSEQ.seek(0)
	for amplicon in SeqIO.parse(alleleSEQ, "fasta"):
		ampFILE = amplicon.description.split()
		ampID = ampFILE[1]
		ampLEN = ampFILE[2]					# Need to check amplicon length to exclude double hits?
		if ampID == 'fHbp':
			fHbp = bxtypeBLAST(amplicon, fHbpDB)
			fHbpCOUNT.append(fHbp)
			fHbpSEQ = amplicon.seq.upper()
			fHbpRECR = SeqRecord(fHbpSEQ, id=f, description='fHbp')
			fHbpSEQS.append(fHbpRECR)
		if ampID == 'NHBA':
			NHBA = bxtypeBLAST(amplicon, NHBADB)
			NHBACOUNT.append(NHBA)
			NHBASEQ = amplicon.seq.upper()
			NHBARECR = SeqRecord(NHBASEQ, id=f, description='NHBA')
			NHBASEQS.append(NHBARECR)
		if ampID == 'NadA':
			NadA = bxtypeBLAST(amplicon, NadADB)
			NadACOUNT.append(NadA)
			NadASEQ = amplicon.seq.upper()
			NadARECR = SeqRecord(NadASEQ, id=f, description='NadA')
			NadASEQS.append(NadARECR)
	if len(fHbpCOUNT) == 0:
		fHbpCOUNT.append('0')
	if len(NHBACOUNT) == 0:
		NHBACOUNT.append('0')
	if len(NadACOUNT) == 0:
		NadACOUNT.append('0')
	return set(fHbpCOUNT), set(NHBACOUNT), set(NadACOUNT)

########## Meningotype main ####################################################

def main():
	parser = argparse.ArgumentParser(
		formatter_class=RawTextHelpFormatter,
		description='In silico typing for Neisseria meningitidis\n'
			'\nPCR Serotyping Ref: Mothershed et al, J Clin Microbiol 2004; 42(1): 320-328\n'
			'PorA and FetA typing Ref: Jolley et al, FEMS Microbiol Rev 2007; 31: 89-96\n'
			'Bexsero antigen sequence typing (BAST) Ref: Brehony et al, Vaccine 2016; 34(39): 4690-4697\n'
			'See also http://www.neisseria.org/nm/typing/',
		usage='\n  %(prog)s [OPTIONS] <fasta1> <fasta2> <fasta3> ... <fastaN>')
	parser.add_argument('fasta', metavar='FASTA', nargs='*', help='input FASTA files eg. fasta1, fasta2, fasta3 ... fastaN')
	# CSV option excluded due to syntax of porA finetype VR1,VR2
	#parser.add_argument('--csv', action='store_true', default=False, help='output comma-separated format (CSV) rather than tab-separated')
	parser.add_argument('--finetype', action='store_true', help='perform porA and fetA fine typing (default=off)')
	parser.add_argument('--bast', action='store_true', help='perform Bexsero antigen sequence typing (BAST) (default=off)')
	parser.add_argument('--db', metavar='DB', help='specify custom directory containing allele databases for porA/fetA typing\n'
		'directory must contain database files: "FetA_VR.fas", "PorA_VR1.fas", "PorA_VR2.fas"\n'
		'for Bexsero typing: "fHbp_peptide.fas", "NHBA_peptide.fas", "NadA_peptide.fas", "BASTalleles.txt"')
	parser.add_argument('--printseq', action='store_true', help='save porA/fetA or BAST allele sequences to file (default=off)')
	parser.add_argument('--updatedb', action='store_true', default=False, help='update allele database from <pubmlst.org>')
	parser.add_argument('--test', action='store_true', default=False, help='run test example')
	parser.add_argument('--version', action='version', version=
		'=====================================\n'
		'%(prog)s v0.6-beta\n'
		'Updated 22-Feb-2017 by Jason Kwong\n'
		'Dependencies: isPcr, BLAST+, BioPython\n'
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
	fHbpalleles = os.path.join( DBpath, 'fHbp_peptide.fas' )
	NHBAalleles = os.path.join( DBpath, 'NHBA_peptide.fas' )
	NadAalleles = os.path.join( DBpath, 'NadA_peptide.fas' )
	BASTalleles = os.path.join( DBpath, 'BASTalleles.txt' )

	seroALLELES = os.path.join( DBpath, 'seroALLELES.fa' )
	allelesDB = os.path.join( DBpath, 'blast', 'seroALLELES' )
	seroPRIMERS = os.path.join( DBpath, 'seroPRIMERS' )
	finetypePRIMERS = os.path.join( DBpath, 'finetypePRIMERS' )
	bxPRIMERS = os.path.join( DBpath, 'bexseroPRIMERS' )

	porA1DB = os.path.join( DBpath, 'blast', 'porA1' )
	porA2DB = os.path.join( DBpath, 'blast', 'porA2' )
	fetDB = os.path.join( DBpath, 'blast', 'fet' )
	fHbpDB = os.path.join( DBpath, 'blast', 'fHbp_peptide' )
	NHBADB = os.path.join( DBpath, 'blast', 'NHBA_peptide' )
	NadADB = os.path.join( DBpath, 'blast', 'NadA_peptide' )
 
	if args.updatedb:
		try:
			msg('Updating "{}" ... '.format(porA1alleles))
			update_db(porA1alleles, porA1URL)
			msg('Updating "{}" ... '.format(porA2alleles))
			update_db(porA2alleles, porA2URL)
			msg('Updating "{}" ... '.format(fetalleles))
			update_db(fetalleles, fetAURL)
			msg('Updating "{}" ... '.format(fHbpalleles))
			update_db(fHbpalleles, fHbpURL)
			msg('Updating "{}" ... '.format(NHBAalleles))
			update_db(NHBAalleles, NHBAURL)
			msg('Updating "{}" ... '.format(NadAalleles))
			update_db(NadAalleles, NadAURL)
			msg('Updating "{}" ... '.format(BASTalleles))
			update_db(BASTalleles, BASTURL)
			msg('Done.')
		except IOError:
			err('ERROR: Unable to update DB at "{}".\nCheck DB directory permissions and connection to http://pubmlst.org.'.format(DBpath))
		except HTTPError:
			err('ERROR: Unable to update DB at "{}". Check connection to http://pubmlst.org.'.format(DBpath))
		sys.exit(0)

	# Check paths and files
	if not os.path.exists(DBpath):
		err('ERROR: Cannot locate "db" directory at "{}"'.format(DBpath))
	check_primer_files(seroPRIMERS)
	check_primer_files(seroALLELES)
	check_primer_files(finetypePRIMERS)
	check_db_files(porA1alleles, porA1URL)
	check_db_files(porA2alleles, porA2URL)
	check_db_files(fetalleles, fetAURL)

	# Setup BLASTDB
	makeblastDB(allelesDB, seroALLELES, 'nucl')
	makeblastDB(porA1DB, porA1alleles, 'prot')
	makeblastDB(porA2DB, porA2alleles, 'prot')
	makeblastDB(fetDB, fetalleles, 'prot')

	if args.bast:
		check_primer_files(bxPRIMERS)
		check_db_files(fHbpalleles, fHbpURL)
		check_db_files(NHBAalleles, NHBAURL)
		check_db_files(NadAalleles, NadAURL)
		check_db_files(BASTalleles, BASTURL)
		makeblastDB(fHbpDB, fHbpalleles, 'prot')
		makeblastDB(NHBADB, NHBAalleles, 'prot')
		makeblastDB(NadADB, NadAalleles, 'prot')

		# Import allele profiles as dictionary
		BAST = {}
		dl = ','
		with open(BASTalleles) as db:
			for line in db:
				if line.strip():
					lines = line.split('\t')
					ST = lines[0]
					alleles = lines[1].rstrip() + dl + lines[2].rstrip() + dl + lines[3].rstrip() + dl + lines[4].rstrip() + dl + lines[5].rstrip('\n')
					BAST[alleles] = ST

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
				err('ERROR: "printseq" folder already exists in this directory.')
		except OSError:
			err('ERROR: Unable to create "printseq" folder in this directory.')

	# Run meningotype
	if len(args.fasta) == 0:
		message = "\033[91mEither use --test or specify at least one FASTA file.\033[0m"
		sys.stderr.write('error: {}\n'.format( message ) )
		parser.print_help()
		parser.exit(1)
	headers = ['SAMPLE_ID', 'SEROGROUP', 'ctrA', 'PorA', 'FetA', 'fHbp', 'NHBA', 'NadA', 'BAST']
	print(sep.join(headers))
	for f in args.fasta:
		# Defaults
		porACOUNT = '-'
		fetACOUNT = '-'
		fHbpCOUNT = '-'
		NHBACOUNT = '-'
		NadACOUNT = '-'
		bxtype = '-'
		# Standard run = serotype + ctrA
		seroCOUNT = '/'.join(seroTYPE(f, seroPRIMERS, allelesDB))
		ctrA_out = ctrA.ctrA_PCR(f, False, DBpath)
		ctrACOUNT = ctrA_out.split('\t')[1]
		# BAST
		if args.bast:
			ftRESULTS = fineTYPE(f, finetypePRIMERS, porA1DB, porA2DB, fetDB)
			bxRESULTS = bxTYPE(f, bxPRIMERS, fHbpDB, NHBADB, NadADB)
			porACOUNT = '/'.join(ftRESULTS[0])
			fetACOUNT = '/'.join(ftRESULTS[1])
			fHbpCOUNT = '/'.join(bxRESULTS[0])
			NHBACOUNT = '/'.join(bxRESULTS[1])
			NadACOUNT = '/'.join(bxRESULTS[2])
			bxallele = fHbpCOUNT + dl + NHBACOUNT + dl + NadACOUNT + dl + porACOUNT
			if bxallele in BAST:
				bxtype = BAST[bxallele]
		# Finetyping (porA, fetA, porB)
		elif args.finetype:
			ftRESULTS = fineTYPE(f, finetypePRIMERS, porA1DB, porA2DB, fetDB)
			porACOUNT = '/'.join(ftRESULTS[0])
			fetACOUNT = '/'.join(ftRESULTS[1])
		results = [f, seroCOUNT, ctrACOUNT, porACOUNT, fetACOUNT, fHbpCOUNT, NHBACOUNT, NadACOUNT, bxtype]
		print(sep.join(results))

	# Print allele sequences to file
	if args.printseq:
		if porASEQS:
			with open('printseq/porA_seqs.fasta', 'w') as output:
				SeqIO.write(porASEQS, output, 'fasta')
		if fetASEQS:
			with open('printseq/fetA_seqs.fasta', 'w') as output:
				SeqIO.write(fetASEQS, output, 'fasta')
		if fHbpSEQS:
			with open('printseq/fHbp_seqs.fasta', 'w') as output:
				SeqIO.write(fHbpSEQS, output, 'fasta')
		if NHBASEQS:
			with open('printseq/NHBA_seqs.fasta', 'w') as output:
				SeqIO.write(NHBASEQS, output, 'fasta')
		if NadASEQS:
			with open('printseq/NadA_seqs.fasta', 'w') as output:
				SeqIO.write(NadASEQS, output, 'fasta')
		msg('Done. Allele sequences saved to "printseq" folder.')

	sys.exit(0)

if __name__ == "__main__":
	main()
