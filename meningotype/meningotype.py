#!/usr/bin/env python

'''
Script by Jason Kwong (refactored by Andreas Stroehlein in v0.9.0)
In silico typing for Neisseria meningitidis
'''

import argparse
from argparse import RawTextHelpFormatter
import sys
import os
import os.path
import re
import urllib
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline, NcbiblastxCommandline
from pkg_resources import resource_filename
from io import StringIO
from . import run_subprocess as rs
from . import menwy, ctrA, finetype, check_deps
from . import __version__ as version

# URLs to update database files (REST API)
porA1URL = 'http://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/PorA_VR1/alleles_fasta'
porA2URL = 'http://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/PorA_VR2/alleles_fasta'
fetAURL = 'http://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/FetA_VR/alleles_fasta'
porBURL = 'http://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/NEIS2020/alleles_fasta'
fHbpURL = 'http://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/fHbp_peptide/alleles_fasta'
NHBAURL = 'http://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/NHBA_peptide/alleles_fasta'
NadAURL = 'http://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/NadA_peptide/alleles_fasta'
BASTURL = 'http://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/53/profiles_csv'

# allele sizes and serotype dictionary
alleleSIZE = {'A': 92, 'B': 169, 'C': 74, 'W': 129, 'X': 65, 'Y': 146}
seroDICT = {'sacB': 'A', 'synD': 'B', 'synE': 'C', 'synG': 'W', 'xcbB': 'X', 'synF': 'Y'}

porASEQS = []
fetASEQS = []
porBSEQS = []
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
    sys.exit(1)

def update_db(db_file, db_url):
    '''
    Update database files
    '''
    if os.path.isfile(db_file):
        os.rename(db_file, db_file+'.old')
    urllib.request.urlretrieve(db_url, db_file)



def check_primer_files(f):
    '''
    Check primer files are present
    '''
    if not os.path.isfile(f):
        err('ERROR: Cannot locate file: "{}"'.format(f))

def check_db_files(f, db_url):
    if not os.path.isfile(f):
        update_db(f, db_url)


def check_fasta(f):
    if not os.path.isfile(f) or os.path.getsize(f) < 1:
        return False
    with open(f, 'r') as fasta:
        if fasta.readline()[0] != '>':
            return False
        for line in fasta:
            line = line.strip()
            if not line or line[0] == '>':
                continue
            if bool(re.search('[^ACTGactgNn-]', line)):
                return False
    return True

def makeblastDB(db, infile, dbtype):
    '''
    Set up BLAST databases if not present
    '''
    command = ['makeblastdb', '-in', infile, '-out', db, '-dbtype', dbtype]
    rs.run_subprocess(command)



def seroTYPE(f, seroprimers, allelesdb, cpus):
    '''
    Main serotyping functions
    '''
    seroCOUNT = [] # List for mixed/multiple hits
    command = ['isPcr', f, seroprimers, 'stdout', '-minPerfect=10']
    pcr_out = rs.run_subprocess(command)
    
    if not pcr_out:
        sero = None
        seroBLAST = NcbiblastnCommandline(query=f,
                                            db=allelesdb,
                                            task='blastn',
                                            perc_identity=90,
                                            evalue='1e-20',
                                            outfmt='"6 sseqid pident length"',
                                            culling_limit='1',
                                            num_threads=cpus)
        stdout, stderr = seroBLAST()
        lenMATCH = 0
        line = stdout.split('\n')[0]
        amp = line.split('\t')
        sero = amp[0]            # Currently only takes top/first BLAST hit
        if not sero:
            seroCOUNT.append('-')
        else:
            if sero == 'W' or sero == 'Y':
                sero = seroWY(f, sero)
            seroCOUNT.append(sero)
    else:
        alleleSEQ = StringIO()
        alleleSEQ.write(pcr_out)
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

def nm_mlst(fasta_input):
    '''
    Run mlst on fasta file
    '''
    command = ['mlst', '--scheme=neisseria', '--quiet', fasta_input]
    mlst_out = rs.run_subprocess(command)
    return mlst_out.split('\t')[2]

def finetypeBLAST(s, db, cpus):
    ft = None
    allele = None
    ftBLAST = NcbiblastxCommandline(query='-', db=db, outfmt='"6 qseqid sseqid pident length slen gaps nident evalue"', seg='no', query_gencode='11', matrix='PAM30', ungapped='true', comp_based_stats='0', evalue='1e-2', num_threads=cpus)        # blastx command to fix finding short sequences
    stdout, stderr = ftBLAST(stdin=s.format('fasta'))
    if stdout:
        BLASTout = stdout.split('\n')
        lenMATCH = 0
        for line in BLASTout:
            if line.strip():
                BLASTline = line.split('\t')
                if BLASTline[5] == '0':                                                        # check no gaps
                    if int(BLASTline[3]) == int(BLASTline[4]) == int(BLASTline[6]):            # check length of db subject = length of match = no. identical matches
                        if int(BLASTline[3]) > int(lenMATCH):                                # select longest match
                            lenMATCH = BLASTline[3]
                            ftRESULT = BLASTline[1]
                            ft = ftRESULT.split('_')[2]
        if not ft:                                                                            # if amplicon detected, but no match in db, assign as new
            ft = 'new'
            
    return str(ft)

def porBTYPE(f, blastdb, cpus):
    porB_result = finetype.porBBLAST(f, blastdb, cpus)
    porBCOUNT = porB_result[2]
    porBSEQS.append(porB_result[5])
    return porBCOUNT

def bxtypeBLAST(s, db, cpus):
    bx = None
    allele = None
    bxBLAST = NcbiblastxCommandline(query='-', db=db, outfmt='"6 qseqid sseqid pident length slen gaps nident evalue"', seg='no', culling_limit='1', evalue='1e-100', query_gencode='11', num_threads=cpus)
    stdout, stderr = bxBLAST(stdin=s.format('fasta'))
    if stdout:
        BLASTout = stdout.split('\n')
        lenMATCH = 0
        for line in BLASTout:
            if line.strip():
                BLASTline = line.split('\t')
                if BLASTline[5] == '0':                                                        # check no gaps
                    if int(BLASTline[3]) == int(BLASTline[4]) == int(BLASTline[6]):            # check length of db subject = length of match = no. identical matches
                        if int(BLASTline[3]) > int(lenMATCH):                                # select longest match
                            lenMATCH = BLASTline[3]
                            bxRESULT = BLASTline[1]
                            bx = bxRESULT.split('_')[2]
        if not bx:                                                                            # if amplicon detected, but no match in db, assign as new
            bx = 'new'
    else:                                                                                    # if no amplicon detected, assign as 0
        bx = '0'
    return str(bx)

def fineTYPE(f, finetypeprimers, poradb, pora1db, pora2db, fetdb, cpus):
    porACOUNT = []                            # Setup list in case there are mixed/multiple hits
    fetACOUNT = []
    global porASEQS
    global fetASEQS

    command = ['isPcr', f, finetypeprimers, 'stdout', '-maxSize=800',
               '-tileSize=10', '-minPerfect=8', '-stepSize=3']
    pcr_out = rs.run_subprocess(command)

    alleleSEQ = StringIO()
    alleleSEQ.write(pcr_out)
    alleleSEQ.seek(0)

    for amplicon in SeqIO.parse(alleleSEQ, "fasta"):
        # msg(amplicon)
        ampFILE = amplicon.description.split()
        ampID = ampFILE[1]
        ampLEN = ampFILE[2]        
        # Need to check amplicon length to exclude double hits?
        if ampID == 'porA':
            A1 = finetypeBLAST(amplicon, pora1db, cpus)
            A2 = finetypeBLAST(amplicon, pora2db, cpus)
            porA = A1 + ',' + A2
            porACOUNT.append(porA)
            porASEQ = amplicon.seq.upper()
            porARECR = SeqRecord(porASEQ, id=f, description='PorA')
            porASEQS.append(porARECR)
        if ampID == 'fetA':
            fet = finetypeBLAST(amplicon, fetdb, cpus)
            fetACOUNT.append(fet)
            fetASEQ = amplicon.seq.upper()
            fetARECR = SeqRecord(fetASEQ, id=f, description='FetA')
            fetASEQS.append(fetARECR)
    
    if len(porACOUNT) == 0:
        porseqBLAST = NcbiblastnCommandline(query=f, db=poradb, perc_identity=90, outfmt='"6 qseq"', culling_limit='1', num_threads=cpus)
        stdout, stderr = porseqBLAST()
        if stdout:
            porAseq = Seq(stdout.strip())
            porArec = SeqRecord(porAseq, id=f, description='PorA')
            A1 = finetypeBLAST(porArec, pora1db, cpus)
            A2 = finetypeBLAST(porArec, pora2db, cpus)
            porA = A1 + ',' + A2
            porACOUNT.append(porA)
            porASEQ = amplicon.seq.upper()
            porARECR = SeqRecord(porASEQ, id=f, description='PorA')
            porASEQS.append(porARECR)
        if len(porACOUNT) == 0:
            porACOUNT.append('-')
    if len(fetACOUNT) == 0:
        fetACOUNT.append('-')
    return porACOUNT, fetACOUNT


def bxTYPE(f, bxPRIMERS, fHbpDB, NHBADB, NadADB, cpus):
    fHbpCOUNT = [] # List for mixed/multiple hits
    NHBACOUNT = []
    NadACOUNT = []
    global fHbpSEQS
    global NHBASEQS
    global NadASEQS
    
    command = ['isPcr', f, bxPRIMERS, 'stdout', '-maxSize=3000', '-tileSize=7', '-minPerfect=8', '-stepSize=3']
    pcr_out = rs.run_subprocess(command)
        
    alleleSEQ = StringIO()
    alleleSEQ.write(pcr_out)
    alleleSEQ.seek(0)
    for amplicon in SeqIO.parse(alleleSEQ, "fasta"):
        ampFILE = amplicon.description.split()
        ampID = ampFILE[1]
        ampLEN = ampFILE[2] # Need to check amplicon length to exclude double hits?
        if ampID == 'fHbp':
            fHbp = bxtypeBLAST(amplicon, fHbpDB, cpus)
            fHbpCOUNT.append(fHbp)
            fHbpSEQ = amplicon.seq.upper()
            fHbpRECR = SeqRecord(fHbpSEQ, id=f, description='fHbp')
            fHbpSEQS.append(fHbpRECR)
        if ampID == 'NHBA':
            NHBA = bxtypeBLAST(amplicon, NHBADB, cpus)
            NHBACOUNT.append(NHBA)
            NHBASEQ = amplicon.seq.upper()
            NHBARECR = SeqRecord(NHBASEQ, id=f, description='NHBA')
            NHBASEQS.append(NHBARECR)
        if ampID == 'NadA':
            NadA = bxtypeBLAST(amplicon, NadADB, cpus)
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

def main():
    '''
    Main function for command line call
    '''
    parser = argparse.ArgumentParser(
        formatter_class=RawTextHelpFormatter,
        description='In silico typing for Neisseria meningitidis\n'
            'Default: Serotyping and ctrA PCR\n'
            '\nPCR Serotyping Ref: Mothershed et al, J Clin Microbiol 2004; 42(1): 320-328\n'
            'PorA and FetA typing Ref: Jolley et al, FEMS Microbiol Rev 2007; 31: 89-96\n'
            'Bexsero antigen sequence typing (BAST) Ref: Brehony et al, Vaccine 2016; 34(39): 4690-4697\n'
            'MenDeVAR Vaccine Reactivity Index Ref: Rodrigues et al. J Clin Microbiol. 2020; 59(1):e02161-20\n'
            'See also https://pubmlst.org/organisms/neisseria-spp',
        usage='\n  %(prog)s [OPTIONS] <fasta1> <fasta2> <fasta3> ... <fastaN>')
    parser.add_argument('fasta', metavar='FASTA', nargs='*', help='input FASTA files eg. fasta1, fasta2, fasta3 ... fastaN')
    parser.add_argument('--finetype', action='store_true', help='perform porA and fetA fine typing (default=off)')
    parser.add_argument('--porB', action='store_true', help='perform porB sequence typing (NEIS2020) (default=off)')
    parser.add_argument('--bast', action='store_true', help='perform Bexsero antigen sequence typing (BAST) and determine MenDeVAR index (default=off)')
    parser.add_argument('--mlst', action='store_true', help='perform MLST (default=off)')
    parser.add_argument('--all', action='store_true', help='perform MLST, porA, fetA, porB, BAST typing and determine MenDeVAR index (default=off)')
    parser.add_argument('--db', metavar='DB', help='specify custom directory containing allele databases for porA/fetA typing\n'
        'directory must contain database files: "FetA_VR.fas", "PorA_VR1.fas", "PorA_VR2.fas"\n'
        'For Bexsero typing and MenDeVAR: "fHbp_peptide.fas", "NHBA_peptide.fas", "NadA_peptide.fas", "BASTalleles.txt"\n')
    parser.add_argument('--printseq', metavar='DIR', help='specify directory to save extracted porA/fetA/porB or BAST allele sequences (default=off)')
    parser.add_argument('--cpus', metavar='CPUS', default=1, help='number of cpus to use in BLAST search (default=1)')
    parser.add_argument('--updatedb', action='store_true', default=False, help='update allele database from <pubmlst.org>')
    parser.add_argument('--test', action='store_true', default=False, help='run test example')
    parser.add_argument('--checkdeps', action='store_true', default=False, help='check dependencies are installed and exit')
    parser.add_argument('--version', action='version', version=
        '%(prog)s v{}\n'.format(version))
    args = parser.parse_args()
    cpus = int(args.cpus)

    # Check dependencies
    dependencies = ['isPcr', 'blastn', 'blastx', 'mlst']
    if args.checkdeps:
        msg("\033[91mChecking dependencies:\033[0m")
        for dep in dependencies:
            if check_deps.which(dep):
                msg(' ....... '.join([dep, '\tFound "{}"'.format(check_deps.which(dep)), '\t[OK]']))
            else:
                msg(' ........ '.join([dep, 'Not found - please check {} is installed and in the path.'.format(dep), '[ERROR]']))
        sys.exit(0)
    dep_errors = 0
    for dep in dependencies:
        if not check_deps.which(dep):
            msg('ERROR: {} not found - please check it is installed and in the path.'.format(dep))
            dep_errors = dep_errors + 1
    if dep_errors != 0:
        sys.exit(1)

    # Set path for database files
    if args.db:
        DBpath = str(args.db).rstrip('/')
    else:
        DBpath = resource_filename(__name__, 'db')

    # Path to database files
    porA1alleles = os.path.join( DBpath, 'PorA_VR1.fas' )
    porA2alleles = os.path.join( DBpath, 'PorA_VR2.fas' )
    fetalleles = os.path.join( DBpath, 'FetA_VR.fas' )
    porBalleles = os.path.join( DBpath, 'PorB.fas' )
    fHbpalleles = os.path.join( DBpath, 'fHbp_peptide.fas' )
    NHBAalleles = os.path.join( DBpath, 'NHBA_peptide.fas' )
    NadAalleles = os.path.join( DBpath, 'NadA_peptide.fas' )
    BASTalleles = os.path.join( DBpath, 'BASTalleles.txt' )

    seroALLELES = os.path.join( DBpath, 'seroALLELES.fa' )
    seroPRIMERS = os.path.join( DBpath, 'seroPRIMERS' )
    finetypePRIMERS = os.path.join( DBpath, 'finetypePRIMERS' )
    bxPRIMERS = os.path.join( DBpath, 'bexseroPRIMERS' )

    allelesDB = os.path.join( DBpath, 'blast', 'seroALLELES' )
    porADB = os.path.join( DBpath, 'blast', 'porA' )
    porA1DB = os.path.join( DBpath, 'blast', 'porA1' )
    porA2DB = os.path.join( DBpath, 'blast', 'porA2' )
    fetDB = os.path.join( DBpath, 'blast', 'fet' )
    fHbpDB = os.path.join( DBpath, 'blast', 'fHbp_peptide' )
    NHBADB = os.path.join( DBpath, 'blast', 'NHBA_peptide' )
    NadADB = os.path.join( DBpath, 'blast', 'NadA_peptide' )
    porBDB = os.path.join( DBpath, 'blast', 'porB' )

    if args.updatedb:
        try:
            msg('Updating "{}" ... '.format(porA1alleles))
            update_db(porA1alleles, porA1URL)
            makeblastDB(porA1DB, porA1alleles, 'prot')
            msg('Updating "{}" ... '.format(porA2alleles))
            update_db(porA2alleles, porA2URL)
            makeblastDB(porA2DB, porA2alleles, 'prot')
            msg('Updating "{}" ... '.format(fetalleles))
            update_db(fetalleles, fetAURL)
            makeblastDB(fetDB, fetalleles, 'prot')
            msg('Updating "{}" ... '.format(porBalleles))
            update_db(porBalleles, porBURL)
            makeblastDB(porBDB, porBalleles, 'nucl')
            msg('Updating "{}" ... '.format(fHbpalleles))
            update_db(fHbpalleles, fHbpURL)
            makeblastDB(fHbpDB, fHbpalleles, 'prot')
            msg('Updating "{}" ... '.format(NHBAalleles))
            update_db(NHBAalleles, NHBAURL)
            makeblastDB(NHBADB, NHBAalleles, 'prot')
            msg('Updating "{}" ... '.format(NadAalleles))
            update_db(NadAalleles, NadAURL)
            makeblastDB(NadADB, NadAalleles, 'prot')
            msg('Updating "{}" ... '.format(BASTalleles))
            update_db(BASTalleles, BASTURL)
            msg('Done.')
        except IOError:
            err('ERROR: Unable to update DB at "{}".\nCheck DB directory permissions and connection to http://pubmlst.org.'.format(DBpath))
        except urllib.error.HTTPError:
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
    check_db_files(porBalleles, porBURL)

    # Setup BLASTDB
    if not os.path.isfile(allelesDB + '.nin'):
        makeblastDB(allelesDB, seroALLELES, 'nucl')

    if args.bast or args.all:
        check_primer_files(bxPRIMERS)
        check_db_files(fHbpalleles, fHbpURL)
        check_db_files(NHBAalleles, NHBAURL)
        check_db_files(NadAalleles, NadAURL)
        check_db_files(BASTalleles, BASTURL)
        if not os.path.isfile(fHbpDB + '.pin'):
            makeblastDB(fHbpDB, fHbpalleles, 'prot')
        if not os.path.isfile(NHBADB + '.pin'):
            makeblastDB(NHBADB, NHBAalleles, 'prot')
        if not os.path.isfile(NadADB + '.pin'):
            makeblastDB(NadADB, NadAalleles, 'prot')

        # Import allele profiles as dictionary
        bexs_dict = {}
        trum_dict = {}
        BAST = {}

        with open(BASTalleles) as db:
            for lines in db:
                line = lines.rstrip('\n').split('\t')
                ST = line[0]
                alleles = tuple(line[1:6])
                BAST[alleles] = ST
                bexs_dict[ST] = line[6]
                trum_dict[ST] = line[7]
                

    # Test example to check meningotype works
    if args.test:
        TESTpath = resource_filename(__name__, 'test')
        testSEQS = [os.path.join( TESTpath, f ) for f in ['A.fna', 'B.fna', 'C.fna', 'W.fna', 'X.fna', 'Y.fna'] ]
        msg('\033[94mRunning meningotype.py on test examples ... \033[0m')
        msg('$ meningotype.py A.fna B.fna C.fna W.fna X.fna Y.fna')
        args.fasta = testSEQS

    # Create folder for output sequences if specified
    # Test if directory exists and is writeable before running meningotype
    if args.printseq:
        try:
            if not os.path.exists(args.printseq):
                os.makedirs(args.printseq)
            else:
                err('ERROR: "{}" already exists.'.format(args.printseq))
        except OSError:
            err('ERROR: Unable to create "{}" in this directory.'.format(args.printseq))

    # Run meningotype
    if len(args.fasta) == 0:
        message = "\033[91mEither use --test or specify at least one FASTA file.\033[0m"
        sys.stderr.write('error: {}\n'.format( message ) )
        parser.print_help()
        parser.exit(1)

    # DEFINING OUTPUT COLUMNS

    # 1       SAMPLE_ID
    # 2       SEROGROUP
    # 3       ctrA
    # 4       MLST
    # 5       PorA
    # 6       FetA
    # 7       PorB
    # 8       fHbp
    # 9       NHBA
    # 10      NadA
    # 11      BAST
    # 12      Bexsero MenDeVAR
    # 13      Trumenba MenDeVAR

    # Default 1-3
    out_cols = {0,1,2}

    # args.finetype add 5,6 columns
    if args.finetype:
        out_cols.update({4,5})

    # args.bast add 5,6, 8-13
    if args.bast:
        out_cols.update({4,5})
        for pos in range(7,13):
            out_cols.update({pos})

    # args.mlst Add 4th columns
    if args.mlst:
        out_cols.update({3})

    # args.all add all columns 4-13
    if args.all:
        for pos in range(3,13):
            out_cols.update({pos})

    # args.porB add 7th column
    if args.porB:
        out_cols.update({6})

    headers = ['SAMPLE_ID', 'SEROGROUP', 'ctrA', 'MLST', 'PorA', 'FetA', 'PorB', 'fHbp', 'NHBA', 'NadA', 'BAST', 'Bexsero MenDeVAR', 'Trumenba MenDeVAR']
    # only show headers of the output that is created
    subs_headers = [headers[i] for i in out_cols]

    print(sep.join(subs_headers))


    for f in args.fasta:
        if check_fasta(f) != True:
            print('{}\tERROR: Check file exists and is in FASTA nucleotide format.'.format(f))
            continue
        # Defaults
        mlst = '-'
        porACOUNT = '-'
        fetACOUNT = '-'
        porBCOUNT = '-'
        fHbpCOUNT = '-'
        NHBACOUNT = '-'
        NadACOUNT = '-'
        bxtype = '-'
        bexsero = '-'
        trumenba = '-'
        # Standard run = serotype + ctrA
        seroCOUNT = '/'.join(seroTYPE(f, seroPRIMERS, allelesDB, cpus))
        ctrA_out = ctrA.ctrA_PCR(f, False, DBpath)
        ctrACOUNT = ctrA_out.split('\t')[1]
        # Optional typing
        if args.porB or args.all:
            porBCOUNT = porBTYPE(f, porBDB, cpus)
        if args.mlst or args.all:
            mlst = nm_mlst(f)
        # BAST
        if args.bast or args.all:
            ftRESULTS = fineTYPE(f, finetypePRIMERS, porADB, porA1DB, porA2DB, fetDB, cpus)
            bxRESULTS = bxTYPE(f, bxPRIMERS, fHbpDB, NHBADB, NadADB, cpus)
            porACOUNT = '/'.join(ftRESULTS[0])
            fetACOUNT = '/'.join(ftRESULTS[1])
            fHbpCOUNT = '/'.join(bxRESULTS[0])
            NHBACOUNT = '/'.join(bxRESULTS[1])
            NadACOUNT = '/'.join(bxRESULTS[2])
            bxallele = tuple([fHbpCOUNT] + [NHBACOUNT] + [NadACOUNT] + porACOUNT.split(","))

            bxtype = BAST.get(bxallele, "-")
            bexsero = bexs_dict.get(bxtype, "insufficient data")
            trumenba = trum_dict.get(bxtype, "insufficient data")


        # Finetyping (porA, fetA, porB)
        elif args.finetype:
            ftRESULTS = fineTYPE(f, finetypePRIMERS, porADB, porA1DB, porA2DB, fetDB, cpus)
            porACOUNT = '/'.join(ftRESULTS[0])
            fetACOUNT = '/'.join(ftRESULTS[1])

        # Print results to stdout

        this_mlst = mlst[2] if isinstance(mlst, list) else mlst
        
        results = [f, seroCOUNT, ctrACOUNT, this_mlst, porACOUNT, fetACOUNT, porBCOUNT, fHbpCOUNT, NHBACOUNT, NadACOUNT, bxtype, bexsero, trumenba]
        # print only the result columns that have been requested via set args
        subs_results = [results[i] for i in out_cols]
        print(sep.join(subs_results))

    # Print allele sequences to file
    if args.printseq:
        if porASEQS:
            with open(os.path.join(args.printseq, 'porA_seqs.fasta'), 'w') as output:
                SeqIO.write(porASEQS, output, 'fasta')
        if fetASEQS:
            with open(os.path.join(args.printseq, 'fetA_seqs.fasta'), 'w') as output:
                SeqIO.write(fetASEQS, output, 'fasta')
        if porBSEQS:
            with open(os.path.join(args.printseq, 'porB_seqs.fasta'), 'w') as output:
                SeqIO.write(porBSEQS, output, 'fasta')
        if fHbpSEQS:
            with open(os.path.join(args.printseq, 'fHbp_seqs.fasta'), 'w') as output:
                SeqIO.write(fHbpSEQS, output, 'fasta')
        if NHBASEQS:
            with open(os.path.join(args.printseq, 'NHBA_seqs.fasta'), 'w') as output:
                SeqIO.write(NHBASEQS, output, 'fasta')
        if NadASEQS:
            with open(os.path.join(args.printseq, 'NadA_seqs.fasta'), 'w') as output:
                SeqIO.write(NadASEQS, output, 'fasta')
        msg('Done. Allele sequences saved in "{}/".'.format(args.printseq))

    sys.exit(0)

if __name__ == "__main__":
    main()
