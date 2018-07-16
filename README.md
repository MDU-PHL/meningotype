# meningotype

*In silico* typing of *Neisseria meningitidis* contigs
- Serotyping
- MLST  
- Finetyping (*porA*, *fetA*, *porB*)  
- Bexsero antigen sequence typing (BAST) (*fHbp*, *NHBA*, *NadA*, *PorA*)

## Quick start

```
# install
$ pip install git+https://github.com/MDU-PHL/meningotype.git

# just serotype
$ meningotype NMA.fasta

SAMPLE_ID  SEROGROUP  ctrA    MLST    PorA    FetA    PorB    fHbp    NHBA    NadA    BAST
NMA.fasta  A          ctrA    -       -       -       -       -       -       -       -

# include all genotypes
$ meningotype --all NMA.fasta

SAMPLE_ID  SEROGROUP  ctrA    MLST    PorA       FetA    PorB            fHbp    NHBA    NadA    BAST
NMA.fasta  A          ctrA    4       7,13-1     F1-5    NEIS2020_28     5       29      0       639

# type lots of files at once
$ meningotype --all *.fna

SAMPLE_ID  SEROGROUP  ctrA    MLST    PorA       FetA    PorB            fHbp    NHBA    NadA    BAST
A.fna      A          ctrA    4       7,13-1     F1-5    NEIS2020_28     5       29      0       639
B.fna      B          ctrA    8       5-2,10-1   F3-6    NEIS2020_12     16      20      8       150
C.fna      C          ctrA    177     21,26-2    F1-5    NEIS2020_3      17      101     9       118
W.fna      W          ctrA    11      5,2        F1-1    NEIS2020_244    623     29      6       141
X.fna      X          ctrA    181     5-1,10-1   F4-23   NEIS2020_509    391     358     0       -
Y.fna      Y          ctrA    23      5-2,10-1   F4-1    NEIS2020_67     25      7       0       228
```

## Installation

### Dependencies

* [Python 2.7.x](https://www.python.org/)
* [BioPython](http://biopython.org/)
* [isPcr v33](http://hgwdev.cse.ucsc.edu/~kent/src/)
* [NCBI BLAST+ >= 2.4](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [mlst](https://github.com/tseemann/mlst)

The simplest way to install dependencies is to use the Brew (MacOS) or
Linuxbrew (Linux) packaging system. 

```
$ brew tap brewsci/bio
$ brew install ispcr blast mlst
```

### Installing

The easiest way of installing `meningotype` is using `pip`:
```
$ pip install --user git+https://github.com/MDU-PHL/meningotype.git
```
 
The `--user` option will install the package locally, rather than in the global `python` directory. 

Thus, by default, this will install the package in `$HOME/.local/`, and the executable in `$HOME/.local/bin/`. 
To install the executable in a custom location (e.g., `$HOME/bin`), use the following:
```
$ pip install --install-option="--install-scripts=$HOME/bin" --user git+https://github.com/MDU-PHL/meningotype.git
```

To upgrade to a newer version: 
```
$ pip install --upgrade --install-option="--install-scripts=$HOME/bin" --user git+https://github.com/MDU-PHL/meningotype.git
```

### Testing

Once installed, you can run the following to ensure `meningotype` is successfully working:

    $ meningotype.py --test

If everything works, you will see the following:

```
$ meningotype.py --test
Running meningotype.py on test examples ... 
$ meningotype.py A.fna B.fna C.fna W.fna X.fna Y.fna
SAMPLE_ID	SEROGROUP	ctrA	MLST	PorA	FetA	PorB	fHbp	NHBA	NadA	BAST
meningotype/test/A.fna	A	ctrA	-	-	-	-	-	-	-	-
meningotype/test/B.fna	B	ctrA	-	-	-	-	-	-	-	-
meningotype/test/C.fna	C	ctrA	-	-	-	-	-	-	-	-
meningotype/test/W.fna	W	ctrA	-	-	-	-	-	-	-	-
meningotype/test/X.fna	X	ctrA	-	-	-	-	-	-	-	-
meningotype/test/Y.fna	Y	ctrA	-	-	-	-	-	-	-	-
```

or to check finetyping:

```
$ meningotype.py --test --finetype
Running meningotype.py on test examples ... 
$ meningotype.py A.fna B.fna C.fna W.fna X.fna Y.fna
SAMPLE_ID	SEROGROUP	ctrA	MLST	PorA	        FetA	PorB	fHbp	NHBA	NadA	BAST
meningotype/test/A.fna	A	ctrA	-	7,13-1		F1-5	-	-	-	-	-
meningotype/test/B.fna	B	ctrA	-	5-2,10-1	F3-6	-	-	-	-	-
meningotype/test/C.fna	C	ctrA	-	21,26-2		F1-5	-	-	-	-	-
meningotype/test/W.fna	W	ctrA	-	5,2		F1-1	-	-	-	-	-
meningotype/test/X.fna	X	ctrA	-	5-1,10-1	F4-23	-	-	-	-	-
meningotype/test/Y.fna	Y	ctrA	-	5-2,10-1	F4-1	-	-	-	-	-
```

or to check finetyping and Bexsero antigen sequence typing:

```
$ meningotype.py --test --all
Running meningotype.py on test examples ... 
$ meningotype.py A.fna B.fna C.fna W.fna X.fna Y.fna
SAMPLE_ID	SEROGROUP	ctrA	MLST	PorA		FetA	PorB		fHbp	NHBA	NadA	BAST
meningotype/test/A.fna	A	ctrA	4	7,13-1		F1-5	NEIS2020_28		5	29	0	639
meningotype/test/B.fna	B	ctrA	8	5-2,10-1	F3-6	NEIS2020_12		16	20	8	150
meningotype/test/C.fna	C	ctrA	177	21,26-2		F1-5	NEIS2020_3		17	101	9	118
meningotype/test/W.fna	W	ctrA	11	5,2		F1-1	NEIS2020_244	623	29	6	141
meningotype/test/X.fna	X	ctrA	181	5-1,10-1	F4-23	NEIS2020_509	391	358	0	-
meningotype/test/Y.fna	Y	ctrA	23	5-2,10-1	F4-1	NEIS2020_67		25	7	0	228
```

## Usage

```
$ meningotype.py -h
usage: 
  meningotype.py [OPTIONS] <fasta1> <fasta2> <fasta3> ... <fastaN>

In silico typing for Neisseria meningitidis
Default: Serotyping, MLST and ctrA PCR

PCR Serotyping Ref: Mothershed et al, J Clin Microbiol 2004; 42(1): 320-328
PorA and FetA typing Ref: Jolley et al, FEMS Microbiol Rev 2007; 31: 89-96
Bexsero antigen sequence typing (BAST) Ref: Brehony et al, Vaccine 2016; 34(39): 4690-4697
See also http://www.neisseria.org/nm/typing/

positional arguments:
  FASTA       input FASTA files eg. fasta1, fasta2, fasta3 ... fastaN

optional arguments:
  -h, --help  show this help message and exit
  --finetype  perform porA and fetA fine typing (default=off)
  --porB      perform porB sequence typing (NEIS2020) (default=off)
  --bast      perform Bexsero antigen sequence typing (BAST) (default=off)
  --mlst      perform MLST (default=off)
  --all       perform MLST, porA, fetA, porB, BAST typing (default=off)
  --db DB     specify custom directory containing allele databases for porA/fetA typing
              directory must contain database files: "FetA_VR.fas", "PorA_VR1.fas", "PorA_VR2.fas"
              for Bexsero typing: "fHbp_peptide.fas", "NHBA_peptide.fas", "NadA_peptide.fas", "BASTalleles.txt"
  --printseq  save porA/fetA or BAST allele sequences to file (default=off)
  --updatedb  update allele database from <pubmlst.org>
  --test      run test example
  --version   show program's version number and exit
```

## Examples

To perform *in silico* serotyping on FASTA files:

	$ meningotype <fasta1> <fasta2> <fasta3> ... <fastaN>`

The serotypes are printed in tab-separated format to `stdout`.
To save results to a tab-separated text file, redirect `stdout`:

	$ meningotype <fasta1> <fasta2> <fasta3> ... <fastaN>  > results.txt

To perform *in silico* serotyping AND finetyping of the porA and fetA genes:

	$ meningotype --finetype <fasta1> <fasta2> <fasta3> ... <fastaN>

To save finetyping sequences of the alleles to a file (eg. for uploading 
"new" sequences to [http://pubmlst.org/neisseria/](http://pubmlst.org/neisseria/)):

	$ meningotype --finetype --printseq <fasta1> <fasta2> <fasta3> ... <fastaN>

These are placed into a folder called `printseq` in the current directory.


## Updating the allele databases

To update the allele databases from http://pubmlst.org/neisseria/

	$ meningotype.py --updatedb

A copy of the original database is saved to `*.old` just in case, 
but is overwritten with each subsequent `--updatedb`.
*Ensure you back up your old databases if you wish to keep them.*

## Citation

Kwong JC, Gonçalves da Silva A, Stinear TP, Howden BP,  Seemann T.  
***meningotype*: *in silico* typing for *Neisseria meningitidis*.**  
GitHub https://github.com/MDU-PHL/meningotype

## Bugs

* Software - submit via the [GitHub issues page](https://github.com/MDU-PHL/meningotype/issues).  
* Database - contact the [pubmlst curator](mailto:keith.jolley@zoo.ox.ac.uk)

## Software Licence

[GPL3](https://github.com/MDU-PHL/meningotype/blob/master/LICENSE)

## Authors

* Jason Kwong (@kwongjc)
* Anders Gonçalves da Silva (@drandersg)
* Torsten Seemann (@torstenseemann)

## References

* [PubMLST Neisseria Database](https://pubmlst.org/neisseria/).  
* [Mothershed et al. J Clin Microbiol. 2004; 42(1):320-328.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC321732/)
* [Jolley et al. FEMS Microbiol Rev. 2007; 31:89-96.](http://onlinelibrary.wiley.com/doi/10.1111/j.1574-6976.2006.00057.x/full)  
* [Brehony et al. Vaccine. 2016; 34(39):4690-4697.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5012890/)
* [Bambini et al. PLoS One. 2013; 8(5):e65043.](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0065043)

