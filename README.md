# meningotype

*In silico* typing of *Neisseria meningitidis* contigs
- Serotyping
- MLST  
- Finetyping (*porA*, *fetA*, *porB*)  
- Bexsero antigen sequence typing (BAST) (*fHbp*, *NHBA*, *NadA*, *PorA*)
- MenDeVAR (Meningococcal Deduced Vaccine Antigen Reactivity) Index

## Quick start

```
# install
$ pip install git+https://github.com/MDU-PHL/meningotype.git

# just serotype
$ meningotype *.fna

SAMPLE_ID  SEROGROUP  ctrA
NMA.fasta  A          ctrA

# type lots of files at once and include all typing results
$ meningotype --all *.fna

SAMPLE_ID       SEROGROUP       ctrA    MLST    PorA      FetA    PorB          fHbp    NHBA    NadA    BAST    Bexsero MenDeVAR        Trumenba MenDeVAR
A.fna           A               ctrA    4       7,13-1    F1-5    NEIS2020_28   5       29      0       639     insufficient data       insufficient data
B.fna           B               ctrA    8       5-2,10-1  F3-6    NEIS2020_12   16      20      8       150     exact match             cross-reactive
C.fna           C               ctrA    177     21,26-2   F1-5    NEIS2020_3    17      101     9       118     insufficient data       insufficient data
W.fna           W               ctrA    11      5,2       F1-1    NEIS2020_244  623     29      6       141     cross-reactive          insufficient data
X.fna           X               ctrA    181     5-1,10-1  F4-23   NEIS2020_509  391     358     0       3042    insufficient data       insufficient data
Y.fna           Y               ctrA    23      5-2,10-1  F4-1    NEIS2020_67   25      7       0       228     insufficient data       cross-reactive
```

## Installation

### Dependencies

* [Python >=3.7.x](https://www.python.org/)
* [BioPython](http://biopython.org/)
* [isPcr v33](http://hgwdev.cse.ucsc.edu/~kent/src/)
* [NCBI BLAST+ >= 2.4](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [mlst](https://github.com/tseemann/mlst)


### Installing

The easiest way of installing `meningotype` is using `conda`, which will install all required dependencies:
```
$ conda create -n meningotype -c bioconda meningotype
```

### Testing

Once installed, you can run the following to ensure `meningotype` is successfully working:

    $ meningotype --version
    meningotype v0.9.0

    $ meningotype --checkdeps
    Checking dependencies:
    isPcr .......   Found "/home/user/bin/isPcr" .......   [OK]
    blastn .......  Found "/home/user/bin/blastn" .......  [OK]
    blastx .......  Found "/home/user/bin/blastx" .......  [OK]
    mlst .......    Found "/home/user/bin/mlst" .......    [OK]


    $ meningotype --test

If everything works, you will see the following:

```
$ meningotype --test
Running meningotype on test examples ... 
$ meningotype A.fna B.fna C.fna W.fna X.fna Y.fna
SAMPLE_ID       SEROGROUP       ctrA
A.fna           A               ctrA
B.fna           B               ctrA
C.fna           C               ctrA
W.fna           W               ctrA
X.fna           X               ctrA
Y.fna           Y               ctrA
```

or to check finetyping:

```
$ meningotype --test --finetype
Running meningotype on test examples ... 
$ meningotype A.fna B.fna C.fna W.fna X.fna Y.fna
SAMPLE_ID       SEROGROUP       ctrA    PorA      FetA
A.fna           A               ctrA    7,13-1    F1-5
B.fna           B               ctrA    5-2,10-1  F3-6
C.fna           C               ctrA    21,26-2   F1-5
W.fna           W               ctrA    5,2       F1-1
X.fna           X               ctrA    5-1,10-1  F4-23
Y.fna           Y               ctrA    5-2,10-1  F4-1
```

or to check finetyping and Bexsero antigen sequence typing and MenDeVAR Vaccine Antigen Reactivity Index:

```
$ meningotype --test --all
Running meningotype.py on test examples ... 
$ meningotype.py A.fna B.fna C.fna W.fna X.fna Y.fna
SAMPLE_ID       SEROGROUP       ctrA    MLST    PorA      FetA    PorB          fHbp    NHBA    NadA    BAST    Bexsero MenDeVAR        Trumenba MenDeVAR
A.fna           A               ctrA    4       7,13-1    F1-5    NEIS2020_28   5       29      0       639     insufficient data       insufficient data
B.fna           B               ctrA    8       5-2,10-1  F3-6    NEIS2020_12   16      20      8       150     exact match             cross-reactive
C.fna           C               ctrA    177     21,26-2   F1-5    NEIS2020_3    17      101     9       118     insufficient data       insufficient data
W.fna           W               ctrA    11      5,2       F1-1    NEIS2020_244  623     29      6       141     cross-reactive          insufficient data
X.fna           X               ctrA    181     5-1,10-1  F4-23   NEIS2020_509  391     358     0       3042    insufficient data       insufficient data
Y.fna           Y               ctrA    23      5-2,10-1  F4-1    NEIS2020_67   25      7       0       228     insufficient data       cross-reactive
```

## Usage

```
$ meningotype -h
usage: 
  meningotype [OPTIONS] <fasta1> <fasta2> <fasta3> ... <fastaN>

In silico typing for Neisseria meningitidis
Default: Serotyping and ctrA PCR

PCR Serotyping Ref: Mothershed et al, J Clin Microbiol 2004; 42(1): 320-328
PorA and FetA typing Ref: Jolley et al, FEMS Microbiol Rev 2007; 31: 89-96
Bexsero antigen sequence typing (BAST) Ref: Brehony et al, Vaccine 2016; 34(39): 4690-4697
MenDeVAR Vaccine Reactivity Index Ref: Rodrigues et al. J Clin Microbiol. 2020; 59(1):e02161-20
See also https://pubmlst.org/organisms/neisseria-spp

positional arguments:
  FASTA       input FASTA files eg. fasta1, fasta2, fasta3 ... fastaN

optional arguments:
  -h, --help      show this help message and exit
  --finetype      perform porA and fetA fine typing (default=off)
  --porB          perform porB sequence typing (NEIS2020) (default=off)
  --bast          perform Bexsero antigen sequence typing (BAST) (default=off)
  --mlst          perform MLST (default=off)
  --all           perform MLST, porA, fetA, porB, BAST typing and MenDeVAR index (default=off)
  --db DB         specify custom directory containing allele databases for porA/fetA typing
                  directory must contain database files: "FetA_VR.fas", "PorA_VR1.fas", "PorA_VR2.fas"
                  for Bexsero typing: "fHbp_peptide.fas", "NHBA_peptide.fas", "NadA_peptide.fas", "BASTalleles.txt"
  --printseq DIR  specify directory to save extracted porA/fetA/porB or BAST allele sequences (default=off)
  --cpus CPUS     number of cpus to use in BLAST search (default=1)
  --updatedb      update allele database from <pubmlst.org>
  --test          run test example
  --checkdeps     check dependencies are installed and exit
  --version       show program's version number and exit
```

## Examples

To perform *in silico* serotyping on FASTA files:

	$ meningotype <fasta1> <fasta2> <fasta3> ... <fastaN>`

The serotypes are printed in tab-separated format to `stdout`.
To save results to a tab-separated text file, redirect `stdout`:

	$ meningotype <fasta1> <fasta2> <fasta3> ... <fastaN>  > results.tsv

To perform *in silico* serotyping AND finetyping of the porA and fetA genes:

	$ meningotype --finetype <fasta1> <fasta2> <fasta3> ... <fastaN>

To save finetyping sequences of the alleles to a file (eg. for uploading 
"new" sequences to [http://pubmlst.org/neisseria/](https://pubmlst.org/organisms/neisseria-spp/)):

	$ meningotype --finetype --printseq <fasta1> <fasta2> <fasta3> ... <fastaN>

These are placed into a folder called `printseq` in the current directory.

To perform *in silico* serotyping AND finetyping AND Bexsero antigen sequence typing AND determine the MenDeVAR Vaccine Antigen Reactivity Index

  `$ meningotype --all --printseq <fasta1> <fasta2> <fasta3> ... <fastaN>`


The MenDeVAR index assigns each record one of four categories:

1. Green / Exact match: isolate contains one or more exact sequence match(es) to antigenic variants found in the vaccine
2. Amber / Cross-reactive: isolate contains one or more antigenic variant deemed cross-reactive to vaccine variants through experimental studies.
3. Red / No match: all the isolate's antigenic variants have been deemed not cross-reactive to vaccine variants through experimental studies.
4. Grey / Insufficient data: isolate contains antigens for which there is insufficient data available or which are yet to be tested in experimental studies.


## Updating the allele databases

To update the allele databases from https://pubmlst.org/organisms/neisseria-spp/

	$ meningotype --updatedb

A copy of the original database is saved to `*.old`, 
but is overwritten with each subsequent `--updatedb`.
*Ensure you back up your old databases if you wish to keep them.*

## Citation

Kwong JC, Gonçalves da Silva A, Stinear TP, Howden BP,  Seemann T.  
***meningotype*: *in silico* typing for *Neisseria meningitidis*.**  
GitHub https://github.com/MDU-PHL/meningotype

## Bugs

* Software - submit via the [GitHub issues page](https://github.com/MDU-PHL/meningotype/issues).  
* Database - contact the [PubMLST curator](mailto:keith.jolley@zoo.ox.ac.uk)

## Software Licence

[GPL3](https://github.com/MDU-PHL/meningotype/blob/master/LICENSE)

## Authors

* Jason Kwong (@kwongjc)
* Andreas Stroehlein (@stroehleina)
* Anders Gonçalves da Silva (@drandersg)
* Torsten Seemann (@torstenseemann)

## References

* [PubMLST Neisseria Database](https://pubmlst.org/neisseria/).  
* [Mothershed et al. J Clin Microbiol. 2004; 42(1):320-328.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC321732/)
* [Jolley et al. FEMS Microbiol Rev. 2007; 31:89-96.](http://onlinelibrary.wiley.com/doi/10.1111/j.1574-6976.2006.00057.x/full)  
* [Brehony et al. Vaccine. 2016; 34(39):4690-4697.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5012890/)
* [Bambini et al. PLoS One. 2013; 8(5):e65043.](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0065043)
* [Rodrigues et al. J Clin Microbiol. 2020; 59(1):e02161-20.](https://journals.asm.org/doi/10.1128/JCM.02161-20)

