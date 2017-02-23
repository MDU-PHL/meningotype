#meningotype

*In silico* typing of *Neisseria meningitidis*  
- Serotyping
- Finetyping (*porA* and *fetA*)  
- Bexsero antigen sequence typing (BAST) (*fHbp*, *NHBA*, *NadA*, *PorA*)

##Authors

* Jason Kwong (@kwongjc)
* Anders Gonçalves da Silva
* Torsten Seemann (@torstenseemann)

##Dependencies

* [Python 2.7.x](https://www.python.org/)
* [BioPython](http://biopython.org/)
* [isPcr v33](http://hgwdev.cse.ucsc.edu/~kent/src/) by Jim Kent
* [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [mlst](https://github.com/tseemann/mlst)

##Installation

The easiest way of installing `meningotype` is using `pip`:

    $ pip install --user git+https://github.com/MDU-PHL/meningotype.git
    
The `--user` option will install the package locally, rather than in the global `python` directory. 

Thus, by default, this will install the package in `$HOME/.local/`, and the executable in `$HOME/.local/bin/`. To install the executable in a custom location (e.g., `$HOME/bin`), use the following:

    $ pip install --install-option="--install-scripts=$HOME/bin" --user git+https://github.com/MDU-PHL/meningotype.git

To upgrade to a newer version: 

    $ pip install --upgrade --install-option="--install-scripts=$HOME/bin" --user git+https://github.com/MDU-PHL/meningotype.git

The simplest way to install dependencies is to use the Brew (Mac OS X) or LinuxBrew (Linux) system. Users who have difficulty installing isPcr from source (eg. Mac OS) may have more success with Brew:
```
$ brew tap homebrew/science
$ brew tap chapmanb/cbl
$ brew tap tseemann/homebrew-bioinformatics-linux
```


### To test installation

Once installed, you can run the following to ensure `meningotype` is successfully working:

    $ meningotype.py --test

If everything works, you will see the following:

```
$ meningotype.py --test
Running meningotype.py on test examples ... 
$ meningotype.py A.fna B.fna C.fna W.fna X.fna Y.fna
SAMPLE_ID	SEROGROUP	ctrA	MLST	PorA	FetA	PorB	fHbp	NHBA	NadA	BAST
/home/jasonk1/scripts/git/meningotype/meningotype/test/A.fna	A	ctrA	4	-	-	-	-	-	-	-
/home/jasonk1/scripts/git/meningotype/meningotype/test/B.fna	B	ctrA	8	-	-	-	-	-	-	-
/home/jasonk1/scripts/git/meningotype/meningotype/test/C.fna	C	ctrA	177	-	-	-	-	-	-	-
/home/jasonk1/scripts/git/meningotype/meningotype/test/W.fna	W	ctrA	11	-	-	-	-	-	-	-
/home/jasonk1/scripts/git/meningotype/meningotype/test/X.fna	X	ctrA	181	-	-	-	-	-	-	-
/home/jasonk1/scripts/git/meningotype/meningotype/test/Y.fna	Y	ctrA	23	-	-	-	-	-	-	-
```

or to check finetyping:

```
$ meningotype.py --test --finetype
Running meningotype.py on test examples ... 
$ meningotype.py A.fna B.fna C.fna W.fna X.fna Y.fna
SAMPLE_ID	SEROGROUP	ctrA	MLST	PorA	FetA	PorB	fHbp	NHBA	NadA	BAST
/home/jasonk1/scripts/git/meningotype/meningotype/test/A.fna	A	ctrA	4	7,13-1		F1-5	NEIS2020_28		-	-	-	-
/home/jasonk1/scripts/git/meningotype/meningotype/test/B.fna	B	ctrA	8	5-2,10-1	F3-6	NEIS2020_12		-	-	-	-
/home/jasonk1/scripts/git/meningotype/meningotype/test/C.fna	C	ctrA	177	21,26-2		F1-5	NEIS2020_3		-	-	-	-
/home/jasonk1/scripts/git/meningotype/meningotype/test/W.fna	W	ctrA	11	5,2			F1-1	NEIS2020_244	-	-	-	-
/home/jasonk1/scripts/git/meningotype/meningotype/test/X.fna	X	ctrA	181	5-1,10-1	F4-23	NEIS2020_509	-	-	-	-
/home/jasonk1/scripts/git/meningotype/meningotype/test/Y.fna	Y	ctrA	23	5-2,10-1	F4-1	NEIS2020_67		-	-	-	-
```

or to check finetyping and Bexsero antigen sequence typing:

```
$ meningotype.py --test --finetype --bast
Running meningotype.py on test examples ... 
$ meningotype.py A.fna B.fna C.fna W.fna X.fna Y.fna
SAMPLE_ID	SEROGROUP	ctrA	MLST	PorA	FetA	PorB	fHbp	NHBA	NadA	BAST
/home/jasonk1/scripts/git/meningotype/meningotype/test/A.fna	A	ctrA	4	7,13-1		F1-5	NEIS2020_28		5	29	0	639
/home/jasonk1/scripts/git/meningotype/meningotype/test/B.fna	B	ctrA	8	5-2,10-1	F3-6	NEIS2020_12		16	20	8	150
/home/jasonk1/scripts/git/meningotype/meningotype/test/C.fna	C	ctrA	177	21,26-2		F1-5	NEIS2020_3		17	101	9	118
/home/jasonk1/scripts/git/meningotype/meningotype/test/W.fna	W	ctrA	11	5,2			F1-1	NEIS2020_244	623	29	6	141
/home/jasonk1/scripts/git/meningotype/meningotype/test/X.fna	X	ctrA	181	5-1,10-1	F4-23	NEIS2020_509	391	358	0	-
/home/jasonk1/scripts/git/meningotype/meningotype/test/Y.fna	Y	ctrA	23	5-2,10-1	F4-1	NEIS2020_67		25	7	0	228
```


##Usage

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
  FASTA               input FASTA files eg. fasta1, fasta2, fasta3 ... fastaN

optional arguments:
  -h, --help          show this help message and exit
  --finetype          perform porA and fetA fine typing (default=off)
  --porB              perform porB sequence typing (NEIS2020) (default=off)
  --bast              perform Bexsero antigen sequence typing (BAST) (default=off)
  --mlst on|off|only  toggles whether MLST is run or not, or the only type run
  --db DB             specify custom directory containing allele databases for porA/fetA typing
                      directory must contain database files: "FetA_VR.fas", "PorA_VR1.fas", "PorA_VR2.fas"
                      for Bexsero typing: "fHbp_peptide.fas", "NHBA_peptide.fas", "NadA_peptide.fas", "BASTalleles.txt"
  --printseq          save porA/fetA or BAST allele sequences to file (default=off)
  --updatedb          update allele database from <pubmlst.org>
  --test              run test example
  --version           show program's version number and exit
```


##Quick start

**To perform *in silico* serotyping on FASTA files:**

`$ meningotype <fasta1> <fasta2> <fasta3> ... <fastaN>`

The serotypes are printed in tab-separated format to `stdout`.

**To save results to a tab-separated text file, redirect `stdout`:**

`$ meningotype <fasta1> <fasta2> <fasta3> ... <fastaN>  > results.txt`

**To perform *in silico* serotyping AND finetyping of the porA and fetA genes:**

`$ meningotype --finetype <fasta1> <fasta2> <fasta3> ... <fastaN>`

**To save finetyping sequences of the alleles to a file (eg. for uploading "new" sequences to [http://pubmlst.org/neisseria/](http://pubmlst.org/neisseria/)):**

`$ meningotype --finetype --printseq <fasta1> <fasta2> <fasta3> ... <fastaN>`

##Updating the allele databases

**To update the allele databases from http://pubmlst.org/neisseria/ :**  
*Warning: Ensure you back up your old databases if you wish to keep them.*

	$ meningotype.py --updatedb

A copy of the old database is saved just in case, but is overwritten with each subsequent   ```--updatedb```.


##Citation

Please cite as:

Kwong JC, Gonçalves da Silva A, Stinear TP, Howden BP and Seemann T.  
***meningotype*: *in silico* typing for *Neisseria meningitidis*.**  
GitHub https://github.com/MDU-PHL/meningotype

##Bugs

Please submit via the [GitHub issues page](https://github.com/MDU-PHL/meningotype/issues).  

Note that the finetyping databases and website are curated and hosted by http://pubmlst.org/neisseria/. For issues with the databases, please contact the [pubmlst curator](mailto:keith.jolley@zoo.ox.ac.uk).

##Software Licence

[GPLv3](https://github.com/MDU-PHL/meningotype/blob/master/LICENSE)

##References

* Mothershed et al. J Clin Microbiol, 2004; 42(1): 320-328.
* Jolley et al. FEMS Microbiol Rev, 2007; 31: 89-96.
* Brehony et al. Vaccine, 2016; 34(39): 4690-4697.
* Bambini et al. PLoS One, 2013; 8(5): e65043.
* See also [http://www.neisseria.org/nm/typing/](http://www.neisseria.org/nm/typing/).
