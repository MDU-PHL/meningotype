#meningotype

*In silico* serotyping and finetyping (*porA* and *fetA*) of *Neisseria meningitidis*  

##Authors

* Jason Kwong (@kwongjc)
* Anders Gonçalves da Silva
* Torsten Seemann (@torstenseemann)

##Dependencies

* [Python 2.7.x](https://www.python.org/)
* [BioPython](http://biopython.org/)
* [isPcr v33](http://hgwdev.cse.ucsc.edu/~kent/src/) by Jim Kent

##Installation

To install:

```
$ git clone https://github.com/MDU-PHL/meningotype
```

### To test installation

Once installed, you can run the following to ensure `meningotype` is successfully working:

    ```$ meningotype.py --test```

If everything works, you will see the following:

```
$ meningotype.py --test
Running meningotype.py on test examples ... 
$ meningotype.py A.fna B.fna C.fna W.fna X.fna Y.fna
SAMPLE_ID	SEROTYPE	PorA_TYPE	FetA_TYPE
test/A.fna	A
test/B.fna	B
test/C.fna	C
test/W.fna	W
test/X.fna	X
test/Y.fna	Y
```

or, to check finetyping:

```
$ meningotype.py --test --finetype
Running meningotype.py on test examples ... 
$ meningotype.py A.fna B.fna C.fna W.fna X.fna Y.fna
SAMPLE_ID	SEROTYPE	PorA_TYPE	FetA_TYPE
test/A.fna	A	7,13-1		F1-5
test/B.fna	B	5-2,10-1	F3-6
test/C.fna	C	21,26-2		F1-5
test/W.fna	W	5,2			F1-1
test/X.fna	X	5-1,10-1	F4-23
test/Y.fna	Y	5-2,10-1	F4-1
```


##Usage

```	$ meningotype.py -h

	usage: 
	  meningotype.py [OPTIONS] <fasta1> <fasta2> <fasta3> ... <fastaN>
	
	In silico typing for Neisseria meningitidis
	
	PCR Serotyping Ref: Mothershed et al, J Clin Microbiol 2004; 42(1): 320-328
	
	porA and fetA typing Ref: Jolley et al, FEMS Microbiol Rev 2007; 31: 89-96
	See also http://www.neisseria.org/nm/typing/
	
	positional arguments:
	  FASTA       input FASTA files eg. fasta1, fasta2, fasta3 ... fastaN
	
	optional arguments:
	  -h, --help  show this help message and exit
	  --finetype  Perform porA and fetA fine typing (default=off)
	  --db DB     Specify custom directory containing allele databases for porA/fetA typing
	              Directory must contain database files "FetA_VR.fas", "PorA_VR1.fas", and "PorA_VR2.fas"
	  --printseq  Save porA, porB and fetA allele sequences to file (default=off)
	  --updatedb  update allele database from <www.ng-mast.net>
	  --test      run test example
	  --version   show program's version number and exit
```

##Quick start

**To perform *in silico* serotyping on FASTA files:**

`$ meningotype.py <fasta1> <fasta2> <fasta3> ... <fastaN>`

The serotypes are printed in tab-separated format to `stdout`.

**To save results to a tab-separated text file, redirect `stdout`:**

`$ meningotype.py <fasta1> <fasta2> <fasta3> ... <fastaN>  > results.txt`

**To perform *in silico* serotyping AND finetyping of the porA and fetA genes:

`$ meningotype.py --finetype <fasta1> <fasta2> <fasta3> ... <fastaN>`

**To save finetyping sequences of the alleles to a file (eg. for uploading "new" sequences to [http://pubmlst.org/neisseria/](http://pubmlst.org/neisseria/)):**

`$ meningotype.py --finetype --printseq <fasta1> <fasta2> <fasta3> ... <fastaN>`

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

Note that the finetyping databases and website are curated and hosted by http://pubmlst.org/neisseria/. For issues with the databases, please contact the pubmlst curator](mailto:keith.jolley@zoo.ox.ac.uk).

##Software Licence

[GPLv3](https://github.com/MDU-PHL/meningotype/blob/master/LICENSE)

##References

* Mothershed et al. J Clin Microbiol, 2004; 42(1): 320-328.
* Jolley et al. FEMS Microbiol Rev, 2007; 31: 89-96.
* See also [http://www.neisseria.org/nm/typing/](http://www.neisseria.org/nm/typing/).
