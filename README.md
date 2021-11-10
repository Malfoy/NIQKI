# NIQKI
NIQKI stand for Next Index to Query K-mer Intersection.

NIQKI is an sketch based software, similar to [Mash](https://github.com/marbl/Mash) or [Dashing](https://github.com/dnbaker/dashing), which can index the largest sequence collections. 

Once the sketch index built, NIQKI can compare  query sequences to indexed sequences and output matches above given threshold.

Using inversed indexes and well designed fingerprint NIQKI can be order of magnitudes faster on large instances than concurent approaches with comparable precision and memory usage.


# Table of contents
1. [Command line options](#cmo)
2. [Installation](#install)
3. [Sample usages](#su)
4. [Credits](#cr)

## Command line options <a name="cmo"></a>

### 1. Indexing  :

NIQKI can handle FASTA, MultiLine FASTA  and FASTQ files gziped or not.

`--index`, `-I <file>`             Input file of files to Index. 

Such file of files should contain one file name or path per line. Each  file will be sketched and inserted in the index as a single entry.

Example:

    niqki --index my_genomes.txt
    
`--indexlines,`, `-i <file>`             Input file where each line is a separate entry to Index.

Each line of this file representing a DNA sequence will be sketched and inserted as a single entry.

Example:

    niqki --indexlines my_genomes.fa
        
### 2. Query  :

`--query,`, `-Q <file>`             Input file of files to Query.

Example:

    niqki --query my_genomes.txt

`--querylines,`, `-q <file>`             Input file where each line is a separate entry to Query.

Example:

    niqki --indexlines my_genomes.fa

### 3. Change sketch parameters

`--kmer`, `-K <int>`         Kmer size (Default value is 31).

Each sequence will be splited into words of length K before to be hashed and potentially inserted to the sequence's sketch.
    
`--sketch,`, `-S <int>`        Set sketch size to 2^S (Default value is 15 meaning 32,768 fingerprints per sketch).

Each sequence will be summarized by a sketch of 2^S fingerprints. Larger sketches improve the estimation accurary but raise linearly the running time and the used memory. Incrementing S may therefore double running time and memory usage. In practice you can expect each sketch to allocate 2^S times 4 Bytes (meaning 16 kilobytes with default value).

### 4. Change sketch advanced parameters

`--word`, `-W <int>`                  Fingerprint size (Default value is  12). 

Modify with caution, larger fingerprints decrease the  false positive  hits rate but increase EXPONENTIALY the overhead as the index intialization allocates S*2^W vectors. 
    
`--Genomes_sizes`, `-G <int>`         Rought expectation of the genome sizes.

Using the expected genome size will allow the index to choose the best Hyperloglog fingerprint size (provinding the lowest false positive rate) by estimating the subsampling rate.

`--HHL`, `-H <int>`                   Size of the hyperloglog fingerprint (4).  Modify with caution and prefer to use -G.

### 5. Change output

`--output`, `-O <filename>`           Output file (Default value is niqkiOutput.gz).

In wich file the hits should be reported.

`--minjac`, `-J <float> `               Minimal jaccard Index to report (Default value is 0.1).

Only the hits with a estimated Jaccard Index above the threshold will be wrote in the output.
                                    
`--pretty`, `-P`                       Print a human-readable outfile.

Using this option the output look like this
```
Query_genome1.fa
Indexed_genome5.fa 31654
Indexed_genome7.fa 20543
Indexed_genome3.fa 10654
Query_genome2.fa
Indexed_genome5.fa 17654
Indexed_genome7.fa 2543
```

Otherwise the output will be displayed in binary following the pattern
```
Query_genome1.fa[ENDL][FOUR BYTES indicating the number of hit to parse]
(Foreach hit)
[FOUR BYTES indicating the indexed genome identifier][FOUR BYTES indicating the number of hits found]
Query_genome2.fa[ENDL]...
```
### 6. Other usage:

`--dump`, `-D <filename>`             Dump the current index to the given file.
  
`--load`, `-L <filename>`             Load an index to the given file.


`--indexdownload`, `-Iddl <filename>` Get a list of NCBI accesion to download and to put it in the index (experimental). 
This this post to get such a list  https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#allcomplete
 
`--matrix`  , `-M <filename>`            Index the input file of file, query the index against itself and output the corresponging matrix distance.
The matrix distance will look like this
```
##Names	ecoli01p.fa.gz	ecoli03p.fa.gz	ecoli07p.fa.gz	ecoli04p.fa.gz	ecoli05p.fa.gz	ecoli02p.fa.gz	ecoli06p.fa.gz	ecoli08p.fa.gz	ecoli09p.fa.gz	ecoli10p.fa.gz	
ecoli01p.fa.gz	1	0.938019	0.836639	0.910797	0.884796	0.967773	0.861206	0.814117	0.791321	0.769226	
ecoli03p.fa.gz	0.938019	1	0.888824	0.970367	0.941711	0.968567	0.915741	0.864014	0.8396	0.815033	
ecoli07p.fa.gz	0.836639	0.888824	1	0.914246	0.941223	0.862457	0.969086	0.970581	0.942017	0.91275	
ecoli04p.fa.gz	0.910797	0.970367	0.914246	1	0.970093	0.940094	0.942688	0.888702	0.863403	0.838287	
ecoli05p.fa.gz	0.884796	0.941711	0.941223	0.970093	1	0.913025	0.971161	0.914337	0.888031	0.861877	
ecoli02p.fa.gz	0.967773	0.968567	0.862457	0.940094	0.913025	1	0.888062	0.838806	0.815247	0.791534	
ecoli06p.fa.gz	0.861206	0.915741	0.969086	0.942688	0.971161	0.888062	1	0.94104	0.913818	0.886292	
ecoli08p.fa.gz	0.814117	0.864014	0.970581	0.888702	0.914337	0.838806	0.94104	1	0.970001	0.939148	
ecoli09p.fa.gz	0.791321	0.8396	0.942017	0.863403	0.888031	0.815247	0.913818	0.970001	1	0.968262	
ecoli10p.fa.gz	0.769226	0.815033	0.91275	0.838287	0.861877	0.791534	0.886292	0.939148	0.968262	1
```

`--logo`                           Print ASCII art logo, then exit.

`--help`, `-h`                    Print usage and exit.


## Installation <a name="install"></a>

### Requirements

NIQKI requires:

* A modern, C++11 ready compiler such as `g++` version 4.9 or higher or `clang` version 3.2 or higher.
* A 64-bit operating system. Either Mac OS X or Linux are currently supported.
* `zlib` to be already installed on your system (on Debian/Ubuntu it corresponds to the packages `zlib1g-dev`).

### Single user installation

To download and install `NIQKI` into some
user local directory (e.g., `${HOME}/local_install`) , use the
following commands:


#### NIQKI

```sh
git clone https://github.com/Malfoy/NIQKI.git ;

cd NIQKI ;

./configure --prefix=$HOME/local_install ;

make ;

make install ;
```


#### Uninstall

To remove `NIQKI`from your system use the following command:

```sh
cd NIQKI && make uninstall
```


## Sample Usages <a name="su"></a>

To test your installation you can go to the resources folder
```sh
cd resources
```


You can generate a distance matrix from the provided file of files.

```sh
niqki -M file_of_file.txt -O matrix.gz
```

And see it with zcat (or gzip -d)

```sh
zcat matrix.gz 
```
By default Niqki generate a hit list for each query file
 
 ```sh
niqki -I file_of_file.txt -Q file_of_file.txt -O hits.gz -P
zcat hits.gz
```

You can change the sketch size depending of our usage (here 1,024 fingerprint per sketch)
 ```sh
niqki -I file_of_file.txt -Q file_of_file.txt -O hitsS10.gz -P -S 10
zcat hitsS10.gz
```

If you have an idea of your mean genomes size you should indicate it to improve the false positive rate
 ```sh
niqki -I file_of_file.txt -Q file_of_file.txt -O hits.gz -P -S 10 -G 5000000
zcat hits.gz
```

You can also tweak the minimal fraction of shared fingerprint necesary to output hits.
You can report all hits with -J 0
 ```sh
niqki -I file_of_file.txt -Q file_of_file.txt -O allhits.gz -P -S 10 -G 5000000 -J 0
zcat allhits.gz
```
Or only hits with 95% shared fingerprints with -J 0.95
 ```sh
niqki -I file_of_file.txt -Q file_of_file.txt -O hits95.gz -P -S 10 -G 5000000 -J 0.95
zcat hits95.gz
```

## Credits <a name="cr"></a>

For further information our preprint "Toward optimal fingerprint indexing for large scale genomics"
 can be consulted in [Biorxiv](https://www.biorxiv.org/content/10.1101/2021.11.04.467355v1)



Authors:
----------------

* Clément AGRET     <clement.agret@lirmm.fr>
* Bastien CAZAUX    <bastien.cazaux@univ-lille.fr>
* Antoine LIMASSET  <antoine.limasset@univ-lille.fr>

Programmers:
-------------------------

* Clément AGRET     <clement.agret@lirmm.fr>
* Antoine LIMASSET  <antoine.limasset@univ-lille.fr>
