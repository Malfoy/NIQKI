# NIQKI
NIQKI stand for Next Index to Query _K_-mer Intersection.

NIQKI is an sketch based software (similar to Mash or Dashing) which can index the largest sequence collections. 

Once the sketch index built, NIQKI can compare  query sequences to indexed sequences and output matches above given threshold.

Using inversed indexes and well designed fingerprint NIQKI can be order of magnitudes faster on large instances than concurent approaches with comparable precision and memory usage.

Documentation
-------------

Usage: `niqki [options]`

1. Indexing  :

`--index`, `-I <file>`             Input file of files to Index.

Example:

    niqki --index my_genomes.txt
    
`--indexlines,`, `-i <file>`             Input fa/fq where each line is a separate entry to Index.

Example:

    niqki --indexlines my_genomes.fa
        
2. Query  :

`--query,`, `-Q <file>`             Input file of files to Query.

Example:

    niqki --query my_genomes.txt

`--querylines,`, `-q <file>`             Input fa/fq where each line is a separate entry to Query.

Example:

    niqki --indexlines my_genomes.fa

3. Change sketch parameters

`--kmer`, `-K <int>`         Kmer size (31).
    
`--sketch,`, `-S <int>`        Set sketch size to 2^S (15).



4. Change sketch advanced parameters

`--word`, `-W <int>`                  Fingerprint size (12). Modify with caution, larger fingerprints enable queries with less false positive but increase
                                    EXPONENTIALY the overhead as the index count S*2^W cells. 
    
`--Genomes_sizes`, `-G <int>`         Rought expectation of the genome sizes.

`--HHL`, `-H <int>`                   Size of the hyperloglog section (4).  Modify with caution and prefer to use -G.

5. Change output

`--output`, `-O <filename>`           Output file (niqkiOutput.gz)

`--minjac`, `-J <int> `               Minimal jaccard Index to report (0.1).
                                    
`--pretty`, `-P`                       Print a human-readable outfile. By default the outfile is in binary.
        
6. Other usage:

`--dumpv`, `-D <filename>`             Dump the current index to the given file.
  
`--load`, `-L <filename>`             Load an index to the given file.


`--indexdownload`, `-Iddl <filename>` Get a list of NCBI accesion to download and to put it in the index (experimental). This this post to get such a list  https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#allcomplete
 
`--matrix`  , `-M <filename>`            Index the inout file of file, query the index against itself and output the corresponging matrix distance


`--logo`                           Print ASCII art logo, then exit.

`--help`, `-h`                    Print usage and exit.


Installation
------------

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

First, clone the `NIQKI` repository:
```sh
git clone https://github.com/Malfoy/NIQKI.git
```

Once cloned, go to the newly created directory and artificially
restore the relative order of creation/modification dates for some
files (see explanation in previous section).

```sh
cd NIQKI
```

Now, run the `configure` script, build and install.
```sh
./configure --prefix=$HOME/local_install
make
make install
```


#### Uninstall

To remove `NIQKI`from your system use the following command:

```sh
cd NIQKI && make uninstall
```


Auteurs/Authors:
----------------

* Clément AGRET     <clement.agret@lirmm.fr>
* Bastien CAZAUX    <bastien.cazaux@univ-lille.fr>
* Antoine LIMASSET  <antoine.limasset@univ-lille.fr>

Programmeurs/Programmers:
-------------------------

* Clément AGRET     <clement.agret@lirmm.fr>
* Antoine LIMASSET  <antoine.limasset@univ-lille.fr>
