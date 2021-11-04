# NIQKI
NIQKI stand for Next index to query kmer intersection using HyperMinhash, NIQKI
is an alignment-free and reference-free software which allows to index
a large collection of similar genomes. NIQKI can also be applied to
reads from unassembled genomes.

Documentation
-------------

Usage: `niqki [options]`

1. Genomes list :

    `--list`, `-l <file>`             Use the content of the given (raw formatted)

    Example:

        niqki --list my_genomes.txt

2. Change k-mer and HyperMinhash size

    `--kmer`, `-k <int>`         Set the value of paramter k to the given value.

    `--w-size`, `-w <int>`       Set the value of parameter w to the given value.

    Example:

        niqki --kmer 30 or -k 30
        niqki --w-size 8 or -w 8


3. Optional usage :

    `--query`, `-Q <file>`        Search the given `<file>` in the index and print the genomes distance with this sequence.

    `--output`, `-o <pattern>`

    Examples:

        niqki --list my_genomes.txt --query my_genomes.txt
        niqki --list my_genomes.txt --output /tmp/outPutNiqki
        
4. Other usage:

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
