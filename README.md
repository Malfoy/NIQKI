# NIHM
NIHM stand for Next index to query kmer intersection using HyperMinhash, NiHM
is an alignment-free and reference-free software which allows to index
a large collection of similar genomes. NiHM can also be applied to
reads from unassembled genomes.

Documentation
-------------

Usage: `nihm [options]`

1. Genomes list :

    `--list`, `-l <file>`             Use the content of the given (raw formatted)

    Example:

        nihm --list my_genomes.txt

2. Change k-mer and HyperMinhash size

    `--kmer`, `-k <int>`         Set the value of paramter k to the given value.

    `--w-size`, `-w <int>`       Set the value of parameter w to the given value.

    Example:

        nihm --kmer 30 or -k 30
        nihm --w-size 8 or -w 8


3. Optional usage :

    `--query`, `-Q <file>`        Search the given `<file>` in the index and print the genomes distance with this sequence.

    `--output`, `-o <pattern>`

    Examples:

        nihm --list my_genomes.txt --query my_genomes.txt
        redoak --list my_genomes.txt --output /tmp/outPutNihm

4. Other usage:

    `--help`, `-h`                    Print usage and exit.


Installation
------------

### Requirements

NiHM requires:

* A modern, C++11 ready compiler such as `g++` version 4.9 or higher or `clang` version 3.2 or higher.
* A 64-bit operating system. Either Mac OS X or Linux are currently supported.
* TODO

### Single user installation

To download and install `NiHM` into some
user local directory (e.g., `${HOME}/local_install`) , use the
following commands:


#### NiHM

First, clone the `NiHM` repository:
```sh
git clone https://github.com/Malfoy/NIHM.git
```

Once cloned, go to the newly created directory and artificially
restore the relative order of creation/modification dates for some
files (see explanation in previous section).

```sh
cd NIHM
touch configure.ac aclocal.m4 Makefile.am */Makefile.am
touch configure Makefile.in */Makefile.in
```

Now, run the `configure` script, build and install.
```sh
./configure --prefix=$HOME/local_install
make
make install
```


#### Uninstall

To remove `NIHM`from your system use the following command:

```sh
cd NIHM && make uninstall
```


Auteurs/Authors:
----------------

* Clément AGRET     <clement.agret@lirmm.fr>
* Bastien CAZAUX    <bastien.casaux@univ-lille.fr>
* Antoine LIMASSET  <antoine.limasset@univ-lille.fr>

Programmeurs/Programmers:
-------------------------

* Clément AGRET     <clement.agret@lirmm.fr>
* Antoine LIMASSET  <antoine.limasset@univ-lille.fr>
