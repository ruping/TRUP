TRUP is a pipeline designed for analyzing RNA-seq data from Tumor samples.


Learn More
---
...


Dependences
---
(i)   Samtools

(ii)  TopHat (version > 1.4)

(iii) GSNAP (version > 2012-01-11)

(iv)  Velvet (https://github.com/dzerbino/velvet)

(v)   Oases (https://github.com/dzerbino/oases).

(vi)  R libraries: feilds, KernSmooth, lattice, RColorBrewer and R2HTML.


Installation
---
No installation is required if pre-compiled binaries in bin/ are compatible with your system.

If pre-compiled binaries are incompatible with your system, you will need to rebuild them from source.

Make sure you have installed Bamtools (https://github.com/pezmaster31/bamtools). Then write down the bamtools_directory where ./lib/ and ./include/ sub-directories are located. Currently bamtools 2.0 has been tested.

run make in following way to install

	make BAMTOOLS_ROOT=/bamtools_directory/

The binaries will be built at bin/. Any existing files will be overwritten.


Usage
---

TRUP can be run with UNIX command-line interface.

Assuming that the two gzipped fastq-files are called as SAMPLE_1.fq.gz and SAMPLE_2.fq.gz (where "\_1" and "\_2" means mate 1 and mate 2 reads, respectively), the path to the directory of the pipeline is called "PD" and the path to tha annotation files are called "AD", one can run the pipeline in the following way:

Run-level (Quality check, mapping of spiked-in read pairs and tell you the insert size etc.):

	perl RTrace.pl --runlevel 1 --lanename SAMPLE --root PD --anno AD

If trimming is needed, add the following two options, ``--readlen RL --trimedlen TL``, where "RL" is the original read length, "TL" is the read length after trimming. If one you want to stop after the quality check, add ``--QC`` in the command call.\\



Options
---
...


Contact
---
Sun Ruping

Dept. Vingron (Computational Molecular Biology)
Max Planck Institute for Molecular Genetics. Ihnestr. 63-73, D-14195 Berlin, Germany

Email: ruping@molgen.mpg.de

Project Website: https://github.com/ruping/TRUP
