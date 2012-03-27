TRUP is a pipeline designed for analyzing RNA-seq data from Tumor samples.


Learn More
---
...


Dependences
---
+   Samtools
+   TopHat (version > 1.4)
+   GSNAP (version > 2012-01-11)
+   Velvet (https://github.com/dzerbino/velvet)
+   Oases (https://github.com/dzerbino/oases).
+   R libraries: feilds, KernSmooth, lattice, RColorBrewer and R2HTML.


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

**Run-level 1** (Quality check, mapping of spiked-in read pairs and tell you the insert size etc.):

	perl RTrace.pl --runlevel 1 --lanename SAMPLE --root PD --anno AD

If trimming is needed, add the following two options, ``--readlen RL --trimedlen TL``, where "RL" is the original read length, "TL" is the read length after trimming. If one you want to stop after the quality check, add ``--QC`` in the command call.

**Run-level 2** (mapping with Tophat, generating reports):

	perl RTrace.pl --runlevel 2 --lanename SAMPLE --$root PD --anno AD --readlen RL --trimedlen TL --seglen SL --threads TH

where "SL" is the segment length used for Tophat, default is 25. "TH" is the number of computing threads. If the reads were trimmed in the Run-lever 1, one should keep the options ``--readlen RL --trimedlen TL``. One can add ``--$WIG`` to generate bigWiggle file for coverage visualization purpose.

**Run-level 3** (select ARP: abnormal read pairs based on second mapping using GSNAP, ARP are fed to denovo assembly):

	perl RTrace.pl --runlevel 3 --lanename SAMPLE --root PD --anno AD --threads TH --SM

where ``--SM`` means second mapping using GSNAP.

**Run-level 4** (Dectecting fusion events using the assembled transcripts):

	perl RTrace.pl --runlevel 4 --lanename SAMPLE --root PD --anno AD --threads TH

If the reads were trimmed in the Run-lever 1, one should keep the options ``--$readlen RL --trimedlen TL``.


The four Run-levels metioned above can also be combined in the same command line, as following,

	perl RTrace.pl --runlevel 1-4 --lanename SAMPLE --root PD --anno AD --readlen RL --trimedlen TL --seglen SL --threads TH --SM

One can also call combinations of different Run-levels, such as ``--$runlevel 1,2`` or ``--runlevel 2,3,4``. But in these combined ways, one should specify all necessary options in the command line.


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
