Copyright (C) 2010-2020, Marx Gomes van der Linden

# MDTools 1.1

MDTools was developed mainly as a tool to generate input files for HmmPred (see references), but some of its features can be useful to accomplish more general tasks when working with protein structure files in the PDB format. 

In particular, MDTools it can be used to generate correct amino acid sequences or burial sequences from PDB files. If there are gaps in the original structure, those gaps will be represented by the appropriate number of Xs in the generated sequence.

## Prerequisites

MDTools was compiled with GCC (g++) and tested in a GNU/Linux
environment, but it should probably work just fine in any common
operating system with a standard C++ compiler. To compile MDTools,
you will need the "program_options" library from Boost, easily available as free software in most Linux distributions.

## Compiling

The easiest way to compile MDTools is simply to run:

    g++ ClashFinder.cpp DatFile.cpp rcntr.cpp Protein.cpp MDTools.cpp Atom.cpp -O3 -lboost_program_options -oMDTools


## Usage

To view all options, just run MDTools without any arguments.

To generate the amino acid sequence for a PDB file, use:

    MDTools seq -fpdb PDBFILE

To generate a sequence of burial levels in 4 layers, use:

    MDTools seq -fpdb -tburial -l4 PDBFILE


## References

MDTools and HmmPred are used in the research published on the following papers:

- ROCHA, J. R. ; VAN DER LINDEN, M. G. ; FERREIRA, D. C. ; AZEVEDO, P. H. ; PEREIRA DE ARAUJO, A. F. . Information-theoretic analysis and prediction of protein atomic burials: on the search for an informational intermediate between sequence and structure. Bioinformatics (Oxford. Print), v. 28, p. 2755-2762, 2012. 

- VAN DER LINDEN, MARX GOMES; FERREIRA, DIOGO CÉSAR ; DE OLIVEIRA, LEANDRO CRISTANTE ; ONUCHIC, JOSÉ N. ; DE ARAÚJO, ANTÔNIO F. PEREIRA . Ab initio protein folding simulations using atomic burials as informational intermediates between sequence and structure. Proteins (Print), v. 82, p. n/a-n/a, 2013. 

- FERREIRA, DIOGO C. ; VAN DER LINDEN, MARX G. ; DE OLIVEIRA, LEANDRO C. ; ONUCHIC, JOSÉ N. ; PEREIRA DE ARAÚJO, ANTÔNIO F. . Information and redundancy in the burial folding code of globular proteins within a wide range of shapes and sizes. Proteins (Print), v. 84, p. 515-531, 2016. 
