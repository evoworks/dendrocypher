# README #

This README is for DendroCypher.  

The primary use of DendroCypher is to assist the process of identify and labelling specific nodes of a phylogenetic tree.  Hypothesis testing in an evolutionary framework often relies on models where known biological events (e.g., gene duplications, niche colonization, LGT) are treated as fixed-effects.  Thus a specific node, or nodes, must be identified and labelled in a way that can be read by software such as CODEML, CODEML_SBA and PROTEUS. The standard text representation (newick format) can be challenging to label by hand. DendroCypher is a simple program to aid this process.

### What is this repository for? ###

DendroCypher is written in C++.  It is based on tree and matrix classes originally developed for the program PROTEUS.

DendroCypher was developed to support a Protocols Unit within the *Current Protocols in Bioinformatics*. The paper describing the DendroCypher:

>*Bielawski, J. P., Baker, J. L., & Mingrone, J. (2016). Inference of episodic changes in natural selection acting on protein coding sequences via CODEML. Current protocols in bioinformatics, 54(1), 6-15.*

Please cite this paper if you use this repository.

This is version 0 of DendroCypher; please be patient as we work out the kinks and make updates as needed.

### How do I get set up? ###

For connivence, a pre-compiled version is supplied for Mac OS X.  This version might not work on all systems.

The program can be compiled from source:  make -f Makefile
