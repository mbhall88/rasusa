---
title: 'Rasusa: Randomly subsample sequencing reads to a specified coverage'
tags:
  - Rust
  - bioinformatics
  - genomics
  - fastq
  - fasta
  - subsampling
  - random
authors:
  - name: Michael B. Hall
    orcid: 0000-0003-3683-6208
    affiliation: 1
affiliations:
 - name: European Molecular Biology Laboratory, European Bioinformatics Institute EMBL-EBI, Hinxton, UK
   index: 1
date: 18 October 2021
bibliography: paper.bib
---

# Summary

A fundamental requirement for many applications in genomics is the sequencing of genetic
material (DNA/RNA). Different sequencing technologies exist, but all aim to accurately
reproduce the sequence of nucleotides (the individual units of DNA and RNA) in the
genetic material under investigation. The result of such efforts is a text file
containing the individual fragments of genetic material - termed "reads" - represented
as strings of letters (A, C, G, and T/U).

The amount of data in one of these read files depends on how much genetic material was
present and how long the sequencing device was operated. Read depth (coverage) is a
measure of the volume of genetic data contained in a read file. For example, coverage of
5x indicates that, on average, each nucleotide in the original genetic material is represented five
times in the read file.

Many of the computational methods employed in genomics are affected by coverage;
counterintuitively, more is not always better. For example, because sequencing devices
are not perfect, reads inevitably contain errors. As such, higher coverage increases the
number of errors and potentially makes them look like alternative sequences.
Furthermore, for some applications, too much coverage can cause a degradation in
computational performance via increased runtimes or memory usage.

We present Rasusa, a software program that randomly subsamples a given read file to a
specified coverage. Rasusa is written in the Rust programming language and is much
faster than current solutions for subsampling read files. In addition, it provides an
ergonomic command-line interface and allows users to specify a desired coverage or a
target number of nucleotides.


# Statement of need

Read subsampling is a useful mechanism for creating artificial datasets, allowing
exploration of a computational method's performance as data becomes more scarce. In
addition, the coverage of a sample can have a significant impact on a variety of
computational methods, such as RNA-seq [@Baccarella2018], taxonomic classification
[@Gweon2019], antimicrobial resistance detection [@Gweon2019], and genome assembly
[@Maio2019] - to name a few.

There is limited available software for subsampling read files. Assumably, most
researchers use custom scripts for this purpose. However, two existing programs for
subsampling are Filtlong [@filtlong] and Seqtk [@seqtk]. Unfortunately, neither of these
tools provides subsampling to a specified coverage "out of the box".

Filtlong is technically a filtering tool, not a subsampling one. It scores each read
based on its length and quality and outputs the highest-scoring subset. Additionally,
minimum and maximum read lengths can be specified, along with the size of the subset
required. Ultimately, the subset produced by Filtlong is not necessarily representative
of the original reads but is biased towards those with the greatest length or quality.
While this may sound like a good thing, in some applications, such as genome assembly,
it has been shown that a random subsample produces superior results to a filtered subset
[@Maio2019].

Seqtk does do random subsampling via the sample subcommand. However, the only option
available is to specify the number of reads required. Thus, it is up to the user to
determine the number of reads required to reach the desired coverage. While this serves
for Illumina sequencing data, which generally have uniform(ish) read lengths, it does
not work for other modalities like PacBio and Nanopore, where read lengths vary
significantly.

Rasusa provides a random subsample of a read file (FASTA or FASTQ), with two ways of
specifying the size of the subset. One method takes a genome size and the desired
coverage, while the other takes a target number of bases (nucleotides). In the genome
size and coverage option, we multiply the genome size by the coverage to obtain the
target number of bases for the subset. As such, the resulting read file will have, on 
average, the amount of coverage requested. In addition, Rasusa allows setting a random seed
to allow reproducible subsampling. Other features include user control over whether the
output is compressed and specifying the compression algorithm and level.

Rasusa is 21 and 1.2 times faster than Filtlong and Seqtk, respectively.

# Availability

Rasusa is open-source and available under an MIT license at
https://github.com/mbhall88/rasusa.

# Acknowledgements

We acknowledge contributions from Pierre Marijon and suggestions from Zamin Iqbal. In
addition, MBH is funded by the EMBL International PhD Programme.

# References

