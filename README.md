# ngshmmalign
David Seifert (david.seifert@bsse.ethz.ch)

## Introduction
In the current sequencing landscape, NGS reads are aligned using such aligners as **bwa** (http://bio-bwa.sourceforge.net) or **bowtie** (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). While these aligners are fast and work well for transcriptomic and exomic data of large eukaryotic genomes, they produce a number of artefacts on genomes experiencing a large number of mutations in the course of their evolution, such as HIV-1 and HCV. These RNA viruses show up as a heterogeneous mixture, with indels and point mutations causing problems in the alignment step. In order to produce more sensitive alignments, **ngshmmalign** implements a profile HMM to produce alignments that are more amenable to studies of such viruses, without resorting to a global, genome-wide multiple sequence alignment.


## Idea
The profile HMM is a well-known probabilistic graphical model, known for instance from HMMER (http://hmmer.org):
<p align="center">
	<img src="img/pHMM.png?raw=true" alt="Profile HMM graphical model"/>
</p>
Here the reference genome consists of the Match (blue) states. Match states can emit only a single base (at conserved loci) or multiple bases (at loci with SNVs). Insertions (red) model technical artefacts, whereas deletions (green) capture both technical artefacts and true biological indels. The flanking states (violet) capture technical problems, such as sequencing into the adapter regions on an Illumina sequencer, and will lead to the clipping of bases.


## Features
**ngshmmalign** currently performs an exhaustive glocal (global-to-local) alignment, which is very similar in nature to the Needleman-Wunsch algorithm. It has the following features:

- Takes either single-end or paired-end reads as input.
- Writes the alignment into a fully compliant SAM file.
- Produces a consensus reference containing ambiguous bases (e.g. A + G = R) and a reference where bases are determined by majority vote.
- Picks a _random_ optimal alignment. All currently known aligners pick an optimal alignment deterministically. While such a strategy makes sense for many situations, it might not be optimal for small viruses. Many viruses for instance have stretches of the same base, also known as homopolymers. Homopolymers can lead to problems, as the deletion probability increases disproportionately in such regions. Imagine having a homopolymeric stretch and some erroneous (lacking one base in the homopolymer) NGS reads. In such cases **bwa** produces the following alignment:
<p align="center">
	<img src="img/bwa_gaps.png?raw=true"/>
</p>
  which is suboptimal, as the frequencies of gaps at the first position of the homopolymer will likely be picked up by SNV callers downstream. In the case of **ngshmmalign** the alignment would look like
<p align="center">
	<img src="img/ngshmmalign_gaps.png?raw=true"/>
</p>
  which is much better, as the _technical_ error is **not** turned into a _systematic_ error, as in the case of all deterministic alignment algorithms.
- The profile HMM allows for introducing gaps due to the inhomogeneous Markov chain, which is not possible with other aligners, due to their costly genome-wide uniform gap open penalties. Assume we have two references, with one having lost a codon (i.e. an indel). We have two reads, both originating from reference 2. An alignment with **bwa** would yield
<p align="center">
	<img src="img/bwa_indel.png?raw=true" alt="Incorrect bwa alignment"/>
</p>
  which is incorrect, whereas an alignment with **ngshmmalign** would yield
<p align="center">
	<img src="img/ngshmmalign_indel.png?raw=true" alt="Correct ngshmmalign alignment"/>
</p>
  due to **ngshmmalign** being liberal on the gap-open penalties around the indel.
- Can filter likely invalid paired-end read configuration, for instance:
<p align="center">
	<img src="img/PE_Modes.png?raw=true" alt="Paired-end alignment outcomes"/>
</p>
  The latter three cases are likely to be a technical artefact and can lead to problems in downstream haplotype assembly.
- Fully parallelised using OpenMP.
- Has a template length cut-off filter. Templates that are too long are removed, using a 3 sigma cut-off threshold.


### Planned Features
- Parameter estimation using a local multiple sequence alignment. Currently, **ngshmmalign** takes a multiple sequence alignment in FASTA format as input reference.
- Lookup using FM index, in order to speed up alignment.
- Banded alignment, further improving speed.
- Writing the CIGAR using either `'M'` for aligned bases, or `'='` and `'X'` for alignment match and mismatch, respectively.
- Writing out `NM:i:` tag, the number of mismatches to the reference.


## Requirements
As **ngshmmalign** is still under heavy development, we will not be making release tarballs yet. If you still wish to give **ngshmmalign** a try, you'll need the following tools and libraries:

1.  A **C++11** compliant compiler. The **ngshmmalign** codebase makes extensive use of C++11 features. A recent GCC or Clang release should suffice.

    It also strongly recommended that you use an **OpenMP**-capable compiler, as you will not be able to utilise full parallelisation otherwise.

2.  **Autoconf**; latest 2.69 release (http://www.gnu.org/software/autoconf/)

    GNU Autoconf produces the ./configure script from configure.ac.

3.  **Automake**; latest 1.15 release (http://www.gnu.org/software/automake/)

    GNU Automake produces the Makefile.in precursor, that is processed with ./configure to yield the final Makefile.

4.  **Autoconf Archive**; latest 2016.03.20 release (http://www.gnu.org/software/autoconf-archive/)

    Our configure.ac requires a number of m4 macros from the Autoconf archive.

5.  **Boost**; latest 1.60 release (http://www.boost.org/)

    Boost provides the necessary abstractions for many different types.

6.  **standard Unix utilities**; such as sed, etc...

    If you cannot execute a command, chances are that you are missing one of the more common utilities we require in addition to the tools listed above.


## Building
If you have all dependencies satisfied, proceed by
```
./autogen.sh
./configure
make
```
You should now have a binary called `ngshmmalign` in the current build directory. You can either install this manually or call
```
make DESTDIR="${D}" install
```
where you specify the destination in `${D}`.