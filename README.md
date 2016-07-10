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

- Parameter estimation using multiple windows with local multiple sequence alignments. To this end, **ngshmmalign** uses **MAFFT** (http://mafft.cbrc.jp/alignment/software/) on a subsample of all reads to estimate the parameters of the profile HMM.
- **ngshmmalign** can also use a multiple sequence alignment (including just a single contig/sequence) in FASTA as an input reference.
- Takes either single-end or paired-end reads as input.
- Writes the alignment into a fully compliant SAM file, that passes Picard (http://broadinstitute.github.io/picard/) validation.
- Produces both a consensus reference containing ambiguous bases (e.g. A + G = R) and a reference where bases are determined by majority vote.
- Serializes the (inferred) profile HMM transition and base emission tables in order to reuse for future alignments.
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
- Writing the `NM:i:` tag, the number of mismatches (or edit distance) to the reference. This takes ambiguous bases into account, for instance
  
  ```
  Ref:   AARAA
  Read:   AAA
  ```
  
  gives an edit distance of 0, since `R` can be either an `A` or `G`, whereas
  
  ```
  Ref:   AARAA
  Read:   ATA
  ```
  
  gives an edit distance of 1, as `R` does not include a `T`.
- Writing the CIGAR using either `'M'` for aligned bases, or `'='` and `'X'` for alignment match and mismatch, respectively. When using `'X'` for alignment mismatches, the number of `'X'` is consistent with the `NM:i:` tag.
- Writing the `MD:Z:` tag correctly, allowing for reference-free analysis.
- Fully parallelised using OpenMP. Allows specifying a seed value for the random number generator, such that the alignment becomes deterministic and reproducible.
- Has a template length cut-off filter. Templates that are too long are removed, using a 3 sigma cut-off threshold.
- Includes a k-mer based index. In practice, on a Haswell-based 2.8 GHz i7-4558U a thread performance of 90-100 reads/thread/s can be achieved (this yields with 4 logical cores an overall performance of 350-400 reads/s) with a 9800 nt HIV-1 genome. Nevertheless, you can always opt out of the indexing, and always perform a globally optimal alignment, with a performance of about 5-6 reads/thread/s

### Planned Features
- Banded alignment, further improving speed.


## Requirements
As **ngshmmalign** is still under heavy development, we will not be making release tarballs yet. If you still wish to give **ngshmmalign** a try, you'll need the following tools and libraries:

1.  A **C++11** compliant compiler. The **ngshmmalign** codebase makes extensive use of C++11 features.

    GCC 4.8 and later have been verified to work, although we recommend you use at least GCC 5. Clang 3.7 and later have been verified and are also recommended, due to Clang introducing OpenMP with 3.7. Versions of Clang before 3.7 will not be able to utilise multi-threading.

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

7.  **MAFFT** (optional); (http://mafft.cbrc.jp/alignment/software/)

    If you wish to align reads and optimize the reference sequence concurrently. MAFFT is used to align a subsample of reads in order to estimate biological indels.

### OS X
We strongly recommend you use MacPorts (http://www.macports.org) to install dependencies. We also recommend you employ Clang from MacPorts, as it is the only OpenMP-capable compiler that is simultaneously ABI-compatible with installed libraries, such as boost. While building with GCC on OS X is possible, it requires an orthogonal toolchain which is far more involved and beyond the scope of this README.

### GNU/Linux
On a GNU/Linux system, the aforementioned recommendations are reversed. Most GNU/Linux distributions are built using GCC/libstdc++, which as of GCC 5.1 is not backwards compatible with Clang, and as such building with Clang produced object files will fail in the final linking step.


## Building
1.  If you have all dependencies satisfied, proceed by generating the build system files
    ```
    ./autogen.sh
    ```

2.  Then, run the configure script. On GNU/Linux, you would do
    ```
    ./configure
    ```
    whereas on OS X, you would also need to specify the OpenMP-capable C++ compiler
    ```
    ./configure CXX=clang++-mp-3.7
    ```
	for instance, if you installed Clang 3.7 from MacPorts.

3.  Then, compile the sources using
    ```
    make -j2
    ```

4.  You should now have a binary called `ngshmmalign` in the current build directory. You can either install this manually or call
    ```
    make DESTDIR="${D}" install
    ```
    where you specify the destination in `${D}`.


## Running
The parameters of **ngshmmalign** can be viewed with the help option `-h`. If you wish to use **MAFFT**, you have two options:

1.  Specify the path of the `mafft` program in the environmental variable `MAFFT_BIN`.

2.  As a fallback, if you ensure that `mafft` is located in your `$PATH`, it will also be found.
