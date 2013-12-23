[metaHGT]
================================================

Real-time in situ Horizontal Gene Transfers (HGTs) detector using time series metagenomes

- Author: Chengwei Luo (luo.chengwei@gatech.edu; cluo@braodinstitute.org);
- Copyright: Chengwei Luo, Konstantinidis Lab, Georgia Institute of Technology, 2013;
- Citation: (we are working on it!!)

[What does it do]
================================================

metaHGT is a breakthrough algorithm designed for HGTs detection using time series metagenomes. Unlike most HGT inferring algorithms, the breakthrough of metaHGT is that it's not designed to infer what had happened in history, but rather, it aims to tell us what's currently happening in real time. More specifically, the current algorithms for HGT detection are based on comparing two sequences (full genome, draft genome, or partial genome sequences), and detecting the incongruence between genomic signal (represent the species evolution) and the target segment's signal (represent the evolution of that segment; in most cases, genes). The limits of such algorithms are obvious: the only signals that would be detected are those that were fixed in the genome by evolution, which usually would date back to tens of thousands of years ago. However, in nowadays microbiology, HGTs detection is essential in many regards including disease control and prevention, environmental protection, natinal security, and law enforcement. All of the above-mentioned require a real time and in-situ predictor that could monitor and eventually predict the genetic flow across species in the microbial world. On the other hand, though much efforts has been carried out in the laboratories to characterize and to quantify HGT, the scientific community still needs to investigate some fundamental aspects of HGT in nature such as the rate, the vehicles, and the dynamics of HGT. MetaHGT is, as far as we know, the first algorithm that would enable real-time and in situ detection of HGTs.

[OK, sounds awesome, what do I need in order to infer real time HGTs?]
================================================

metaHGT requires time-series metagenomes, which means multiple sampling points at different times on the same community. This being said, we do not offer advice on how long the whole longitudinal series should be, neither could we comment on how dense the sampling points should be. This is still a pioneering project and you just have to try to know. As of our experience with the Lanier Project, we sampled the freshwater lake at roughly monthly basis for three years, each sample had about 3-5 Gbp Illumina paired-end (GA-II and HiSeq) reads, and were able to track HGT between populations. 

In terms of selecting sequencing technologies, currently metaHGT requires paired-end reads, preferably Illumina; but in theory, other sequencing platforms such as Ion Torrent and AB SOLid would work as well. As other sequencing technologies that would produce significantly longer sequences (e.g., Oxford Nanopore's platforms) are getting huge momentum and will likely to be proven to be a breakthrough in DNA sequencing, metaHGT will evolve to use also longer single reads to better resolve HGTs.

You will need to assemble your metagenomes into contigs, and bin them into population bins. How to do that is out of the scope of metaHGT, but you can in general use the hybrid protocol presented in my previous work (Luo et al, ISME J, 2012; Luo et al, PLoS ONE, 2012) to assemble a metagenome, and you can use BinGeR (Luo et al, in prep; https://github.com/luo-chengwei/BinGeR) to bin those contigs into population bins. After that you can run metaHGT by following the next section.

[How to install]
================================================

Install is simple, first you clone the repository to your local machine by:

$ github clone git://github.com/luo-chengwei/metaHGT

then you should "cd" to the local clone and run:

$ python setup.py install

You can alway supply --prefix to alter the installation directory:

$ python setup.py install --prefix = user/defined/installation/directory

You should be all set after this.

There are some libraries that are needed for running metaHGT:

- Python 2.7
- Pysam >= 0.6 
- BioPython >= 1.58
- Numpy >= 1.5.1 
- Scipy >= 0.12.0
- BWA >= 0.6 (http://bio-bwa.sourceforge.net)
- Samtools >= 0.1.18 (http://samtools.sourceforge.net)

NOTE: BWA and Samtools might not be needed if you have already produced all required indexed BAM files (detailed in the next section).

Packages with older versions might work, but not tested. You can certainly manually install all of the packages, however, I recommend using Anaconda, which is a Python distribution for big data processing, it bundles many packages for data scientists. For more information, go to:

https://store.continuum.io/cshop/anaconda/

Anaconda includes all the dependencies to run metaHGT.

[How to run it]
================================================

** Input:

   You will absolutely need:

   - list of samples. They can be listed in a file in which each sample occupy a line, and replicates can be on the same line, separated by comma; you can also supply this list by typing them in the command line for the -l option;

   - all assembled contig files in fastA format.
  
  Then you will need either:
   
- all sorted and indexed BAM mapping files put in one directory. For instance, if you have 3 samples: A, B, and C, in longitudinal order; you will need cross-mapping and self-mapping .bam and .bai files of the adjacent samples (in this case: A-B, B-C pairs). For instance, for the A-B pair, you will need readA.vs.assemblyA.bam, readA.vs.assemblyB.bam, readB.vs.assemblyB.bam, readB.vs.assemblyA.bam; and all the associated .bai files from the indexing by Samtools.
  
  OR, if you don't have them ready, you will need:
   - reads of each sample in the one directory. 
   - BWA binary
   - Samtools binary

** Usage:

  You can run it as simple as:
  
  $ metaHGT.py -l A,B,C -o <out_dir> -b <bam_dir>

  if you don't have BAM files ready:
  
  $ metaHGT.py -l A,B,C -o <out_dir> -r <reads_dir> --bwa <BWA_binary> --samtools <Samtools_binary>

  Certainly, if you have BWA and Samtools already in your $ENV, and named as "bwa" and "samtools", respectively, you don't need to specify --bwa or --samtools. You can simply check if you have them ready by typing:

$ bwa

and:

$ samtools

For detailed usage, please run:

$ metaHGT.py --help

This will print out a detailed usage message.

[Output files]
================================================

From running metaHGT.py, an output file will be generated (specified by -o, or by default: ./HGTs_info.txt)

- HGTs_info.txt is the general information about HGTs inferred by metaHGT. It is a tab-delimited text file. Each line is a record, which contains the following field:
  
  1. time point 1 (sample ID 1);
  2. time point 2 (sample ID 2);
  3. bin ID 1;
  4. contig A from bin 1;
  5. contig A breakpoint location;
  6. contig A breakpoint orientation ('>' or '<');
  7. bin ID 2;
  8. contig B from bin2;
  9. contig B breakpoint location;
  10. contig B breakpoint orientation ('>' or '<');
  11. percentage of bin 1 involved;
  12. percentage of bin 2 involved;
  13. raw p-value;
  14. FDR-corrected p-value using Benjamini-Hochberg method.

  Note on orientation: '>' means that the fragment to the left of breakpoint of the current contig is involved in HGT. For instance, if a sequence can be presented as '++++++++++++++^--------------------', where '^' denotes the breakpoint. With '>' orientation, it means the '++++++++' part was linked to some sequences in other populations by HGT; the opposite ('-------' part was involved) is true with a '<' orientation.



