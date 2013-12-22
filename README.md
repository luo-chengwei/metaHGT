[metaHGT]

Real-time in situ Horizontal Gene Transfers (HGTs) detector using time series metagenomes

Author: Chengwei Luo (luo.chengwei@gatech.edu; cluo@braodinstitute.org)
Copyright: Chengwei Luo, Konstantinidis Lab, Georgia Institute of Technology, 2013
Citation: (we are working on it!!)

[What does it do]

metaHGT is a breakthrough algorithm designed for HGTs detection using time series metagenomes. Unlike most HGT inferring algorithms, the breakthrough of metaHGT is that it's not designed to infer what had happened in history, but rather, in real time. More specifically, the current algorithms for HGT detection are based on comparing two sequences (full genome, draft genome, or partial genome sequences), and detecting the incongruence between genomic signal (represent the species evolution) and the target segment's signal (represent the evolution of that segment; in most cases, genes). The limits of such algorithms are obvious: the only signals that would be detected are those that were fixed in the genome by evolution, which usually would date to tens of thousands of years. However, in nowadays microbiology, HGTs detection is essential in many regards including disease control and prevention, environmental protection, natinal security, and law enforcement. All of the above-mentioned require a real time and in-situ predictor that could monitor and eventually predict the genetic flow across species in the microbial world. On the other hand, though much efforts has been carried out in the laboratories to characterize and to quantify HGT, the scientific community still needs to investigate some fundamental aspects of HGT in nature such as the rate, the vehicles, and the dynamics of HGT. MetaHGT is, as far as we know, the first algorithm that would enable real-time and in situ detection of HGTs.

[OK, sounds awesome, what do I need in order to infer real time HGTs?]

metaHGT requires time-series metagenomes, which means multiple sampling points at different time of the same community. This being said, we do not offer advice on how long the whole longitudinal series should be, neither could we comment on how dense the sampling points should be. This is still a pioneering project and you just have to try to know. As of our experience with the Lanier Project, we sampled the freshwater lake at roughly monthly basis for three years, each sample had about 3-5 Gbp Illumina paired-end (GA-II and HiSeq) reads, and were able to track HGT between populations. 

In terms of selecting sequencing technologies, currently metaHGT requires paired-end reads, preferably Illumina; but in theory, other sequencing platforms such as Ion Torrent and AB SOLid would work as well. As other sequencing technologies that would produce significantly longer sequences (e.g., Oxford Nanopore's platforms) are getting huge momentum and will likely to be proven to be a breakthrough in DNA sequencing, metaHGT will evolve to use also longer single reads to better resolve HGTs.

[How to install]

Install is simple, first you clone the repository to your local machine by:

$ github clone git://github.com/luo-chengwei/metaHGT

then you should "cd" to the local clone and run:

$ python setup.py install

You should be all set.

There are some libraries that are needed for running metaHGT:



