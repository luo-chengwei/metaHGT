#! /usr/bin/env python

__author__ = 'Chengwei Luo (luo.chengwei@gatech.edu)'
__version__ = '0.0.1'
__date__ = 'November 2013'


"""
metaHGT (meta-community Horizontal Gene Transfer tracker): 
	in-situ and real time HGT tracker for series metagenomes

Copyright(c) 2013 Chengwei Luo (luo.chengwei@gatech.edu)

	This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>

https://github.com/luo-chengwei/metaHGT

for help, type:
python metaHGT.py --help
"""

USAGE = \
"""Usage: %prog <required_parameters> [options]

metaHGT: in-situ and real time HGT tracker for series metagenomes

metaHGT is a platform for in-situ and real time HGT tracker series metagenomes. 
It utilizes the power of PE mapping and infers HGTs from mapping contrasts across
different time points.
It is written in Python, therefore it should run on Mac OS, and Unix/Linux. 

Add --help to see a full list of required and optional
arguments to run metaHGT.

Additional information can also be found at:
https://github.com/luo-chengwei/metaHGT/wiki

If you use metaHGT in your work, please cite it as:
<metaHGT citation here>

Copyright: Chengwei Luo, Konstantinidis Lab, Georgia Institute of Technology, 2013
"""

import sys
import os
import re
import glob
from optparse import OptionParser, OptionGroup
from subprocess import PIPE, Popen
import cPickle

from time import ctime, time
from datetime import timedelta

import calHGT
import bams

####################### PROJECT INFORMATION CLASS #############################
class ProjectInfo:
	def __init__(self):
		self.samples = []
		self.timepairs = []
		self.num_proc = 1
		self.quiet = False
		self.bam_dir = None
		self.assembly_dir = None
		self.reads_dir = None
		self.interleaved = []
		self.outdir = None
	
	def getBAMFiles(self, sample1, sample2):
		BAMs = []
		meshedSamplePairs = [(sample1, sample1), (sample2, sample2), 
							(sample1, sample2), (sample2, sample1)]
		for sampleA, sampleB in meshedSamplePairs:
			files = glob.glob(self.bam_dir + '/' + sampleA + '.vs.' + sampleB + '*bam')
			if len(files) == 0:
				sys.stderr.write('FATAL: Eror in fetching the BAM file for samples: %s and %s\n' % (sampleA, sampleB))
				exit(0)
			if len(files) > 1:
				sys.stderr.write('FATAL: Ambiguous naming for BAM files for samples: %s and %s\n' % (sampleA, sampleB))
				sys.stderr.write('       The following files are found:\n')
				for file in files:
					sys.stderr.write('       %s\n' % file)
				exit(1)
			else:
				BAMs.append(os.path.realpath(files[0]))
				
		return BAMs
		
	def getReadsFile(self, sample):
		files1 = glob.glob(self.reads_dir + '/' + sample + '*1.fa')
		if len(files1) == 0:
			files1 = glob.glob(self.reads_dir + '/' + sample + '*1.fastq')
		if len(files1) == 0:
			files1 = glob.glob(self.reads_dir + '/' + sample + '*1.fasta')
			
		files2 = glob.glob(self.reads_dir + '/' + sample + '*2.fa')
		if len(files2) == 0:
			files2 = glob.glob(self.reads_dir + '/' + sample + '*2.fastq')
		if len(files2) == 0:
			files2 = glob.glob(self.reads_dir + '/' + sample + '*2.fasta')
			
		if len(files1) == 0 and len(files2) == 0:
			files = glob.glob(self.reads_dir + '/' + sample + '*fa')
			if len(files) == 0:
				files = glob.glob(self.reads_dir + '/' + sample + '*fastq')
			if len(files) == 0:
				files = glob.glob(self.reads_dir + '/' + sample + '*fasta')
			
			if len(files) == 0:
				sys.stderr.write('FATAL: Eror in fetching the reads file for sample: %s\n' % sample)
				exit(0)
			
			if len(files) > 1:
				sys.stderr.write('FATAL: Ambiguous naming for reads file for sample: %\n' % sample)
				sys.stderr.write('       The following files are found:\n')
				for file in files:
					sys.stderr.write('       %s\n' % file)
				exit(0)
			else:
				return files
			
		else:	
			if len(files1) == 0:
				sys.stderr.write('FATAL: Eror in fetching the 5\' reads file for sample: %s\n' % sample)
				exit(0)
			if len(files) > 1:
				sys.stderr.write('FATAL: Ambiguous naming for reads file for sample: %\n' % sample)
				sys.stderr.write('       The following files are found:\n')
				for file in files1:
					sys.stderr.write('       %s\n' % file)
				exit(0)
			
			if len(files2) == 0:
				sys.stderr.write('FATAL: Eror in fetching the 3\' reads file for sample: %s\n' % sample)
				exit(0)
			if len(files2) > 1:
				sys.stderr.write('FATAL: Ambiguous naming for reads file for sample: %\n' % sample)
				sys.stderr.write('       The following files are found:\n')
				for file in files2:
					sys.stderr.write('       %s\n' % file)
				exit(0)
			
			return files1 + files2
	
	
	def getAssemblyFile(self, sample):
		files = glob.glob(self.assembly_dir + '/' + sample + '*fa')
		if len(files) == 0:
			sys.stderr.write('FATAL: Eror in fetching the assembly file for sample: %s\n' % sample)
			exit(0)
		if len(files) > 1:
			sys.stderr.write('FATAL: Ambiguous naming for assembly file for sample: %\n' % sample)
			sys.stderr.write('       The following files are found:\n')
			for file in files:
				sys.stderr.write('       %s\n' % file)
			exit(1)
		else:
			return os.path.realpath(files[0])
	
		
	def initProject(self, options):
		if os.path.exists(options.sample_list):
			for timepoint in open(options.sample_list, 'r'):
				self.samples.append(tuple(timepoint.rstrip('\n').split(';')))
		elif options.sample_list.count(':') > 0:
			for timepoint in options.sample_list.split(','):
				self.samples.append(tuple(timepoint.split(':')))
		else:
			sys.stderr.write('FATAL: Error in extracting samples, please check your input.\n')
		
		# init outdir
		self.outdir = options.out_dir
		if not os.path.exists(self.outdir):
			os.mkdir(self.outdir)
		
		# generate timepairs
		for timepoint1, timepoint2 in zip(self.samples[:-1], self.samples[1:]):
			for sample1 in timepoint1:
				for sample2 in timepoint2:
					self.timepairs.append((sample1, sample2))
		
		self.num_proc = options.num_proc
		self.quiet = options.quiet
		
		# you supply either BAM dir or assembly dir + reads dir
		# if you supply BAMs, then it skips the BWA+Samtoools step to generate the BAM files,
		# it would rely on you for generate the correct BAMs (sorted and indexed).
		# otherwise, this will try to generate all the BAMs files needs.
		
		if options.bam_dir and os.path.exists(options.bam_dir):
				self.bam_dir = options.bam_dir
				for sample1, sample2 in self.timepairs:
					bamfiles = self.getBAMFiles(sample1, sample2)
				
		else:
			if options.assembly_dir == None or options.reads_dir == None:
				sys.stderr.write('FATAL: You need to either supply the BAM file directory or the assembly and reads directory.\n')
				exit(0)
				
			if not os.path.exists(options.assembly_dir):
				sys.stderr.write('FATAL: cannot locate the assembly fasta files directory, you supplied: %s\n' % options.assembly_dir)
				eixt(0)	
			else:
				self.assembly_dir = options.assembly_dir
			
			if not os.path.exists(options.reads_dir):
				sys.stderr.write('FATAL: cannot locate the reads directory, you supplied: %s\n' % options.reads_dir)
				eixt(0)	
			else:
				self.reads_dir = options.reads_dir
			
			# test files	
			for timepoint in self.samples:
				for sample in timepoint:
					readsfile = self.getReadsFile(sample)
					assemblyfile = self.getAssemblyFile(sample)
					
			# test samtools and bwa
			bwaTest = Popen(options.bwa, shell=True, stdout=PIPE).stdout.read()
			if bwaTest == None or len(bwaTest) == 0:
				sys.stderr.write("FATAL: BWA not found in path!")
				exit(0)
				
			samtoolsTest = Popen(options.samtools, shell=True, stdout=PIPE).stdout.read()
			if samtoolsTest == None or len(samtoolsTest) == 0:
				sys.stderr.write("FATAL: samtools not found in path!")
				exit(0)
				
	# End of initProject

################################### MAIN ######################################
def main(argv = sys.argv[1:]):

	parser = OptionParser(usage = USAGE, version="Version: " + __version__)
	
	# Required arguments
	requiredOptions = OptionGroup(parser, "Required options",
								"These options are required to run BinGeR, and may be supplied in any order.")
	
	requiredOptions.add_option("-l", "--sample_list", type = "string", metavar = "FILE/STRING",
							help = "Text file containing all sample names, one per line, in longitudinal order; \
							replicates of the same time point should be separated by colon.\n\
							Alternatively, you can directly supply the sample names in longitudinal order, \
							timepoints separated by comma and replicates separated by colon.")

	requiredOptions.add_option("-o", "--out_dir", type = "string", metavar = "DIR",
							help = "Directory where intermediary and output files be reported to.")

	parser.add_option_group(requiredOptions)

	# Optional arguments that need to be supplied if not the same as default
	
	optOptions = OptionGroup(parser, "Optional parameters",
						"There options are optional, and may be supplied in any order.")
	
	optOptions.add_option("-b", "--bam_dir", type = "string", default = "BAMs", metavar = "DIR",
							help = "Directory where sorted bam files (reads versus assembly, same sample) are, \
							the naming should follow \"sample1.vs.sample2.*.bam\" convention. \
							NOTE: if you specify this option, metaHGT will assume that you have performed the BWA")
	
	optOptions.add_option("-a", "--assembly_dir", type = "string", default = "Assemblies", metavar = "DIR",
							help = "Directory where assembly files ni fasta format are, \
							the naming should follow \"sample.*.fa\" convention.\n\
							The tags of contigs should follow: binID.contigXX.* fashion, where binID is the \
							identifier of bins, and contigXX is the identifier of the contigs belong to the bin. \
							Unclassified (unbinned) contigs should be renamed as ID.contigXX.*, where ID can be \
							anything but the claimed binIDs.")
							
	optOptions.add_option("-r", "--reads_dir", type = "string", default = "Reads", metavar = "DIR",
							help = "Directory where reads are. They should be in fastq format, \
							and can be in both interleaved, or separate two files. \
							The naming should follow \"sample.*.fastq\" (interleaved) or \"sample.*.1.fastq\" \
							and \"sample.*.2.fastq\" (separate) convention.")
							
	optOptions.add_option("--bwa", type = "string", default = "bwa", metavar = "STRING",
							help = "Location of BWA (Li et al, Bioinformatics, 2009). Only needed if you haven't \
							generated the BAM files yourself; otherwise, please specify BAM_dir (-d/--BAM_dir).")
							
	optOptions.add_option("--samtools", type = "string", default = "samtools", metavar = "STRING",
							help = "Location of the Samtools binary (Li et al, Bioinformatics, 2009).Only needed \
							if you haven't generated the BAM files yourself; otherwise, please specify BAM dir (-d/--bam_dir).")

	optOptions.add_option("-t", "--num_proc", type = "int", default = 1, metavar = 'INT',
							help = "Number of processor for BinGeR to use [default: 1].")

	parser.add_option_group(optOptions)
	
	# runtime settings that could affect the file saving and message printing
	runtimeSettings = OptionGroup(parser, "Runtime settings",
						"There options are optional, and may be supplied in any order.")
						
	runtimeSettings.add_option("-q", "--quiet", default = False, action = "store_true",
								help = "Suppress printing detailed runtime information, only important messages will show [default: False].")

	parser.add_option_group(runtimeSettings)

	(options, args) = parser.parse_args(argv)
	
	if options.sample_list is None:
		parser.error("A list of samples is required!")
		exit(0)
		
	if options.out_dir is None:
		parser.error("An output directory is required to supply!")
		exit(0)
	
	# kickstart
	sys.stdout.write("metaHGT started at %s\n"%(ctime()))
	sys.stdout.flush()
	
	# if necessary, run bwa + samtools to generate sorted and indexed BAM files
	
	# check sanity of the files in required directories
	projInfo = ProjectInfo()
	projInfo.initProject(options)
	
	
	# end
	sys.stdout.write("metaHGT finished at %s\n"%(ctime()))
	sys.stdout.flush()


if __name__ == '__main__':
	main()
	
