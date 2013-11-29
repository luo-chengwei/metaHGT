#! /usr/bin/env python

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

import os
import sys
import multiprocessing as mp
from subprocess import PIPE, Popen, call

def bwaIndex(x):
	bwa, assembly, db, log = x
	cmd = [bwa, "index", "-p", db, assembly]
	logfh = open(log, 'a')
	retcode = call(cmd, stdout = logfh, stderr = logfh)
	logfh.close()

def mapping2bam(x):
	bwa, samtools, num_proc, reads, db, bam, log = x
	
	logfh = open(log, 'a')
	
	if len(reads) == 1:
		readfile = reads[0]
		cmd = " ".join([bwa, "mem", "-t", str(num_proc), db, readfile,
					"|", samtools, "view", "-bS", "-",
					"|", samtools, "sort", "-", bam])

	elif len(reads) == 2:
		readfile1, readfile2 = reads
		cmd = " ".join([bwa, "mem", "-t", str(num_proc), db, readfile1, readfile2,
					"|", samtools, "view", "-bS", "-",
					"|", samtools, "sort", "-", bam])
	
	sam2bam = call(cmd, shell = True, stdout = logfh, stderr = logfh)
	index_bam = call([samtools, "index", bam + '.bam'], stdout = logfh, stderr = logfh)
		
	logfh.close()

def genBAMs(projInfo):
	if not projInfo.quiet:
		sys.stdout.write('Now generating sorted BAM files...\n')
		
	projInfo.bam_dir = projInfo.outdir + '/BAMs/'
	if not os.path.exists(projInfo.bam_dir):
		os.mkdir(projInfo.bam_dir)
		
	# index assemblies first
	if not projInfo.quiet:
		sys.stdout.write('  Indexing assemblies...\n')
	indexCMDs = []
	for sampleTuple in projInfo.samples:
		for sample in sampleTuple:
			assemblyFile = projInfo.getAssemblyFile(sample)
			indexDB = projInfo.bam_dir + sample
			log = projInfo.bam_dir + sample + '.index.log'
			indexCMDs.append((projInfo.bwa, assemblyFile, indexDB, log))
	
	pool = mp.Pool(projInfo.num_proc)
	pool.map_async(bwaIndex, indexCMDs)
	pool.close()
	pool.join()
	
	
	if not projInfo.quiet:
		sys.stdout.write('  Indexing done.\n')
	
	# mapping using BWA and pipe through samtools
	if not projInfo.quiet:
		sys.stdout.write('  Cross mapping reads onto assemblies...\n')
	mappingCMDs = []
	for timepair in projInfo.timepairs:
			# B-vs-A mapping
			indexDB = projInfo.bam_dir + timepair[0]
			readfiles = projInfo.getReadsFile(timepair[1])
			bam = projInfo.bam_dir + timepair[1] + '.vs.' + timepair[0] + '.sorted'
			log = projInfo.bam_dir + timepair[1] + '.vs.' + timepair[0] + '.mapping2bam.log'
			cmd = (projInfo.bwa, projInfo.samtools, projInfo.num_proc, 
					readfiles, indexDB, bam, log)
			mappingCMDs.append(cmd)
			
			# A-vs-B mapping
			indexDB = projInfo.bam_dir + timepair[1]
			readfiles = projInfo.getReadsFile(timepair[0])
			bam = projInfo.bam_dir + timepair[0] + '.vs.' + timepair[1] + '.sorted'
			log = projInfo.bam_dir + timepair[0] + '.vs.' + timepair[1] + '.mapping2bam.log'
			cmd = (projInfo.bwa, projInfo.samtools, projInfo.num_proc, 
					readfiles, indexDB, bam, log)
			mappingCMDs.append(cmd)
			
			# A-vs-A mapping
			indexDB = projInfo.bam_dir + timepair[0]
			readfiles = projInfo.getReadsFile(timepair[0])
			bam = projInfo.bam_dir + timepair[0] + '.vs.' + timepair[0] + '.sorted'
			log = projInfo.bam_dir + timepair[0] + '.vs.' + timepair[0] + '.mapping2bam.log'
			cmd = (projInfo.bwa, projInfo.samtools, projInfo.num_proc, 
					readfiles, indexDB, bam, log)
			mappingCMDs.append(cmd)
			
			# B-vs-B mapping
			indexDB = projInfo.bam_dir + timepair[1]
			readfiles = projInfo.getReadsFile(timepair[1])
			bam = projInfo.bam_dir + timepair[1] + '.vs.' + timepair[1] + '.sorted'
			log = projInfo.bam_dir + timepair[1] + '.vs.' + timepair[1] + '.mapping2bam.log'
			cmd = (projInfo.bwa, projInfo.samtools, projInfo.num_proc, 
					readfiles, indexDB, bam, log)
			mappingCMDs.append(cmd)
			
	pool = mp.Pool(projInfo.num_proc)
	pool.map_async(mapping2bam, mappingCMDs)
	pool.close()
	pool.join()
	
	
	if not projInfo.quiet:
		sys.stdout.write('  Cross mapping done.\n')
	
	# finished
	if not projInfo.quiet:
		sys.stdout.write('BAM files done.\n')
		