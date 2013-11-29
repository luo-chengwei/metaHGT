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

import sys
import os
import re
import pysam

class location:
	def __init__(self):
		self.contig = None
		self.start

class HGT:
	def __init__(self):
		self.genomeA = None
		self.genomeB = None
		self.locationA = location()
		self.locationB = location()
		self.pvalue = 1
		self.percentage = 0.
		

def calHGT(projInfo):
	if not projInfo.quiet:
		sys.stdout.write('Now calculating HGT events...\n')
	
	for timepair in projInfo.timepairs:
		if not projInfo.quiet:
			sys.stdout.write('[%s.vs.%s] Initiating...\n' % (timepair[0], timepair[1]))
		
		HGTs = []
		# get all the BAMs
		if not projInfo.quiet:
			sys.stdout.write('[%s.vs.%s] Interpreting BAMs...\n' % (timepair[0], timepair[1]))
		bamfile00, bamfile11, bamfile01, bamfile10 = projInfo.getBAMFiles(timepair[0], timepair[1])
		
		bam00 = pysam.Samfile(bamfile00, 'rb')
		bam11 = pysam.Samfile(bamfile11, 'rb')
		bam01 = pysam.Samfile(bamfile01, 'rb')
		bam10 = pysam.Samfile(bamfile10, 'rb')
		
		references0 = bam00.references
		lengths0 = bam00.lengths
		references1 = bam11.references
		lengths1 = bam11.lengths
		
		# find hotspots
		for reference, length in zip(references0, lengths0):
			reads = bam00.fetch(reference)
			for read in reads:
				if read.is_paired:
					print read.name
				
			break
			
		# estimate library size
		
		bam00.close()
		bam11.close()
		bam10.close()
		bam01.close()
		
	# finished
	if not projInfo.quiet:
			sys.stdout.write('HGT events calculated.\n')
			
			