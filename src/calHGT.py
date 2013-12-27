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
import math
from scipy import stats
import numpy as np
from scipy.stats import norm
from itertools import product
from operator import itemgetter, attrgetter
		
class Span:
	def __init__(self):
		self.reference = None
		self.left = 0
		self.right =  0

class BreakPoint:
	def __init__(self):
		self.reference = None
		self.direction = '>'
		self.loc = 0

class HGT:
	def __init__(self):
		self.genome_a = None
		self.genome_b = None
		self.reference_a = None
		self.reference_b = None
		self.breakpoint_a = BreakPoint()
		self.breakpoint_b = BreakPoint()
		self.percentage_a = 0.
		self.percentage_b = 0.
		self.pvalue = 1.
		self.adj_pvalue = 1.
		
	def printHGT(self):
		sys.stdout.write('%s\t%s\t%i\t%s\t%s\t%s\t%i\t%s\t' % (self.genome_a, self.reference_a, \
							self.breakpoint_a.loc, self.breakpoint_a.direction, self.genome_b, \
							self.reference_b, self.breakpoint_b.loc, self.breakpoint_b.direction))
		
		sys.stdout.write('%.4f\t%.4f\t%.4f\t%.4f\n' % (self.percentage_a, self.percentage_b, \
														self.pvalue, self.adj_pvalue))
														
	def strHGT(self):
		return '%s\t%s\t%i\t%s\t%s\t%s\t%i\t%s\t%.4f\t%.4f\t%.4f\t%.4f' % (self.genome_a, self.reference_a, \
							self.breakpoint_a.loc, self.breakpoint_a.direction, self.genome_b, \
							self.reference_b, self.breakpoint_b.loc, self.breakpoint_b.direction, \
							self.percentage_a, self.percentage_b, self.pvalue, self.adj_pvalue)

def BH_correction(HGTs):
	m = len(HGTs)
	sorted_HGTs	= sorted(HGTs, key = lambda x: x.pvalue)
	
	for k, x in enumerate(sorted_HGTs):
		sorted_HGTs[k].adj_pvalue = min(1., x.pvalue * (float(m)/(k+1)))
	
	return sorted_HGTs

# End of BH_correction

def GtoP(g,df):
	return 1-stats.chi2.cdf(g, df)
# End of GtoP
    
def flnf(f):
	return f*math.log(f) if f >0.5 else 0
# End of flnf

def gtest(table):
	a, b, c, d = table
	row1 = a+b
	row2 = c+d
	col1 = a+c
	col2 = b+d
     
	total = flnf(a+b+c+d)
	celltotals = flnf(a)+flnf(b)+flnf(c)+flnf(d)
	rowtotals = flnf(row1)+flnf(row2)
	coltotals = flnf(col1)+flnf(col2)
     
	gscore = 2*(celltotals + total-(rowtotals+coltotals))
	
	return GtoP(gscore, 1)   # return gscore

# End of gtest

def estimate_library(bamfh, ref_length, projInfo):
	# library size estimate
	library_sizes = []
	size_limit = 1e5
	for reference in ref_length:
		length = ref_length[reference]
		if length < projInfo.contig_length:
			continue
		try:
			reads = bamfh.fetch(reference)
		except ValueError:
			continue
			
		for read in reads:
			if read.mapq < projInfo.mapq:
				continue
			if float(read.alen)/read.rlen < projInfo.align_perc:
				continue
			if read.tid == read.mrnm and (read.is_reverse ^ read.mate_is_reverse):
				library_sizes.append(abs(read.pos - read.mpos) + 1+ read.rlen)
		if len(library_sizes) > size_limit:
			break
	
	mean = np.mean(np.array(library_sizes))
	std = np.std(np.array(library_sizes))
	
	return mean, std
	
# End of estimate_library

def infer_breakpoint(reference, reads10):
	bk = BreakPoint()
	bk.reference = reference
	
	positions = []
	for read in reads10:
		if read.is_secondary or read.mapq < 30:
			continue
		if read.tid == read.mrnm:
			continue
		if read.is_read1 and not read.is_reverse:
			positions.append((read.pos + read.qlen, '>'))
		elif read.is_read1 and read.is_reverse:
			positions.append((read.pos, '<'))
		elif read.is_read2 and not read.is_reverse:
			positions.append((read.pos + read.qlen, '>'))
		elif read.is_read2 and read.is_reverse:
			positions.append((read.pos, '<'))
		
		
	dir = map(itemgetter(1), positions)
	pos = map(itemgetter(0), positions)
	
	if dir.count('>') > dir.count('<') and dir.count('<') <= 1:
		bk.direction = '>'
		bk.loc = max(pos)
	elif dir.count('>') < dir.count('<') and dir.count('>') <= 1:
		bk.direction = '<'
		bk.loc = min(pos)
	else:
		return None
	
	return bk
	
# End of infer_breakpoint

def test_HGT(mappings0, breakpoint_a, mappings1, breakpoint_b):
	# init basic parameters
	ref0_self = [0, 0, 0, 0]
	ref1_self = [0, 0, 0, 0]
	crossing0 = [0, 0, 0, 0]
	crossing1 = [0, 0, 0, 0]
	
	# render ref0
	for index, mapping in enumerate(mappings0):
		for read in mapping:
			if read.mapq < 30 or read.is_secondary:
				continue
			if read.tid == read.mrnm:
				if read.pos < read.pnext:
					start, end = read.pos, read.pnext
				else:
					start, end = read.pnext, read.pos
				if start <= breakpoint_a.loc and end >= breakpoint_a.loc:
					ref0_self[index] += 1
				else:
					continue
			else:
				if breakpoint_a.direction == '>':
					if read.pos <= breakpoint_a.loc and not read.is_reverse:
						crossing0[index] += 1
				elif breakpoint_a.direction == '<':
					if read.pos >= breakpoint_a.loc and read.is_reverse:
						crossing0[index] += 1
				
	# render ref1
	for index, mapping in enumerate(mappings1):
		for read in mapping:
			if read.mapq < 30 or read.is_secondary:
				continue
			if read.tid == read.mrnm:
				if read.pos < read.pnext:
					start, end = read.pos, read.pnext
				else:
					start, end = read.pnext, read.pos
				if start <= breakpoint_b.loc and end >= breakpoint_b.loc:
					ref1_self[index] += 1
				else:
					continue
			else:
				if breakpoint_b.direction == '>':
					if read.pos <= breakpoint_b.loc and not read.is_reverse:
						crossing1[index] += 1
				elif breakpoint_b.direction == '<':
					if read.pos >= breakpoint_b.loc and read.is_reverse:
						crossing1[index] += 1
	
	# perform hypothesis testing
	if ref0_self[-1] + crossing0[-1] > 0:
		percentage_a = (100. * crossing0[-1]) / (ref0_self[-1] + crossing0[-1])
	else:
		percentage_a = None
	if ref1_self[-1] + crossing1[-1] > 0:
		percentage_b = (100. * crossing1[-1]) / (ref1_self[-1] + crossing1[-1])
	else:
		percentage_b = None
	
	# first need to satisfy that self-mapping has no crossing
	if  sum([crossing0[0], crossing1[0]]) > 0:
		return None, None, None
		
	table = [ref0_self[0] + ref1_self[0], crossing0[0] + crossing1[0], 
			ref0_self[-1] + ref1_self[-1], crossing0[-1] + crossing1[-1]]
	
	pvalue = gtest(table)
	
	return pvalue, percentage_a, percentage_b
	
# End of test_HGT
	
def infer_HGT(spans0, spans1, bamfhs):
	# bamfhs order: 00, 11, 01, 10
	HGTs = []
	for span0, span1 in product(spans0, spans1):
		try:
			ref0_reads10 = bamfhs[3].fetch(span0.reference, span0.left, span0.right)
			ref1_reads10 = bamfhs[3].fetch(span1.reference, span1.left, span1.right)
		except ValueError:
			continue
			
		breakpoint0 = infer_breakpoint(span0.reference, ref0_reads10)
		breakpoint1 = infer_breakpoint(span1.reference, ref1_reads10)
		
		if breakpoint0 == None or breakpoint1 == None:
			continue
		
		e = HGT()
		e.genome_a = span0.genome
		e.genome_b = span1.genome
		e.reference_a = span0.reference
		e.reference_b = span1.reference
		e.breakpoint_a = breakpoint0
		e.breakpoint_b = breakpoint1
		
		try:
			ref0_reads00 = bamfhs[0].fetch(span0.reference, span0.left, span0.right)
			ref0_reads11 = bamfhs[1].fetch(span0.reference, span0.left, span0.right)
			ref0_reads01 = bamfhs[2].fetch(span0.reference, span0.left, span0.right)
			ref0_reads10 = bamfhs[3].fetch(span0.reference, span0.left, span0.right)
			ref0_mappings = [ref0_reads00, ref0_reads11, ref0_reads01, ref0_reads10]
		
			ref1_reads00 = bamfhs[0].fetch(span1.reference, span1.left, span1.right)
			ref1_reads11 = bamfhs[1].fetch(span1.reference, span1.left, span1.right)
			ref1_reads01 = bamfhs[2].fetch(span1.reference, span1.left, span1.right)
			ref1_reads10 = bamfhs[3].fetch(span1.reference, span1.left, span1.right)
			ref1_mappings = [ref1_reads00, ref1_reads11, ref1_reads01, ref1_reads10]
		except ValueError:
			continue
			
		pvalue, percentage_a, percentage_b = test_HGT(ref0_mappings, e.breakpoint_a, \
													ref1_mappings, e.breakpoint_b)
		
		if percentage_a == None or percentage_b == None or pvalue == None:
			continue
		else:					
			e.percentage_a = percentage_a
			e.percentage_b = percentage_b
			e.pvalue = pvalue
		
		HGTs.append(e)
		
	return HGTs
		
# end of infer_HGT	

def HGT_search(lib_stats, bamfhs, projInfo):
	sampleIDs = {}
	for sample in projInfo.samples:
		for x in sample: sampleIDs[x] = 1
	
	
	## bamfhs: [bam00, bam11, bam01, bam10]
	references0 = bamfhs[0].references
	lengths0 = bamfhs[0].lengths
	refs0 = {}
	for ref, lgth in zip(references0, lengths0):
		refs0[ref] = lgth
			
	references1 = bamfhs[1].references
	lengths1 = bamfhs[1].lengths
	refs1 = {}
	for ref, lgth in zip(references1, lengths1):
		refs1[ref] = lgth
		
	mean0, std0 = lib_stats[0]
	mean1, std1 = lib_stats[1]
	
	G = {}  # Graph that holds the hotspot candidates
	
	# scan for hotspot candidate
	if not projInfo.quiet:
		sys.stdout.write('Scanning for hotspots...\n')

	for read in bamfhs[3].fetch(until_eof = True):
		if read.mapq < projInfo.mapq:
			continue
		if float(read.alen)/read.rlen < projInfo.align_perc:
			continue
		if read.is_paired and not read.is_secondary:
			if read.mrnm == read.tid:
				continue
		if read.is_reverse:
			strand0 = '-'
		else:
			strand0 = '+'
		if read.mate_is_reverse:
			strand1 = '-'
		else:
			strand1 = '+'
				
		ref1 = bamfhs[3].getrname(read.mrnm)
		ref0 = bamfhs[3].getrname(read.tid)
		
		if ref0 not in refs0 or ref1 not in refs0:
			continue
			
		# length filter
		if refs0[ref0] < projInfo.contig_length or refs0[ref1] < projInfo.contig_length:
			continue
		
		(ref0, pos0, strand0), (ref1, pos1, strand1) = \
			sorted([(ref0, read.pos, strand0), (ref1, read.mpos, strand1)], key = lambda x: x[0])
		binID0 = re.search('(.+)\.contig\d+', ref0).group(1)
		binID1 = re.search('(.+)\.contig\d+', ref1).group(1)
		
		# cross bin filter
		if binID0 in sampleIDs or binID1 in sampleIDs or binID1 == binID0:
			continue
		
		if (binID0, binID1) not in G:
			G[(binID0, binID1)] = {}
		if (ref0, ref1) not in G[(binID0, binID1)]:
			G[(binID0, binID1)][(ref0, ref1)] = []

		G[(binID0, binID1)][(ref0, ref1)].append((pos0, strand0, pos1, strand1))
		
	# scan through all candidate hotspot and calculate new genotype 
	# and corresponding percentage and p-value	
	if not projInfo.quiet:
		sys.stdout.write('Calculating HGT precise locations and likelihood...\n')
	
	HGTs = []
	for (binID0, binID1) in G:
		for (ref0, ref1) in G[(binID0, binID1)]:
			lgth0 = refs0[ref0]
			lgth1 = refs0[ref1]
			refined_positions = []
			for index, (pos0, strand0, pos1, strand1) in enumerate(G[(binID0, binID1)][(ref0, ref1)]):
				if pos0 <= mean1 or pos0 >= lgth0 - mean1 or pos1 <= mean1 or pos1 >= lgth1 - mean1:
					continue
				refined_positions.append((pos0, strand0, pos1, strand1))
			if len(refined_positions) < projInfo.min_links:
				continue
			
			span_links = []
			sorted_refined_positions_0 = sorted(refined_positions, key = lambda x: x[0])
			sorted_refined_positions_1 = sorted(refined_positions, key = lambda x: x[2])
			
			## get the spans in ref0
			old_i = 0
			total_ele = 0
			link_groups = []
			spans0 = []
			for i in range(len(sorted_refined_positions_0) - 1):
				gap = sorted_refined_positions_0[i+1][0] - sorted_refined_positions_0[i][0]
				if gap > mean1 + std1:
					link_groups.append(sorted_refined_positions_0[old_i:i+1])
					total_ele += len(sorted_refined_positions_0[old_i:i+1])
					old_i = i+1
			if total_ele < len(sorted_refined_positions_0):
				link_groups.append(sorted_refined_positions_0[old_i:])
			
			for link_group in link_groups:
				s = Span()
				s.genome = binID0
				s.reference = ref0
				s.left = max(int(link_group[0][0] - mean1), 0)
				s.right = min(int(link_group[-1][0] + mean1), lgth0)
				spans0.append(s)
			
			## get the spans in ref1
			old_i = 0
			total_ele = 0
			link_groups = []
			spans1 = []
			for i in range(len(sorted_refined_positions_1) - 1):
				gap = sorted_refined_positions_1[i+1][2] - sorted_refined_positions_1[i][2]
				if gap > mean1 + std1:
					link_groups.append(sorted_refined_positions_1[old_i:i+1])
					total_ele += len(sorted_refined_positions_1[old_i:i+1])
					old_i = i+1
			if total_ele < len(sorted_refined_positions_1):
				link_groups.append(sorted_refined_positions_1[old_i:])
			
			for link_group in link_groups:
				s = Span()
				s.genome = binID1
				s.reference = ref1
				s.left = max(int(link_group[0][2] - mean1), 0)
				s.right = min(int(link_group[-1][2] + mean1), lgth1)
				spans1.append(s)
				
			## examine between any combo of two spans
			results = infer_HGT(spans0, spans1, bamfhs)
			HGTs += results
		
	if not projInfo.quiet:
		sys.stdout.write('Done HGT inferences.\n')
				
	return HGTs
	
# End of HGT_search

def HGTscan(bamfhs, lib_stats, projInfo):
	
	# select hotspots and calculate HGT likelihood
	HGTs = HGT_search(lib_stats, bamfhs, projInfo)
	
	# FDR adjustment
	if not projInfo.quiet:
		sys.stdout.write('Correcting for FDR using Benjamini-Hochberg methods...\n')
	adj_HGTs = BH_correction(HGTs)
	
	return adj_HGTs
	
# End of HGTscan

def calHGT(projInfo):
	if not projInfo.quiet:
		sys.stdout.write('Now calculating HGT events...\n')
	
	all_HGTs = {}
	
	for timepair in projInfo.timepairs:
		if not projInfo.quiet:
			sys.stdout.write('[%s.vs.%s] Initiating...\n' % (timepair[0], timepair[1]))
		
		# get all the BAMs
		if not projInfo.quiet:
			sys.stdout.write('[%s.vs.%s] Interpreting BAMs...\n' % (timepair[0], timepair[1]))
		bamfile00, bamfile11, bamfile01, bamfile10 = projInfo.getBAMFiles(timepair[0], timepair[1])
		
		# file handles
		bam00 = pysam.Samfile(bamfile00, 'rb')
		bam11 = pysam.Samfile(bamfile11, 'rb')
		bam01 = pysam.Samfile(bamfile01, 'rb')
		bam10 = pysam.Samfile(bamfile10, 'rb')
		bamfhs = [bam00, bam11, bam01, bam10]
		
		# estimate library size and std
		if not projInfo.quiet:
			sys.stdout.write('Now estimating the read library size...\n')
	
		refs0 = {}
		for ref, length in zip(bam00.references, bam00.lengths):
			refs0[ref] = length
		mean0, std0 = estimate_library(bam00, refs0, projInfo)
		if not projInfo.quiet:
			sys.stdout.write('[Library 1] mean +/- std: %.3f +/- %.3f\n' % (mean0, std0))
	
		refs1 = {}
		for ref, length in zip(bam11.references, bam11.lengths):
			refs1[ref] = length
		mean1, std1 = estimate_library(bam11, refs1, projInfo)
		if not projInfo.quiet:
			sys.stdout.write('[Library 2] mean +/- std: %.3f +/- %.3f\n' % (mean1, std1))
		
		lib_stats = [(mean0, std0), (mean1, std1)]
		
		# run HGT calculation
		HGTs = HGTscan(bamfhs, lib_stats, projInfo)

		for bfh in bamfhs:
			bfh.close()
	
		all_HGTs[timepair] = HGTs

	if not projInfo.quiet:
		sys.stdout.write('All HGT predication finished.\n')

	return all_HGTs				
	
			