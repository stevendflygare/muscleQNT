#!/opt/local/bin/python2.7

import sys
import re
from PIL import Image
import scipy.ndimage as si
import numpy as np
import random
import scipy.stats as ss
import argparse
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import itertools

#command line parsing 
parser = argparse.ArgumentParser()
parser.add_argument("fiber_file", help="fiber size input file")
parser.add_argument("image_metrics", help="image metrics input file")
parser.add_argument("-c","--conversion",help="pixel to micro-meter conversion (multiplicative factor)",type=float)
parser.add_argument("-rf","--relative",help="use relative frequency instead of raw counts in histogram",action="store_true")
parser.add_argument("-bs","--bin_sizes",help="comma separated list of pixel break values for bins",type=str)
parser.add_argument("-hb",help="number of histogram bins",type=int,default=20)
parser.add_argument("-xl",help="x-axis label",type=str,default="")
parser.add_argument("-yl",help="y-axis label",type=str,default="")
args = parser.parse_args()

def load_fiber_sizes(fiber_file,mice_hash): 
	f = open(fiber_file,'r')
	for line in f:
		data = re.split(",",line.strip())
		m = re.search("([0-9a-zA-Z]+)_([0-9a-zA-Z]+)_(\d+)\.",data[0])
		if m:
			if m.group(1) in mice_hash:
				mice_hash[m.group(1)][0].append(map(int,data[1:]))
			else:
				print "problem, mouse not found going to exit: " + m.group(1)
				sys.exit(0)
		else:
			print "problem, image name is not recognized, please check naming scheme of images: " + data[0]
			print "exiting..."
			sys.exit(0)
	f.close()

def compute_ks_test(mice_hash):
	image_names = []
	for k in mice_hash:
		image_names.append([])
		for im in mice_hash[k][1]:
			image_names[-1].append((k,im))
	for p in itertools.product(*image_names):	
		mutant = []
		wild_type = []
		#print p
		for i in p:
			if mice_hash[i[0]][2] == 0:
				wild_type.extend(mice_hash[i[0]][0][mice_hash[i[0]][1].index(i[1])])
			else:
				mutant.extend(mice_hash[i[0]][0][mice_hash[i[0]][1].index(i[1])])
		ks_stat, p_value =  ss.ks_2samp(wild_type, mutant)
		if p_value > .05:
			print ks_stat, p_value
			print p				


def show_histogram(wt_histograms, mt_histograms, bins):
	wt_hist = np.mean(wt_histograms,axis=0)
	wt_std = np.std(wt_histograms,axis=0)
	mt_hist = np.mean(mt_histograms,axis=0)
	mt_std = np.std(mt_histograms,axis=0)
			
	#print "making histogram"		
	#use relative frequencies instead of the density
	ind = np.arange(len(wt_hist))  # the x locations for the groups
	width = 0.35       # the width of the bars
	rects1 = plt.bar(ind, wt_hist, width, color='b', yerr=wt_std, error_kw=dict(elinewidth=4, ecolor='green'))
	rects2 = plt.bar(ind+width, mt_hist, width, color='r', yerr=mt_std, error_kw=dict(elinewidth=4, ecolor='pink'))
	a = list(wt_hist)
	a.extend(mt_hist)
	if not args.relative:
		plt.ylim([0, int(max(a)*1.10)])
	else:
		plt.ylim([0, max(a)+.1])
	tick_vals = []
	print "\nbins start and end pixel values: "
	for i in range(len(bins)):
		tick_vals.append(i+1)
		if i > 0:
			print "\tbin "+str(i)+": "+str(bins[i-1])+"-"+str(bins[i])				
	#for i in range(len(bins)):
	#	if i < len(bins)-1:
	#		if not conversion:
				#tick_vals.append(int(np.mean(bins[i:i+2])/100))
	#			tick_vals.append(int(np.mean(bins[i:i+2])))
	#		else:
	#			tick_vals.append(round(np.mean(bins[i:i+2])*conversion,2))
	#plt.xticks(ind,tick_vals)
	plt.xticks(ind+width,tick_vals)
	plt.xlabel(args.xl)
	
	plt.ylabel(args.yl)
	plt.legend( (rects1[0], rects2[0]), ('Wild Type', 'Mutant') )
	plt.show()
	#plt.savefig(hist_name)
	plt.clf()	


def compile_rf_histogram(mice_hash,hist_name,hb,user_bins):
	#do permutations to get all possibilities and then let the final histogram be the average
	image_names = []
	mt_histograms = []
	wt_histograms = []
	p_index = 0
	bins = []
	if user_bins:
		bins = map(int,re.split(",",user_bins))
		hb = len(bins)-1
	for k in mice_hash:
		image_names.append([])
		for im in mice_hash[k][1]:
			image_names[-1].append((k,im))
	for p in itertools.product(*image_names):	
		mutant = []
		wild_type = []
		#print p
		for i in p:
			if mice_hash[i[0]][2] == 0:
				wild_type.extend(mice_hash[i[0]][0][mice_hash[i[0]][1].index(i[1])])
			else:
				mutant.extend(mice_hash[i[0]][0][mice_hash[i[0]][1].index(i[1])])
		if p_index == 0:
			p_index+=1
			wt_hist, tmpbins = np.histogram(wild_type,hb)
			if not bins:
				bins = tmpbins[:]
			wt_relfreq = wt_hist/float(sum(wt_hist))
			wt_histograms.append(wt_relfreq)
			mt_hist, tmpbins = np.histogram(mutant,bins)
			mt_relfreq = mt_hist/float(sum(mt_hist))
			mt_histograms.append(mt_relfreq)
		else:
			wt_hist, tmpbins = np.histogram(wild_type,bins)
			wt_relfreq = wt_hist/float(sum(wt_hist))
			wt_histograms.append(wt_relfreq)
			mt_hist, tmpbins = np.histogram(mutant,bins)
			mt_relfreq = mt_hist/float(sum(mt_hist))
			mt_histograms.append(mt_relfreq)			
	show_histogram(wt_histograms, mt_histograms, bins)
	#print np.matrix(wt_histograms)

def compile_count_histogram(mice_hash,hist_name,hb,user_bins):
	#do permutations to get all possibilities and then let the final histogram be the average
	image_names = []
	mt_histograms = []
	wt_histograms = []
	p_index = 0
	bins = []
	if user_bins:
		bins = map(int,re.split(",",user_bins))
		hb = len(bins)-1
	for k in mice_hash:
		image_names.append([])
		for im in mice_hash[k][1]:
			image_names[-1].append((k,im))
	#need to equalize number of fibers between groups, so find smallest number among all the groups
	#and randomly select this number of fibers
	counts = []
	for p in itertools.product(*image_names):	
		mutant = []
		wild_type = []
		#print p
		for i in p:
			if mice_hash[i[0]][2] == 0:
				wild_type.extend(mice_hash[i[0]][0][mice_hash[i[0]][1].index(i[1])])
			else:
				mutant.extend(mice_hash[i[0]][0][mice_hash[i[0]][1].index(i[1])])
		counts.append(len(mutant))
		counts.append(len(wild_type))
	
	sample_size = min(counts)
	print "number of fibers to sample from each collection of mutant/wild type: " + str(sample_size)
	for p in itertools.product(*image_names):	
		mutant = []
		wild_type = []
		#print p
		for i in p:
			if mice_hash[i[0]][2] == 0:
				wild_type.extend(mice_hash[i[0]][0][mice_hash[i[0]][1].index(i[1])])
			else:
				mutant.extend(mice_hash[i[0]][0][mice_hash[i[0]][1].index(i[1])])	
		mutant = random.sample(mutant,sample_size)
		wild_type = random.sample(wild_type,sample_size)
		if p_index == 0:
			p_index+=1
			wt_hist, tmpbins = np.histogram(wild_type,hb)
			if not bins:
				bins = tmpbins[:]
			wt_histograms.append(wt_hist)
			mt_hist, tmpbins = np.histogram(mutant,bins)
			mt_histograms.append(mt_hist)
		else:
			wt_hist, tmpbins = np.histogram(wild_type,bins)
			wt_histograms.append(wt_hist)
			mt_hist, tmpbins = np.histogram(mutant,bins)
			mt_histograms.append(mt_hist)	
	show_histogram(wt_histograms, mt_histograms, bins)		


def create_decile_plot(mice_hash,plot_name):
	xvals = range(0,105,5)
	for m in mice_hash:
		if mice_hash[m][2] == 0:
			y = np.percentile(mice_hash[m][0],xvals)
			plt.plot(xvals,y,"bo-")
	for m in mice_hash:
		if mice_hash[m][2] == 1:
			y = np.percentile(mice_hash[m][0],xvals)
			plt.plot(xvals,y,"ro-")
	plt.axis([0,110,0,2000])				
	plt.savefig(plot_name)
	plt.clf()	

f = open(args.image_metrics,'r')
index = 0
mice = {} #key is mouse id, then the value of each key is [[fibers],[images],mutant/wildtype]
m_slow = {}
m_fast = {} 
for line in f:
	if index > 0:
		data = re.split(",",line.strip())
		m = re.search("([0-9a-zA-Z]+)_([0-9a-zA-Z]+)_(\d+)\.",data[0])
		if m:
			if m.group(1) in mice:
				mice[m.group(1)][1].append(m.group(2))
				mice[m.group(1)][2] = int(m.group(3))
				m_slow[m.group(1)][1].append(m.group(2))
				m_slow[m.group(1)][2] = int(m.group(3))
				m_fast[m.group(1)][1].append(m.group(2))
				m_fast[m.group(1)][2] = int(m.group(3))
			else:
				mice[m.group(1)] = [[],[],0]
				mice[m.group(1)][1].append(m.group(2))
				mice[m.group(1)][2] = int(m.group(3))
				m_slow[m.group(1)] = [[],[],0]
				m_slow[m.group(1)][1].append(m.group(2))
				m_slow[m.group(1)][2] = int(m.group(3))
				m_fast[m.group(1)] = [[],[],0]
				m_fast[m.group(1)][1].append(m.group(2))
				m_fast[m.group(1)][2] = int(m.group(3))
		else:
			print "problem, image name is not recognized, please check naming scheme of images: " + data[0]
			print "exiting..."
			sys.exit(0)
	index+=1	
f.close()

load_fiber_sizes(args.fiber_file,mice)
if args.relative:
	compile_rf_histogram(mice,"fiber_histogram.svg",args.hb,args.bin_sizes)
else:
	compile_count_histogram(mice,"fiber_histogram.svg",args.hb,args.bin_sizes)
	
#chi2, p, dof, ex = ss.chi2_contingency([wt_hist,mt_hist])
#print "chi2 test p-value: " + str(p)

#plotarray = []
#plotarray.append(mutant)
#plotarray.append(wild_type)
#plt.boxplot(plotarray)
#plt.xticks(np.arange(3),("","mutant","wild type"))
#plt.savefig("box_plots.png")

#plt.hist(plotarray,20,histtype='bar',color=['red','blue'],label=['Mutant','Wild Type'])
#plt.legend()
#plt.savefig("fiber_histogram.svg")



