#!/opt/local/bin/python2.7
import analyze_all
from analyze_all import *
import analyze_slow
from analyze_slow import *

from PIL import Image
import numpy as np
import scipy.ndimage as si
import re
import os
import sys
import time
import argparse

#command line parsing 
parser = argparse.ArgumentParser()
parser.add_argument("-bs","--box_size",help="indicate the box size to use",type=int,default=7,choices=range(3,25,2))
parser.add_argument("-off","--offset",help="indicate the offset to use (integer)",type=int,default=9)
parser.add_argument("-es","--erosion_steps",help="indicate the number of times to erode (integer)",type=int,default=2,choices=range(7))
parser.add_argument("-co","--colocate",help="indicates whether to colocate main channel with another",type=str,choices=["r","g","b"])
parser.add_argument("-ol","--outline",help="indicates whether to find the outline when computing between feature space",action="store_true")
parser.add_argument("-mxfb","--max_fiber_size",help="specifies maximum fiber size (in pixels)",type=int,default=2000)
parser.add_argument("-mnfb","--min_fiber_size",help="specifies minimum fiber size (in pixels)",type=int,default=10)
parser.add_argument("-s","--save",help="file name to save quantified and colored image",type=str)
parser.add_argument("image",type=str,help="image file")
parser.add_argument("chan",type=str,help="specify an image channel",choices=["r","g","b"])
args = parser.parse_args()

if args.offset < -50 or args.offset > 50:
	print "offset parameter too extreme, please use an integer between -50 and 50"
	sys.exit(0)

channel = 0 #default channel is red
if args.chan == "g":
	channel = 1
elif args.chan == "b":
	channel = 2

co_channel = -1
if args.colocate != None:
	if args.colocate == args.chan:
		print "can't colocate with the same channel, going to exit"
		sys.exit(0)
	co_channel = 0
	if args.colocate == "g":
		co_channel = 1
	elif args.colocate == "b":
		co_channel = 2

start = time.time()
	
#extract the specified channel and make 8 bit image
im_rgb = Image.open(args.image) 
if im_rgb.mode != "RGB":
	print "Image needs to be an RGB, this is not.  Going to skip image " + tmp_img
	sys.exit(0)
rgb_pix = im_rgb.getdata()
vis_rgb = list(rgb_pix)
bit_pix = []
for p in rgb_pix:
	bit_pix.append(p[channel])
im = Image.new("L",im_rgb.size)
im.putdata(bit_pix)

#quantify all fibers
vis_pix = np.zeros(im.size[0]*im.size[1],dtype=np.uint8)
all = analyze_all(vis_pix, args, im, vis_rgb)
slow = None

if co_channel >= 0: #quantify slow fibers
	co_pix = []
	for p in rgb_pix:
		co_pix.append(p[co_channel])
	im_co = Image.new("L",im_rgb.size)
	im_co.putdata(co_pix)
	#im_co.show()
	slow = analyze_slow(vis_pix, args, im_co, vis_rgb)	
	#need to get just fast fibers (look for 1 in vis_pix)
	fast_pix = vis_pix*(vis_pix == 1)
	labeled_array, num_features = si.label(np.array(fast_pix).reshape(im_co.size[1],im_co.size[0]),structure=[[1,1,1],[1,1,1],[1,1,1]])	
	label_hist = si.measurements.histogram(labeled_array,1,num_features,num_features)
	condition = (label_hist > args.min_fiber_size) & (label_hist < args.max_fiber_size)
	fiber_sizes = label_hist[condition]
	all.sizes = fiber_sizes[:] #overwrite the all fiber class sizes to hold just the fast fibers

#visualize
im_rgb.putdata(vis_rgb)
im_rgb.show()
if args.save:
	im_rgb.save(args.save+".png")

fiber_sizes = []
num_fast, num_slow = 0, 0
if slow:
	fiber_sizes.extend(slow.sizes)
	num_fast = len(all.sizes)
	num_slow = len(slow.sizes)
fiber_sizes.extend(all.sizes)
fiber_number = len(fiber_sizes)
mean_fiber = np.mean(fiber_sizes)
fiber_variance = np.var(fiber_sizes)
fiber_area = sum(fiber_sizes)
subtracted_area = 0
#subtracted_area = main_structure_area - fiber_area
print "number of fibers: " + str(fiber_number)
print "mean fiber size: " + str(mean_fiber)
print "fiber size variance: " + str(fiber_variance)
print "total fiber area: " + str(fiber_area)

print "time to process: " + str(time.time() - start)
	
	
	
	