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
import argparse
import time

#command line parsing 
parser = argparse.ArgumentParser()
parser.add_argument("-bs","--box_size",help="indicate the box size to use",type=int,default=7,choices=range(3,25,2))
parser.add_argument("-off","--offset",help="indicate the offset to use (integer)",type=int,default=9)
parser.add_argument("-es","--erosion_steps",help="indicate the number of times to erode (integer)",type=int,default=2,choices=range(7))
parser.add_argument("-co","--colocate",help="indicates whether to colocate main channel with another",type=str,choices=["r","g","b"])
parser.add_argument("-ol","--outline",help="indicates whether to find the outline when computing between feature space",action="store_true")
parser.add_argument("-mxfb","--max_fiber_size",help="specifies maximum fiber size (in pixels)",type=int,default=2000)
parser.add_argument("-mnfb","--min_fiber_size",help="specifies minimum fiber size (in pixels)",type=int,default=10)
parser.add_argument("dir",type=str,help="image file directory")
parser.add_argument("chan",type=str,help="specify an image channel",choices=["r","g","b"])
args = parser.parse_args()

if args.offset < -50 or args.offset > 50:
	print "offset parameter too extreme, please use an integer between -50 and 50"
	sys.exit(0)

images = []
path = (args.dir).rstrip("/") + "/"		
files = os.listdir(path)
for f in files:
	if not re.search("^\.",f):
		images.append(path+f)	

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

out = open("image_metrics.csv",'w')
col_names = ["image file", "number of fibers", "mean fiber size (pixels)", "fiber variance (pixels)", "total fiber area (pixels)", "subtracted area (pixels)", "number of fast fibers", "number of slow fibers" ]
out.write(",".join(col_names)+"\n")
out2 = open("fiber_sizes.csv",'w')
fast_out = open("fast_fibers.csv",'w')
slow_out = open("slow_fibers.csv",'w')	

start = time.time()
			
print "\n\n---------------------------------------------"
print "Running all images, please be patient, progress will be printed out as images are processed"

img_count = 1
for tmp_img in images:
	print "processing image " + str(img_count) + " of " + str(len(images))
	img_count+=1
	
	#extract the specified channel and make 8 bit image
	im_rgb = Image.open(tmp_img) 
	if im_rgb.mode != "RGB":
		print "Image needs to be an RGB, this is not.  Going to skip image " + tmp_img
		continue
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
	main_structure_area = 0
	if args.outline: #go through each column and add number of pixels between first and last non-zero terms
		tmp_vis_pix = vis_pix.reshape(im.size[1],im.size[0])
		for i in range(im.size[0]):
			min = 0
			max = 0
			for j in range(im.size[1]):
				if tmp_vis_pix[j,i] > 0:
					if min == 0:
						min = j
					if j > max:
						max = j
			main_structure_area += (max - min)											
	
	slow = None
	
	if co_channel >= 0: #quantify slow fibers
		co_pix = []
		for p in rgb_pix:
			co_pix.append(p[co_channel])
		im_co = Image.new("L",im_rgb.size)
		im_co.putdata(co_pix)
		slow = analyze_slow(vis_pix, args, im_co, vis_rgb)	
		#need to get just fast fibers (look for 1 in vis_pix)
		fast_pix = vis_pix*(vis_pix == 1)
		labeled_array, num_features = si.label(np.array(fast_pix).reshape(im_co.size[1],im_co.size[0]),structure=[[1,1,1],[1,1,1],[1,1,1]])	
		label_hist = si.measurements.histogram(labeled_array,1,num_features,num_features)
		condition = (label_hist > args.min_fiber_size) & (label_hist < args.max_fiber_size)
		fiber_sizes = label_hist[condition]
		all.sizes = fiber_sizes[:] #overwrite the all fiber class sizes to hold just the fast fibers
		#write fiber sizes to separate files
		fast_out.write(tmp_img+","+",".join([str(n) for n in all.sizes])+"\n")	
		slow_out.write(tmp_img+","+",".join([str(n) for n in slow.sizes])+"\n")		
	
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
	subtracted_area = im.size[0]*im.size[1] - fiber_area
	if main_structure_area > 0:
		subtracted_area = main_structure_area - fiber_area
	out_array = [tmp_img,str(fiber_number),str(mean_fiber),str(fiber_variance),str(fiber_area),str(subtracted_area),str(num_fast),str(num_slow)]
	out.write(",".join(out_array)+"\n")
	fs = [str(x) for x in fiber_sizes]
	#fs = map(str,fiber_sizes)
	out2.write(tmp_img+","+",".join(fs)+"\n")
	
print "time to process (seconds): " + str(time.time() - start)	
	
	
	