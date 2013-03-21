from PIL import Image
import sys
import cv
import numpy as np
import scipy.ndimage as si
import scipy.spatial as sp
import argparse
import mahotas as mh


class analyze_all:

	def __init__(self, vis_pix, args, im, rgb_pix):
		self.adaptive_threshold(im, args)
		self.sizes = self.quantify_fibers(args, vis_pix, rgb_pix)
				
	def adaptive_threshold(self, im, args): #do adaptive thresholding to isolate features
		cv_im = cv.CreateImageHeader(im.size, cv.IPL_DEPTH_8U, 1)
		cv.SetData(cv_im,im.tostring(),im.size[0])	
		cv.AdaptiveThreshold(cv_im,cv_im,255.0,adaptive_method=cv.CV_ADAPTIVE_THRESH_MEAN_C,thresholdType=cv.CV_THRESH_BINARY_INV, blockSize=args.box_size, param1=args.offset) 
		cv.ConvertScale(cv_im,cv_im,-1,255)
		cv.Erode(cv_im,cv_im,iterations=int(args.erosion_steps))
		cv.SaveImage("tmp.png",cv_im)		
		
	def quantify_fibers(self, args, vis_pix, rgb_pix): #label feature with number (1)
		im2 = Image.open("tmp.png")
		#im2.show()
		#sys.exit(0)
		pix = im2.getdata()
		labeled_array, num_features = si.label(np.array(pix).reshape(im2.size[1],im2.size[0]),structure=[[1,1,1],[1,1,1],[1,1,1]])
		label_hist = si.measurements.histogram(labeled_array,1,num_features,num_features)
		condition = (label_hist > args.min_fiber_size) & (label_hist < args.max_fiber_size)
		fiber_sizes = label_hist[condition]
		fiber_labels = np.arange(1,num_features+1)[condition]
		keep_fibers = dict.fromkeys(fiber_labels,0)
		for index,i in enumerate(labeled_array.flatten()):
			if i in keep_fibers:
				vis_pix[index] = 1
				rgb_pix[index] = (102, 204, 255)
		return fiber_sizes