from PIL import Image
import sys
import cv
import numpy as np
import scipy.ndimage as si
import scipy.spatial as sp
import argparse
import mahotas as mh


class analyze_slow:

	def __init__(self, vis_pix, args, im_co, rgb_pix):
		thresholded_copix = self.global_threshold(im_co, args)
		self.sizes = self.quantify_fibers(args, vis_pix, thresholded_copix, im_co, rgb_pix)		
				
	def global_threshold(self, im_co, args): #do global thresholding to isolate slow fibers
		co_pix = np.array(im_co.getdata(),dtype=np.uint8)
		T_otsu = mh.otsu(co_pix.reshape(im_co.size[1],im_co.size[0]))
		thresholded_copix = (co_pix*(co_pix > T_otsu))
		thresholded_copix = si.grey_erosion(np.array(thresholded_copix).reshape(im_co.size[1],im_co.size[0]), size=(3,3))
		thresholded_copix = si.grey_closing(thresholded_copix, size=(10,10))
		#thresholded_copix = si.grey_closing(thresholded_copix, size=(3,3))
		return thresholded_copix
		
	def quantify_fibers(self, args, vis_pix, thresholded_copix, im_co, rgb_pix): #label feature with number (2)
		co_labeled_array, co_num_features = si.label(thresholded_copix.reshape(im_co.size[1],im_co.size[0]),structure=[[1,1,1],[1,1,1],[1,1,1]])
		together_labeled_array = co_labeled_array.flatten()*(co_labeled_array.flatten() > 0)
		co_label_hist = si.measurements.histogram(together_labeled_array,1,np.max(together_labeled_array),np.max(together_labeled_array))
		condition = (co_label_hist > args.min_fiber_size) & (co_label_hist < args.max_fiber_size)
		co_fiber_sizes = co_label_hist[condition]
		co_fiber_labels = np.arange(1,np.max(together_labeled_array)+1)[condition]
		colocated_hash = dict.fromkeys(co_fiber_labels,0)
		for index,i in enumerate(co_labeled_array.flatten()):
			if i in colocated_hash:
				vis_pix[index] = 2
				rgb_pix[index] = (255, 255, 0)
		return co_fiber_sizes