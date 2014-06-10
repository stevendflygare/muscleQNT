from PIL import Image
import sys
import cv
import numpy as np
import scipy.ndimage as si
import scipy.spatial as sp
import argparse
import mahotas as mh
import collections


class analyze_slow:

	def __init__(self, vis_pix, args, im_co, rgb_pix):
		thresholded_copix = self.global_threshold(im_co, vis_pix, args)
		#tmp = Image.new("L",im_co.size)
		#tmp.putdata(thresholded_copix)
		#tmp.show()
		#self.sizes = self.quantify_fibers(args, vis_pix, thresholded_copix, im_co, rgb_pix)		
		self.sizes = self.quantify_fibers2(args, vis_pix, thresholded_copix, im_co, rgb_pix)		
						
				
	def global_threshold(self, im_co, vis_pix, args): #do global thresholding to isolate slow fibers
		co_pix = np.array(im_co.getdata(),dtype=np.uint8)
		co_pix = co_pix*vis_pix
		T_otsu = mh.otsu(co_pix.reshape(im_co.size[1],im_co.size[0]))
		thresholded_copix = (co_pix*(co_pix > T_otsu))
		#thresholded_copix = si.grey_erosion(np.array(thresholded_copix).reshape(im_co.size[1],im_co.size[0]), size=(3,3))
		thresholded_copix = si.grey_closing(np.array(thresholded_copix).reshape(im_co.size[1],im_co.size[0]), size=(10,10))
		#thresholded_copix = si.grey_closing(thresholded_copix, size=(3,3))
		return thresholded_copix
		
	def quantify_fibers(self, args, vis_pix, thresholded_copix, im_co, rgb_pix): #label feature with number (2)
		co_labeled_array, co_num_features = si.label(thresholded_copix.reshape(im_co.size[1],im_co.size[0]),structure=[[1,1,1],[1,1,1],[1,1,1]])
		#im = Image.new("L", im_co.size)
		#im.putdata(thresholded_copix.flatten())
		#im.show()
		together_labeled_array = co_labeled_array.flatten()*(co_labeled_array.flatten() > 0)
		co_label_hist = si.measurements.histogram(together_labeled_array,1,np.max(together_labeled_array),np.max(together_labeled_array))
		condition = (co_label_hist > args.min_fiber_size) & (co_label_hist < args.max_fiber_size)
		co_fiber_sizes = co_label_hist[condition]
		#print co_fiber_sizes
		#print len(co_fiber_sizes)
		co_fiber_labels = np.arange(1,np.max(together_labeled_array)+1)[condition]
		colocated_hash = dict.fromkeys(co_fiber_labels,0)
		for index,i in enumerate(co_labeled_array.flatten()):
			if i in colocated_hash:
				vis_pix[index] = 2
				rgb_pix[index] = (255, 255, 0)
		return co_fiber_sizes
		
	def quantify_fibers2(self, args, vis_pix, thresholded_copix, im_co, rgb_pix): #label feature with number (2)
		vis_pix_matrix = vis_pix.reshape(im_co.size[1],im_co.size[0])
		all_array, num_features = si.label(vis_pix_matrix,structure=[[1,1,1],[1,1,1],[1,1,1]])
		co_fiber_sizes = []
		for fiber_slice in si.find_objects(all_array):
			fiber_pix = (thresholded_copix[fiber_slice].flatten()) > 0
			fiber_types = collections.Counter(all_array[fiber_slice].flatten())
			max_type = None
			if len(fiber_types) > 2:
				max_count = 0
				for ft in fiber_types:
					if ft > 0 and fiber_types[ft] > max_count:
						max_type = ft
						max_count = fiber_types[ft]
			if float(np.sum(fiber_pix))/float(fiber_pix.size) > .3: #if > 30% of the pixels are colocalized call entire fiber slow
				slow_fiber_size = 0
				for i in xrange(fiber_slice[0].start,fiber_slice[0].stop):
					for j in xrange(fiber_slice[1].start,fiber_slice[1].stop):
						if vis_pix_matrix[i][j] == 1:
							if max_type:
								if all_array[i][j] == max_type:
									vis_pix_matrix[i][j] = 2
									rgb_pix[i*im_co.size[0] + j] = (255, 255, 0)
									slow_fiber_size += 1
							else:
								vis_pix_matrix[i][j] = 2
								rgb_pix[i*im_co.size[0] + j] = (255, 255, 0)
								slow_fiber_size += 1
				co_fiber_sizes.append(slow_fiber_size)				
		vis_pix = vis_pix_matrix.flatten()
		return np.array(co_fiber_sizes)
						