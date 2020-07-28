import os
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
from cv2 import cv2
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import pandas as pd
from scipy.optimize import curve_fit


def func_k(indep_var, k):
	"""
	Return the function used for the estimation of the error item of k.

	It has the form:

	$gap = 2\pi \times l - 2\pi \times [(R - \frac{k}{2\pi}) \times (sin(\pi - \frac{l}{R}))]$

	Arguments:

		indep_var (numpy.ndarray): shape is [n,2]. The 1st column is the array of the concentric circle radius. The 2nd column is the array of the ground truth radius `r`.

		k (variable): the variable to be estimated.
	"""
	l,r = indep_var[:,0], indep_var[:,1]
	return 2*np.pi*l - 2*np.pi*((r-k/(2*np.pi))*np.sin(np.pi - l/r))


def func_r(indep_var, r):
	"""
	Return the function used for the estimation of the radius.

	It has the form:
	
	$gap = 2\pi \times l - 2\pi \times [(R - \frac{k}{2\pi}) \times (sin(\pi - \frac{l}{R}))]$

	Arguments:

		indep_var (numpy.ndarray): shape is [n,2]. The 1st column is the array of the concentric circle radius. The 2nd column is the array of the error term `k`.

		r (variable): the variable to be estimated.
	"""
	l,k = indep_var[:,0], indep_var[:,1]
	return 2*np.pi*l - 2*np.pi*((r-k/(2*np.pi))*np.sin(np.pi - l/r))


def get_k_value(df_file):
	"""
	Return the estimate of the k value.

	Arguments:

		df_file (str): path to the 'csv' table file containing the estimated k values.

	Returns:

		float: mean of the k values estimated from all the sample points.

		float: median of the k values estimated from all the sample points.

		float: mean of the k values estimated from the middle sample points.

		float: median of the k values estimated from the middle sample points.

		float: micrometers per pixel.
	"""
	df = pd.read_csv(df_file)
	ks_all = df['Fitted_k_pixel_all'].to_numpy()
	ks_mid = df['Fitted_k_pixel_middle'].to_numpy()
	return np.mean(ks_all), np.median(ks_all), np.mean(ks_mid), np.median(ks_mid)


def get_files_by_format(fpath, format):
	"""
	Get all the file names with the specific format `format` in directory `fpath`.
	
	Arguments:

		fpath (str): the files' folder path.

		format (str): the specific format of the interested file without the prefixed dot.

	Returns:

		list of str: all the file names.

		int: the total number of the interested files.
	"""
	s_files = os.listdir(fpath)
	files = [x for x in s_files if x.split('.')[-1].upper()==format.upper()]
	fnum = len(files)
	return files, fnum


def make_sure_dir_exist(dir_):
	"""
	Make sure the directory exists.

	Arguments:

		dir_ (str): path to the directory.
	"""
	if not os.path.isdir(dir_):
		os.makedirs(dir_)


def get_bmask(img, num_open):
	"""
	Generate tissue mask by thresholding the 'S' channel in the HSV space.

	Arguments:

		img (numpy.ndarray): RGB image.

		num_open (int): number of the open operations.

	Returns:

		mask (numpy.ndarray): binary mask of the tissue region.
	"""
	thresh = 30
	# RGB
	img_hsv = cv2.cvtColor(img, code=cv2.COLOR_RGB2HSV)
	# find non-white pixels by thresholding on the channel 'S'
	ind = img_hsv[:,:,1]>thresh
	# construct binary mask
	h, w, _ = img.shape
	mask = np.zeros((h, w)).astype(np.uint8)
	mask[ind] = 255
	# remove noise
	for i in range(num_open):
		mask = cv2.morphologyEx(mask, op=cv2.MORPH_OPEN,
								kernel=cv2.getStructuringElement(shape=cv2.MORPH_RECT, ksize=(3,3)),
								iterations=5)
		mask = cv2.morphologyEx(mask, op=cv2.MORPH_CLOSE,
								kernel=cv2.getStructuringElement(shape=cv2.MORPH_RECT, ksize=(3,3)),
								iterations=7)
	return mask


def get_tissue_inds(img):
	"""Get the indices limits of the tissue region without the white margin.
	
	Arguments:

		img (numpy.ndarray): RGB image.

	Returns:

		row_min, row_max, col_min, col_max (int): indices limits of the tissue region.
	"""
	mask = get_bmask(img, 1)
	rows,cols = np.nonzero(mask)
	row_min, row_max = np.amin(rows), np.amax(rows)
	col_min, col_max = np.amin(cols), np.amax(cols)
	return row_min, row_max, col_min, col_max


def get_coords(ann_file, tag_txt):
	"""
	Return annotation coordinates.

	Arguments:

		ann_file (str): path to annotation file.

		tag_txt (str): label/tag name of the annotation.

	Returns:

		pts (dict): dictionary of coordinates. Keys are the Id of the annotation, values are the coordinates 2D list of (x,y).
	"""
	tree = ET.parse(ann_file)
	pts = {} 
	for ann in tree.iter('Annotation'):
		if ann.get('Name')==tag_txt:
			for reg in ann.iter('Region'):
				Id = int(float(reg.get('Id')))
				for vertex in reg.iter('Vertex'):
					pt = [int(float(vertex.get('X'))), int(float(vertex.get('Y')))]
				pts[Id] = pt
	return pts


def get_circle_coords(origin, r):
	"""
	Return the coordinates of the circle.

	Arguments:

		origin (list of int): coordinate of the origin, a list, [x,y].

		r (int): radius of the circle.

	Returns:

		xyu (list of tuple): coordinates of the circle. Each tuple is of the format (x,y).
	"""
	phi = np.linspace(0, 2*np.pi, round(r)*10)
	x = np.round(r * np.cos(phi) + origin[0])
	y = np.round(r * np.sin(phi) + origin[1])
	xy = [(i,j) for i,j in zip(x,y)]
	# not ordered unique coordinates of the circle
	xyu = list(set(xy))
	return xyu


def get_gap(mask, origin, l):
	"""
	Return the length of the gap between each lobe.

	Arguments:

		mask (numpy.ndarray): binary mask of the tissue.

		origin (list of int): coordinate of the origin, a list, [x,y].

		l (list of int): list of the concentric circle radius.

	Returns:

		numpy.ndarray: list of the gap length. It has the same shape as `l`.
	"""
	# 255 indicates gap
	mask = 255-mask
	h, w = mask.shape
	gap = []
	for r in l:
		xy = get_circle_coords(origin, r)
		xy = np.array(xy, dtype=int)
		# in case the concentric circle is outside of the image
		if np.amax(xy[:,1])>=h:
			xy[:,1][xy[:,1]>=h] = h-1
		if np.amax(xy[:,0])>=w:
			xy[:,0][xy[:,0]>=w] = w-1
		gap_pix = mask[xy[:,1], xy[:,0]]
		gap.append(np.sum(gap_pix)/255.)
	return np.array(gap)


def prepare_data(img_file, ann_file):
	"""
	Prepare data.
	
	Arguments:

		img_file (str): path to the image file.

		ann_file (str): path to the annotation file.

	Returns:

		l_array (numpy.ndarray): list of the concentric circle radius.

		gap (numpy.ndarray): list of the gap length. It has the same shape as `l_array`.

		l_for_theta (int): the concentric circle radius used to estimate the `theta`.

		origin (list of int): coordinate of the origin, a list, [x,y].
	"""
	# get tissue binary mask
	img = cv2.imread(img_file)[:,:,::-1]
	bmask = get_bmask(img, 1)
	# read annotation coordinates
	median_pts = get_coords(ann_file, 'Median_Points')
	origin = get_coords(ann_file, 'Origin')[1]
	# compute lobe median lenght
	lobe_medians = []
	for _,val in median_pts.items():
		distance = np.linalg.norm(np.array(val) - np.array(origin))
		lobe_medians.append(distance)
	# select concentric circle radius used to estimate theta
	l_for_theta = np.amax(np.array(lobe_medians))
	# list of concentric circle radius
	l_array = np.linspace(20, np.amin(np.array(lobe_medians)), 500)
	# list of gap length
	gap = get_gap(bmask, origin, l_array)
	return l_array, gap, l_for_theta, origin


def sel_mid_data(rm_first_percent=0.1, rm_last_percent=0.1, rm_zero_gap=True):
	"""
	Select the middle data points from the input data arrays.

	Arguments:

		rm_first_percent (float): Percentage of the first sample points to be removed. Default: 0.1

		rm_last_percent (float): Percentage of the last sample points to be removed. Default: 0.1

		rm_zero_gap (bool): If true, remove the sample points with the zero gap. Default: True

	Call arguments:

		l_array (numpy.ndarray): list of the concentric circle radius.

		gap (numpy.ndarray): list of the gap length. It has the same shape as `l_array`.

	Returns:

		numpy.ndarray: list of the selected concentric circle radius.

		numpy.ndarray: list of the selected gap length. It has the same shape as the 1st return.

	Raises:

		ValueError: if the value of `rm_first_percent` and `rm_last_percent` are not between 0 and 1.
	"""
	if rm_last_percent<0 or rm_last_percent>=1:
		raise ValueError("Invalid value for `rm_last_percent`, must be between 0 and 1.")
	if rm_first_percent<0 or rm_first_percent>=1:
		raise ValueError("Invalid value for `rm_first_percent`, must be between 0 and 1.")
	def get_data(l_array,gap):
		pt_num = len(gap)
		# remove zero gaps
		if rm_zero_gap:
			sel_index = gap>0
		else:
			sel_index = np.array([True for _ in range(pt_num)])
		# remove the first tenth and the last fifth points
		sel_index[:int(pt_num*rm_first_percent)] = 0
		sel_index[int(pt_num*(1-rm_last_percent)):] = 0
		return l_array[sel_index], gap[sel_index]
	return get_data


def show_fit_both_r(im_file, func, indep_var_all, indep_var_mid, r_all, r_mid, gap_all, gap_mid, theta_all, theta_mid, origin, um_per_pixel, save_file=None):
	"""
	Show the fitting curves for both the all-sample-points and the middle-sample-points during the estimation of the radius.
	
	Arguments:

		im_file (str): path to the image file.

		func (function): the function to estimate.

		indep_var_all (numpy.ndarray): independent variables of the function for the all-sample-points. Refer to the `func_r` for more informations.

		indep_var_mid (numpy.ndarray): independent variables of the function for the middle-sample-points. Refer to the `func_r` for more informations.

		r_all (float): the estimated radius using the all-sample-points.

		r_mid (float): the estimated radius using the middle-sample-points.

		gap_all (numpy.ndarray): list of the measured gap length for the all-sample-points.

		gap_mid (numpy.ndarray): list of the measured gap length for the middle-sample-points.

		theta_all (float): estimated theta using the all-sample-points.

		theta_mid (float): estimated theta using the middle-sample-points.

		origin (list of int): coordinate of the origin, a list, [x,y].

		um_per_pixel (float): micrometer per pixel.

		save_file (str): path to the the save image file. If None, do not save. Default: None.
	"""
	# fitted gap length
	y_all = func(indep_var_all, r_all)
	y_mid = func(indep_var_mid, r_mid)
	# estimated diameter
	diameter_all = r_all*2*um_per_pixel/1000
	diameter_mid = r_mid*2*um_per_pixel/1000
	# concentric circle radius
	l_array_all = indep_var_all[:,0]
	l_array_mid = indep_var_mid[:,0]
	# plot figure
	fig = plt.figure(constrained_layout=True)
	gs = GridSpec(2, 2, figure=fig)
	fig.set_size_inches(16,8)
	name = Path(im_file).stem
	fig.suptitle(name, size='xx-large')
	# Axis0: original image with annotations indicating the largest concentric circle and the origin
	ax0 = fig.add_subplot(gs[:,0])
	img = cv2.imread(im_file)[:,:,::-1]
	row_min, row_max, col_min, col_max = get_tissue_inds(img)
	ax0.imshow(img[row_min:row_max+1, col_min:col_max+1])
	xy = get_circle_coords(origin, l_array_all[-1])
	X = [i[0]-col_min for i in xy]
	Y = [i[1]-row_min for i in xy]
	ax0.scatter(X, Y, c='r', s=1)
	ox, oy = origin
	ox, oy = ox-col_min, oy-row_min
	ax0.scatter(ox, oy, c='r', s=100)
	ax0.set_axis_off()
	# Axis1: fitting curves using the all-sample-points
	ax1 = fig.add_subplot(gs[0,1:])
	ax1.plot(l_array_all, gap_all, 'b-', label='Data')
	ax1.plot(l_array_all, y_all, 'r-', label='Fit: diameter={:.2f} mm'.format(diameter_all))
	ax1.set_xlabel('Radius of Concentric Circles', size='x-large')
	ax1.set_ylabel('Gaps', size='x-large')
	ax1.set_title('Estimated sphere diameter from all points: {:.2f} mm, Estimated $\\theta$: {:.1f} degree'.format(diameter_all, theta_all/np.pi*180), size='x-large')
	ax1.legend()
	# Axis2: fitting curves using the middle-sample-points
	ax2 = fig.add_subplot(gs[1,1:])
	ax2.plot(l_array_all, gap_all, 'b-', label='Data')
	ax2.plot(l_array_mid, y_mid, 'r-', label='Fit: diameter={:.2f} mm'.format(diameter_mid))
	ax2.set_xlabel('Radius of Concentrical Circles', size='x-large')
	ax2.set_ylabel('Gaps', size='x-large')
	ax2.set_title('Estimated sphere diameter from middle points: {:.2f} mm, Estimated $\\theta$: {:.1f} degree'.format(diameter_mid, theta_mid/np.pi*180), size='x-large')
	ax2.legend()
	# save figure
	if save_file!=None:
		make_sure_dir_exist(str(Path(save_file).parent))
		fig.savefig(save_file, dpi=100)
	return


def estimate_paras_r(img_dir, res_file, k, func, um_per_pixel, draw_fit_curve=True, save_fig_file=None, sel_pt_method=sel_mid_data()):
	"""
	Estimate the parameters of the radius.
	
	Arguments:

		img_dir (str): path to the directory storing all the source images and their annotations files.

		res_file (str): path to the table file saving the fitted results. Its format is csv.

		k (float): the error term `k`. Refer to `get_k_value` for more information.

		func (function): the function to estimate.

		um_per_pixel (float): micrometer per pixel.

		draw_fit_curve (bool): whether to draw the fitting curves. Default: True.

		save_fig_file (str): path to the the save image file. If None, do not save. Default: None.

		sel_pt_method (function): the method to select partial sample points. Default: sel_mid_data.

	Raises:

		FileNotFoundError: if the image file ('tif' or 'jpg') does not exist.
	"""
	# get annotation names
	ann_names, ann_num = get_files_by_format(img_dir, 'xml')
	# table to store results
	res = []
	for count, ann_name in enumerate(ann_names, start=1):
		name = Path(ann_name).stem
		imFileTiff = os.path.join(img_dir, name+'.tif')
		imFileJpg = os.path.join(img_dir, name+'.jpg')
		if os.path.isfile(imFileTiff):
			im_file = imFileTiff
		elif os.path.isfile(imFileJpg):
			im_file = imFileJpg
		else:
			raise FileNotFoundError("File {} or {} not found.".format(imFileTiff, imFileJpg))
		ann_file = os.path.join(img_dir, name+'.xml')
		l_array_all, gap_all, l_for_theta, origin = prepare_data(im_file, ann_file)
		l_array_mid, gap_mid = sel_pt_method(l_array_all, gap_all)
		# independent variables for the function
		k_array_all = np.ones_like(l_array_all)*k
		indep_var_all = np.vstack((l_array_all,k_array_all)).T
		# to fit radius, must set a initial value `p0` much greater than 0, becuase the fitting has a local minimum near 0
		popt_all, pcov_all = curve_fit(f=func, xdata=indep_var_all, ydata=gap_all, p0=1000)
		# independent variables for the function
		k_array_mid = np.ones_like(l_array_mid)*k
		indep_var_mid = np.vstack((l_array_mid,k_array_mid)).T
		# to fit radius, must set a initial value `p0` much greater than 0, becuase the fitting has a local minimum near 0
		popt_mid, pcov_mid = curve_fit(f=func, xdata=indep_var_mid, ydata=gap_mid, p0=1000)
		# estimates
		r_all = popt_all[0]
		r_mid = popt_mid[0]
		perr_all = np.sqrt(np.diag(pcov_all))
		perr_mid = np.sqrt(np.diag(pcov_mid))
		theta_all = np.pi - l_for_theta/r_all
		theta_mid = np.pi - l_for_theta/r_mid
		res.append([name, r_all*2*um_per_pixel/1000, r_mid*2*um_per_pixel/1000, theta_all/np.pi*180, theta_mid/np.pi*180])
		if draw_fit_curve:
			show_fit_both_r(im_file, func, indep_var_all, indep_var_mid, r_all, r_mid, gap_all, gap_mid, theta_all, theta_mid, origin, um_per_pixel, save_file=save_fig_file)
		print("{}/{}: {}".format(count, ann_num, name))
	df = pd.DataFrame(data=res, columns=['Name', 'Diameter-All(mm)', 'Diameter-Middle(mm)', 'Theta-All(degree)', 'Theta-Middle(degree)'])
	make_sure_dir_exist(str(Path(res_file).parent))
	df.to_csv(res_file, index=False)
	print("Complete.\n")
	return


def show_fit_both_k(im_file, func, indep_var_all, indep_var_mid, k_all, k_mid, gap_all, gap_mid, r_av_pix, origin, um_per_pixel, save_file=None):
	"""
	Show the fitting curves for both the all-sample-points and the middle-sample-points during the estimation of the error item 'k'.

	Arguments:

		im_file (str): path to the image file.

		func (function): the function to estimate.

		indep_var_all (numpy.ndarray): independent variables of the function for the all-sample-points. Refer to the `func_k` for more informations.

		indep_var_mid (numpy.ndarray): independent variables of the function for the middle-sample-points. Refer to the `func_k` for more informations.

		k_all (float): the estimated error item using the all-sample-points.

		k_mid (float): the estimated error item using the middle-sample-points.

		gap_all (numpy.ndarray): list of the measured gap length for the all-sample-points.

		gap_mid (numpy.ndarray): list of the measured gap length for the middle-sample-points.

		r_av_pix (float): the ground-truth radius in pixels.

		origin (list of int): coordinate of the origin, a list, [x,y].

		um_per_pixel (float): micrometer per pixel.

		save_file (str): path to the the save image file. If None, do not save. Default: None.
	"""
	# fitted gap length
	y_all = func(indep_var_all, k_all)
	y_mid = func(indep_var_mid, k_mid)
	# concentric circle radius
	l_array_all = indep_var_all[:,0]
	l_array_mid = indep_var_mid[:,0]
	# error item in pixels
	k_div_2pi_all = k_all/2/np.pi*um_per_pixel/1000
	k_div_2pi_mid = k_mid/2/np.pi*um_per_pixel/1000
	# ground truth diameter in mm
	diameter = r_av_pix*2*um_per_pixel/1000
	# plot figure
	fig = plt.figure(constrained_layout=True)
	gs = GridSpec(2, 2, figure=fig)
	fig.set_size_inches(16,8)
	name = Path(im_file).stem
	fig.suptitle(name, size='xx-large')
	# Axis0: original image with annotations indicating the largest concentric circle and the origin
	ax0 = fig.add_subplot(gs[:,0])
	img = cv2.imread(im_file)[:,:,::-1]
	row_min, row_max, col_min, col_max = get_tissue_inds(img)
	ax0.imshow(img[row_min:row_max+1, col_min:col_max+1])
	xy = get_circle_coords(origin, l_array_all[-1])
	X = [i[0]-col_min for i in xy]
	Y = [i[1]-row_min for i in xy]
	ax0.scatter(X, Y, c='r', s=1)
	ox, oy = origin
	ox, oy = ox-col_min, oy-row_min
	ax0.scatter(ox, oy, c='r', s=100)
	ax0.set_axis_off()
	# Axis1: fitting curves using the all-sample-points
	ax1 = fig.add_subplot(gs[0,1:])
	ax1.plot(l_array_all, gap_all, 'b-', label='Data')
	ax1.plot(l_array_all, y_all, 'r-', label='Ground truth diameter: {:.2f} mm, fitted $k/2\pi$: {:.5f}'.format(diameter, k_div_2pi_all))
	ax1.set_xlabel('Radius of Concentrical Circles', size='x-large')
	ax1.set_ylabel('Gaps', size='x-large')
	ax1.set_title('Estimated error $k/2\pi$ from all points: {:.5f} mm'.format(k_div_2pi_all), size='x-large')
	ax1.legend()
	# Axis2: fitting curves using the middle-sample-points
	ax2 = fig.add_subplot(gs[1,1:])
	ax2.plot(l_array_all, gap_all, 'b-', label='Data')
	ax2.plot(l_array_mid, y_mid, 'r-', label='Ground truth diameter: {:.2f} mm, fitted $k/2 \pi$: {:.5f}'.format(diameter, k_div_2pi_mid))
	ax2.set_xlabel('Radius of Concentrical Circles', size='x-large')
	ax2.set_ylabel('Gaps', size='x-large')
	ax2.set_title('Estimated error $k/2\pi$ from middle points: {:.5f} mm'.format(k_div_2pi_mid), size='x-large')
	ax2.legend()
	# save figure
	if save_file!=None:
		make_sure_dir_exist(str(Path(save_file).parent))
		fig.savefig(save_file, dpi=100)
	return


def estimate_paras_k(gt_file, res_file, func, draw_fit_curve=True, sel_pt_method=sel_mid_data()):
	"""
	Estimate the parameters.
	
	Arguments:

		gt_file (str): path to the table file containing all the ground truth of the measured radius and the image path.

		res_file (str): path to the table file saving the fitted results.

		func (function): the function to estimate.

		draw_fit_curve (bool): whether to draw the fitting curves. Default: True.

		sel_pt_method (function): the method to select partial sample points. Default: sel_mid_data.
	"""
	# results table
	res_all = []
	res_mid = []
	res_all_pix = []
	res_mid_pix = []
	# read ground truth table
	df_gt = pd.read_csv(gt_file)
	df_gt = df_gt.dropna(axis='index', how='any', subset=['Av', 'FilePath'], inplace=False)
	samp_num = len(df_gt)
	for count, ind in enumerate(df_gt.index, 1):
		img_file = df_gt.at[ind, 'FilePath']
		name = Path(img_file).stem
		ann_file = str(Path(img_file).with_suffix('.xml'))
		if not os.path.isfile(ann_file):
			continue
		# image preeprocessing to get `l` and `gap`
		l_array_all, gap_all, l_for_theta, origin = prepare_data(img_file, ann_file)
		l_array_mid, gap_mid = sel_pt_method(l_array_all, gap_all)
		# ground truth radius
		um_per_pixel = df_gt.at[ind, 'um_per_pixel']
		diameter_av = df_gt.at[ind, 'Av']
		r_av_pix = diameter_av*1000/2/um_per_pixel
		# fit the function to estimate `k`
		r_array_all = np.ones_like(l_array_all)*r_av_pix
		indep_var_all = np.vstack((l_array_all,r_array_all)).T
		popt_all, pcov_all = curve_fit(func, indep_var_all, gap_all)
		r_array_mid = np.ones_like(l_array_mid)*r_av_pix
		indep_var_mid = np.vstack((l_array_mid,r_array_mid)).T
		popt_mid, pcov_mid = curve_fit(func, indep_var_mid, gap_mid)
		# get `k` and the estimation error
		k_all = popt_all[0]
		k_mid = popt_mid[0]
		perr_all = np.sqrt(np.diag(pcov_all))
		perr_mid = np.sqrt(np.diag(pcov_mid))
		# k in mm
		k_div_2pi_all = k_all/2/np.pi*um_per_pixel/1000
		k_div_2pi_mid = k_mid/2/np.pi*um_per_pixel/1000
		# save results in table
		res_all.append(k_div_2pi_all)
		res_mid.append(k_div_2pi_mid)
		res_all_pix.append(k_all)
		res_mid_pix.append(k_mid)
		# print output
		print('{}/{}    '.format(count, samp_num), name+': ')
		print(
			'Fitted (k/2pi) by all points is: {:.5f} mm (standard deviation: {:.5f} mm), '.format(k_div_2pi_all, perr_all[0]/2/np.pi*um_per_pixel/1000) + 
			'Fitted (k/2pi) by middle points is: {:.5f} mm (standard deviation: {:.5f} mm)'.format(k_div_2pi_mid, perr_mid[0]/2/np.pi*um_per_pixel/1000)
		)
		# Draw figures
		if draw_fit_curve:
			show_fit_both_k(img_file, func, indep_var_all, indep_var_mid, k_all, k_mid, gap_all, gap_mid, r_av_pix, origin, um_per_pixel, save_file=None)
	# insert results into the table
	df_gt.insert(loc=len(df_gt.columns), column='Fitted_k_div_2pi_all(mm)', value=res_all, allow_duplicates=False)
	df_gt.insert(loc=len(df_gt.columns), column='Fitted_k_div_2pi_middle(mm)', value=res_mid, allow_duplicates=False)
	df_gt.insert(loc=len(df_gt.columns), column='Fitted_k_pixel_all', value=res_all_pix, allow_duplicates=False)
	df_gt.insert(loc=len(df_gt.columns), column='Fitted_k_pixel_middle', value=res_mid_pix, allow_duplicates=False)
	# save results into spreadsheet file
	make_sure_dir_exist(str(Path(res_file).parent))
	df_gt.to_csv(res_file, index=False)
	print("Complete.\n")
	return

