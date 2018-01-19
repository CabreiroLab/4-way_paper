import os
import sys
import tifffile as tiff
import textwrap as tw
import csv
import pandas
import matplotlib

matplotlib.use('Agg') #TkAgg
#matplotlib.use('MacOSX') #TkAgg

import matplotlib.pyplot as plt


#fig, axes = plt.subplots(nrows=8, ncols=12, figsize=(72, 47), dpi=300)

#%matplotlib osx
#plt.interactive(False)


import numpy as np
import itertools as IT
import time


import skimage
from skimage import io as skio
from skimage import color
from skimage import data
from skimage import img_as_float
from skimage import filters
from skimage.util.dtype import dtype_range
from skimage.util import img_as_ubyte
from skimage.feature import canny
from skimage.filters import sobel
from skimage.morphology import disk, opening, dilation, square, watershed
from skimage.morphology import erosion, white_tophat, black_tophat, closing
from skimage import exposure


from skimage.filters import rank
from skimage.filters.rank import median, mean
from skimage import measure
from skimage.measure import label, regionprops

#from skimage.filters import threshold_otsu, threshold_adaptive

#from skimage.feature import peak_local_max


#from scipy import misc


import scipy.optimize as opt
from scipy import ndimage as ndi


def filename(ifile):
	if ifile.split('.')[0] == '':
		ipat = ''
		iname = ''
		itype = ifile.split('.')[1]
	else:
		if "\\" in ifile.split('.')[0]:
			sep = "\\"
		elif "/" in ifile.split('.')[0]:
			sep = "/"
		else:
			ipat = ''
			iname = ifile.split('.')[0]
			itype = ifile.split('.')[1]
			return ipat, iname, itype
		allpath = ifile.split('.')[0]
		iname = allpath.split(sep)[-1]
		ipath = allpath.split(sep)[:-1]
		ipat = '/'.join(ipath)
		itype = ifile.split('.')[1]
	return ipat, iname, itype


def slice_list(input, size):
	input_size = len(input)
	slice_size = input_size / size
	remain = input_size % size
	result = []
	iterator = iter(input)
	for i in range(size):
		result.append([])
		for j in range(slice_size):
			result[i].append(iterator.next())
		if remain:
			result[i].append(iterator.next())
			remain -= 1
	return result


class NestedDict(dict):
	def __getitem__(self, key):
		if key in self: return self.get(key)
		return self.setdefault(key, NestedDict())


def numerize(s):
	try:
		if s == 'NAN':
			return s
		float(s)
		if float(s).is_integer():
			return int(float(s))
		elif float(s) == 0:
			return float(s)
		else:
			return float(s)

	except ValueError:
		return s


def indx2well(ind, start=0, rowln=12):
	indt = ind - start
	rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
	row = indt // rowln
	col = indt - row * rowln + 1
	well = '{}{}'.format(rows[row], col)
	return well


def readmet(ifile):
	print ifile
	nutrients = NestedDict()
	rdr = csv.reader(open(ifile, 'r'), delimiter=',')
	data = [ln for ln in rdr]
	headers = data[0]
	headin = {hd: headers.index(hd) for hd in headers}
	nec = ['Metabolite', 'EcoCycID', 'Plate', 'Well', 'Index']
	if all(n in headers for n in nec):
		print 'Necessary headers found!'
	else:
		print 'Missing essential headers in metabolites file!'
		print headers
		sys.exit(0)
	for ln in data[1:]:
		pl = str(numerize(ln[headin['Plate']])).strip().encode('ascii', 'ignore')
		wl = str(numerize(ln[headin['Well']])).strip().encode('ascii', 'ignore')
		for hd in headin.keys():
			nutrients[pl][wl][hd] = str(numerize(ln[headin[hd]])).strip().encode('ascii', 'ignore')
	return nutrients


def readdel(ifile):
	print ifile

	deletions = NestedDict()
	rdr = csv.reader(open(ifile, 'r'), delimiter=',')
	data = [ln for ln in rdr]
	headers = data[0]
	headin = {hd: headers.index(hd) for hd in headers}
	nec = ['Replicate', 'Plate', 'Index', 'File', 'X1', 'Y1', 'X2', 'Y2']  # ,'Well'
	if all(n in headers for n in nec):
		print 'Necessary headers found!'
	else:
		print 'Missing essential headers in deletions file!'
		print headers
		sys.exit(0)
	for ln in data[1:]:
		# print ln
		# Reading
		rep = ln[headin['Replicate']].strip().encode('ascii', 'ignore')
		pl = ln[headin['Plate']].strip().encode('ascii', 'ignore')
		fl = ln[headin['File']].strip().encode('ascii', 'ignore')
		indx = numerize(ln[headin['Index']])
		coords = [numerize(ln[headin[hdr]]) for hdr in ['X1', 'Y1', 'X2', 'Y2']]
		# Storing
		deletions[rep][pl][indx]['File'] = fl

		if 'Worms' in deletions[rep][pl][indx].keys():
			worms = deletions[rep][pl][indx]['Worms']
			worms.append(coords)
			deletions[rep][pl][indx]['Worms'] = worms
		else:
			deletions[rep][pl][indx]['Worms'] = [coords]

	return deletions



def wormdel(img, worms):
	for worm in delworms:
		x1, y1, x2, y2 = worm
		xs = [x1, x2]
		ys = [y1, y2]
		img[min(ys):max(ys), min(xs):max(xs), :] = 0
	return img


def readman(ifile):
	print ifile
	manual = NestedDict()
	rdr = csv.reader(open(ifile, 'r'), delimiter=',')
	data = [ln for ln in rdr]
	headers = data[0]
	headin = {hd: headers.index(hd) for hd in headers}
	nec = ['Replicate', 'Plate', 'Index', 'Threshold']  # ,'Well'
	if all(n in headers for n in nec):
		print 'Necessary headers found!'
	else:
		print 'Missing essential headers in deletions file!'
		print headers
		sys.exit(0)

	for ln in data[1:]:
		# print ln
		# Reading
		rep = ln[headin['Replicate']].strip().encode('ascii', 'ignore')
		pl = ln[headin['Plate']].strip().encode('ascii', 'ignore')
		indx = numerize(ln[headin['Index']])
		thrs = numerize(ln[headin['Threshold']])
		thrs = thrs if thrs != "" else 0
		manual[rep][pl][indx] = thrs if thrs != 0.0 else 0.0
	return manual




def freqtable(vector):
	my_series = pandas.Series(vector)
	counts = my_series.value_counts()
	ftable = [[key, value] for key, value in dict(counts).iteritems()]
	return ftable


def writecsv(data, ofile, sep=','):
	f = open(ofile, "wb")
	ofile = csv.writer(f, delimiter=sep)  # dialect='excel',delimiter=';', quotechar='"', quoting=csv.QUOTE_ALL
	for row in data:
		# row = [item.encode("utf-8") if isinstance(item, unicode) else str(item) for item in row]
		ofile.writerow(row)
	f.close()



def labeling3(imghsv, hthres,cthres,size):
	h = imghsv[:, :, 0]
	# s=imghsv[:,:,1]
	v = imghsv[:, :, 2]

	# Filtering
	dh = rank.median(h, disk(3)) / 255.0
	dv = rank.median(v, disk(3)) / 255.0
	# np.max(dh)
	# hmax=np.percentile(h,95)

	hf = dh.copy()
	hf[(hf < hthres).astype(bool)] = 0
	hf = hf / hthres
	# plt.imshow(hf)

	comb_r = hf * dv
	comb = opening(comb_r, selem=disk(3))

	markers = np.zeros_like(comb)
	# Mark background
	markers[comb == 0] = 1
	markers[comb > cthres] = 2

	elevation_map = sobel(v)
	segmentation = watershed(elevation_map, markers)
	segmentation = ndi.binary_fill_holes(segmentation - 1)
	labeled_worms, _ = ndi.label(segmentation)

	for w in list(np.unique(labeled_worms)):
		# print labeled_worms[labeled_worms==w].shape[0]
		if labeled_worms[labeled_worms == w].shape[0] < size:
			labeled_worms[labeled_worms == w] = 0

	return labeled_worms




def map2table(data):
	table = []
	for a in data.keys():
		for b in data[a].keys():
			for c in data[a][b].keys():
				row = [a, b, c, data[a][b][c]]
				table.append(row)
	return table


def gauss(x, p):  # p[0]==mean, p[1]==stdev
	return 1.0 / (p[1] * np.sqrt(2 * np.pi)) * np.exp(-(x - p[0]) ** 2 / (2 * p[1] ** 2))

def fitgauss(xr, yr):

	if xr[0]==0.0:
		xr=xr[1:]
		yr=yr[1:]

	diffs = [j - i for i, j in zip(xr[:-1], xr[1:])]
	# Renormalize to a proper PDF
	Y = yr / np.sum(yr[1:] * diffs)
	X = xr
	# Fit a guassian
	p0 = [0.3, 0.05]  # Inital guess is a normal distribution
	errfunc = lambda p, x, y: gauss(x, p) - y  # Distance to the target function
	p1, success = opt.leastsq(errfunc, p0[:], args=(X, Y))
	fit_mu, fit_stdev = p1
	FWHM = 2 * np.sqrt(2 * np.log(2)) * fit_stdev
	print "FWHM: {}, mu: {}, SD: {}".format(FWHM, fit_mu, fit_stdev)
	return X, Y, fit_mu, fit_stdev



os.chdir("/home/pnorv/Dropbox/Projects/2015-Metformin/Biolog_Met_NGM/Celegans")
matplotlib.rcParams.update({'font.size': 12})

sourceloc = "/home/pnorv/Dropbox/Projects/Biolog_Celegans_Metformin"
replicate = "Rep1"
plate = "PM1"

odir = '.'

nutrients = readmet('/home/pnorv/Dropbox/Projects/2015-Metformin/Biolog/Biolog_metabolites_EcoCyc.csv')

alldeletions = readdel('{}/Deletions.csv'.format(odir))


thresholds_man = readman('Thresholds_all_3.csv')

# ==============================================================================
# #@#96 well previews
# ==============================================================================
# rows=['A','B','C','D','E','F','G','H']

ttllen = 24
ttlsize = 24

minimage = 1
maximage = 96

maxpix = 0

data = []
thresholds=[]
# figures=['Original','Labeling']
#figures = ['Original']

figures = ['Labeling']

#plates = ['PM1']
plates=['PM1','PM2A','PM3B','PM4A']
#replicates = ['Rep2', 'Rep3', 'Rep4','Rep5','Rep6']
#replicates = ['Rep1']
#replicates = ['Rep5', 'Rep6']

replicates = ['Rep1', 'Rep2', 'Rep3', 'Rep4','Rep5', 'Rep6']

#np.set_printoptions(formatter={'float': lambda x: "{0:.3f}".format(x)})
np.set_printoptions(precision=3)

step = 0.001
levels = ["{0:.3f}".format(lvl) for lvl in np.arange(0, 1 + step, step)]


printfigures = True
indo = 0.0

sthres=1000
cthres=0.02

label = 'auto_threshold'

for replicate in replicates:

	platesel = plates[:]
	if replicate == 'Rep4':
		platesel=['PM1','PM3B','PM4A']
	elif replicate == 'Rep5':
		platesel=['PM2A']
	elif replicate == 'Rep6':
		platesel=['PM1','PM2A']

	for plate in platesel:
		mthres = thresholds_man[replicate[3]][plate]

		files = ["{}_{}.tif".format(plate, str(im).zfill(3)) for im in range(1, 96 + 1)]
		total = len(plates) * len(files) * len(replicates)
		figall = []
		axall = []
		if printfigures:
			for fi, fign in enumerate(figures):
				fig, axes = plt.subplots(nrows=8, ncols=12, figsize=(72, 47), dpi=300)
				fig.suptitle('{}-{}'.format(replicate, plate), fontsize=40)
				plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=2, wspace=0.5, hspace=0.1)
				figall.append(fig)
				axall.append(axes)
		row = 0
		col = 0

		for flin, fln in enumerate(files):
			fpat, fname, ftype = filename(fln)
			# print plate, fln, row, col
			location = "{}/{}/{}/{}".format(sourceloc, replicate, plate, fln)
			indx = flin + 1
			well = indx2well(flin)
			ttl = '{}-{} | {}'.format(flin + 1, well, nutrients[plate][well]['Metabolite'])
			ttlw = '\n'.join(tw.wrap(ttl, ttllen))

			if indx in alldeletions[replicate][plate].keys():
				print 'Deleting some worms!'
				delworms = alldeletions[replicate][plate][indx]['Worms']
			else:
				delworms = []

			image = tiff.imread(location)
			imgclean = wormdel(image, delworms)

			imghsv = color.rgb2hsv(imgclean)
			imgrgb = color.hsv2rgb(imghsv)
			rhead = [replicate, plate, well, fln]

			if 'Labeling' in figures:
				h = imghsv[:, :, 0]
				v1D = np.around(h.ravel(), 3)
				ftable = np.array(freqtable(v1D))
				ftable = ftable[np.argsort(ftable[:, 0])]
				X, Y, mu, sd = fitgauss(ftable[:, 0], ftable[:, 1])
				# What hue threshold o use

				mthr=mthres[indx] if indx in mthres.keys() else ''

				if replicate in ['Rep1','Rep2','Rep3','Rep4']:
					hthr=mu * 0.995318 + sd * 0.427193 + 0.020434
				else:
					hthr = mu * 0.995318 + sd * 0.427193 + 0.03

				if mthr!='' and mthr!=0:
					hthres = mthr
				else:
					hthres = hthr

				labeled_worms = labeling3(imghsv, hthres,cthres,sthres)


				thresholds.append(rhead+[mu,sd,mthr,hthr])

				wormind = list(np.unique(labeled_worms))
				worms = {w: labeled_worms[labeled_worms == w].shape[0] for w in wormind}
				worms = {wk: wp for wk, wp in worms.items() if wp < 1000000}


				#Make white background
				extract = imgrgb.copy()
				extract[labeled_worms == 0, :] = [1, 1, 1]
				contours = measure.find_contours(labeled_worms, 0.8)

				for w in worms.keys():
					rowhead = [replicate, plate, well, fln, w]
					vw = imghsv[:, :, 2].copy()
					vw[labeled_worms != w] = 0
					bright1D = np.around(vw.ravel(), 3)
					ftable = np.array(freqtable(bright1D))
					ftable = ftable[np.argsort(ftable[:, 0])]
					fdict = {"{0:.3f}".format(freq[0]): int(freq[1]) for freq in ftable}
					values = [fdict[key] if key in fdict.keys() else 0 for key in levels]
					data.append(rowhead + values)


			if printfigures:
				for fi, fign in enumerate(figures):
					fig = figall[fi]
					axes = axall[fi]

					# plt.figure(fi)
					# print row,col
					try:
						plt.sca(axes[row, col])
					except Exception as e:
						print e
						sys.exit(1)
					try:
						ax = axes[row, col]
					except Exception as e:
						print e
						sys.exit(1)
					#ax.axis('off')
					ax.set_adjustable('box-forced')

					if fign == 'Original':
						plt.imshow(imgrgb)
					elif fign == 'Labeling':
						plt.imshow(extract)
						for n, contour in enumerate(contours):
							plt.plot(contour[:, 1], contour[:, 0], linewidth=1)

						for region in regionprops(labeled_worms):
							minr, minc, maxr, maxc = region.bbox
							plt.text(maxc, maxr, region.label)
					else:
						plt.imshow(imgrgb)

					plt.title(ttlw, fontsize=ttlsize)
					plt.setp(ax.get_yticklabels(), visible=False)
					plt.setp(ax.get_xticklabels(), visible=False)

			col = col + 1
			if col == 12:
				col = 0
				row = row + 1

			indo = indo + 1
			prc = (flin + 1) * 100.0 / len(files)
			prco = (indo) * 100.0 / total

			print '{:} {:} {:}:{:6.1f}% | {:6.1f}%'.format(replicate, plate, fname, prc, prco)

		if printfigures:
			for fi, fign in enumerate(figures):
				fig = figall[fi]
				axes = axall[fi]
				# plt.figure(fi)
				fig.tight_layout()

				if fign == 'Original':
					ofignm = '{}/{}_{}_{}.pdf'.format(odir, replicate, plate, fign)
				elif fign == 'Labeling':
					ofignm = '{}/{}_{}_{}_{}.pdf'.format(odir, replicate, plate, label, fign)
				else:
					ofignm = '{}/{}_{}_{}.pdf'.format(odir, replicate, plate, fign)
				fig.savefig(ofignm, bbox_inches='tight')
				plt.close(fig)  # header=['Replicate','Plate','Well','File','Worm','Brightness_value','Frequency']




header = ['Replicate', 'Plate', 'Well', 'File', 'Worm'] + levels
data.insert(0, header)
ofname = '{}/Summary_{}_{}.csv'.format(odir, 'all_new', label)
writecsv(data, ofname, '\t')



header = ['Replicate', 'Plate', 'Well', 'File', 'Worm'] + ['Mu','SD','Manual_t','Estimated_t']
thresholds.insert(0, header)
otname = '{}/Summary_{}_{}.csv'.format(odir, 'thresholds_new', label)
writecsv(thresholds, otname, '\t')


