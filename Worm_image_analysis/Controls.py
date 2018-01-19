
ttllen = 24
ttlsize = 24

minimage = 1
maximage = 96

maxpix = 0


# figures=['Original','Labeling']
#figures = ['Original']

figures = ['Original','Labeling']

#plates = ['PM1']
plates=['PM1','PM2A','PM3B','PM4A']
#replicates = ['Rep1', 'Rep2', 'Rep3', 'Rep4','Rep5', 'Rep6']
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



thresholds_man = readman('Thresholds_Control_manual.csv')



data = []
thresholds=[]


replicates = ['Rep1', 'Rep2', 'Rep3', 'Rep4','Rep5', 'Rep6']





#replicates = ['Rep2']



for replicate in replicates:

	platesel = plates[:]
	if replicate == 'Rep2':
		platesel=['Controls_1930','Controls_2240']
	elif replicate == 'Rep3':
		platesel=['Control_1900']
	elif replicate == 'Rep4':
		platesel = ['Controls_1930']
	elif replicate == 'Rep5':
		platesel = ['NGM_Control']
	else:
		platesel=[]



	for plate in platesel:

		mthres = thresholds_man[replicate[3]][plate]

		#files = ["{}_{}.tif".format(plate, str(im).zfill(3)) for im in range(1, 96 + 1)]
		files = os.listdir("{}/{}/{}".format(sourceloc, replicate, plate))
		files = [fl for fl in files if '.tif' in fl]

		total = len(plates) * len(files) * len(replicates)
		figall = []
		axall = []
		if printfigures:
			for fi, fign in enumerate(figures):
				fig, axes = plt.subplots(nrows=2, ncols=8, figsize=(60, 30), dpi=300)
				fig.suptitle('{}-{}'.format(replicate, plate), fontsize=40)
				plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=1, wspace=1, hspace=0.5)
				figall.append(fig)
				axall.append(axes)
		row = 0
		col = 0

		for flin, fln in enumerate(files):
			fpat, fname, ftype = filename(fln)
			# print plate, fln, row, col
			location = "{}/{}/{}/{}".format(sourceloc, replicate, plate, fln)

			med,tp,indx=fname.split('_')

			print tp, indx

			row=0 if tp=='NoMetf' else 1
			col= int(indx)-1

			#indx = flin + 1
			well = indx2well(flin)
			ttl = '{}'.format(fname)
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

				mthr=mthres[fln] if fln in mthres.keys() else ''

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

			# col = col + 1
			# if col == 8:
			# 	col = 0
			# 	row = row + 1

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
ofname = '{}/Summary_controls_{}_{}.csv'.format(odir, 'All', label)
writecsv(data, ofname, '\t')






header = ['Replicate', 'Plate', 'Well', 'File', 'Worm'] + ['Mu','SD','Manual_t','Estimated_t']
thresholds.insert(0, header)
otname = '{}/Summary_controls_{}_{}.csv'.format(odir, 'thresholds_56', label)
writecsv(thresholds, otname, '\t')






manual = readman('Thresholds_Control.csv')
tmap = manual

sthres = 1000
cthres = 0.02


# Get manual threshold values

nfig = 2
figure, axes = plt.subplots(ncols=nfig, figsize=(8 * nfig, 6), sharex=True,
                            sharey=True)
for fin in [0, 1]:
	ax = axes[fin]
	ax.axis('off')
	ax.set_adjustable('box-forced')

plt.tight_layout()

for rep in tmap.keys():
	# if rep in ['5','6']:
	#	continue
	# a, b, cmin = settings[replicate]
	for pl in tmap[rep].keys():
		mthres = tmap[rep][pl]
		for indx in tmap[rep][pl].keys():

			mthr = mthres[indx] if indx in mthres.keys() else ''

			print 'Manual: {}'.format(mthr)
			if mthr != '' and mthr != 0:
				continue

			print "{} {} {}".format(rep, pl, indx)

			thrs = tmap[rep][pl][indx]

			tfile = "{}/Rep{}/{}/{}".format(sourceloc, rep, pl, indx)

			replicate = 'Rep{}'.format(rep)
			if indx in alldeletions[replicate][plate].keys():
				delworms = alldeletions[replicate][plate][indx]['Worms']
			else:
				delworms = []

			image = tiff.imread(tfile)
			image_clean = wormdel(image, delworms)
			imghsv = color.rgb2hsv(image_clean)
			imgrgb = color.hsv2rgb(imghsv)
			h = imghsv[:, :, 0]
			# s = imghsv[:, :, 1]
			# v = imghsv[:, :, 2]


			v1D = np.around(h.ravel(), 3)
			ftable = np.array(freqtable(v1D))
			ftable = ftable[np.argsort(ftable[:, 0])]
			X, Y, mu, sd = fitgauss(ftable[:, 0], ftable[:, 1])

			if rep in ['1', '2', '3', '4']:
				hthr = mu * 0.995318 + sd * 0.427193 + 0.020434
			else:
				hthr = mu * 0.995318 + sd * 0.427193 + 0.03

			hthres = hthr

			nfig = 2
			figure.suptitle("Rep{} {}-{}".format(rep, pl, indx))
			ax = axes[0]
			ax.imshow(imgrgb)
			ax.set_title("Original")
			while 1:
				print 'Rendering: {} {} {}...'.format(rep, pl, indx)

				# h = imghsv[:, :, 0]


				labeled_worms = labeling3(imghsv, hthres, cthres, sthres)

				# labeled_worms, hmax, hthres = labeling2(imghsv, a, b, cmin, hmin=hmin)
				contours = measure.find_contours(labeled_worms, 0.2)
				extract = imgrgb.copy()
				extract[labeled_worms == 0, :] = [1, 1, 1]
				ax = axes[1]
				ax.clear()
				ax.axis('off')
				ax.set_adjustable('box-forced')
				ax.imshow(extract)
				for n, contour in enumerate(contours):
					ax.plot(contour[:, 1], contour[:, 0], linewidth=1)
				ax.set_title("H_peak={0:.3f} Hthres={1:.3f}".format(mu, hthres))
				time.sleep(0.05)
				plt.pause(0.0001)
				print 'H_peak={0:.3f}, Hthres={1:.3f}'.format(mu, hthres)
				newthr = raw_input("Good threshold/retry?: ")
				if newthr == "":
					tmap[rep][pl][indx] = hthres
					# ofile.writerow([rep,pl,indx,hthres])
					break
				elif newthr == "q":
					raise StopIteration
				hthres = float(newthr)

# csvout.close()








thresholds = map2table(tmap)
header = ['Replicate', 'Plate', 'Index', 'Threshold']

otname = '{}/Thresholds_Control_manual.csv'.format(odir)
thresholds.insert(0, header)
writecsv(thresholds, otname, ',')
