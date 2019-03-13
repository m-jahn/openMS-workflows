from __future__ import division
import pandas, numpy, csv, re
from cStringIO import StringIO
from scipy.stats import pearsonr
import os
import collections

import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from sklearn.neighbors import KNeighborsRegressor
import seaborn as sb
import matplotlib.patches as mpatches #for the plot legend


#consensus_txt is the consensus features map generated from FeatureLinkerUnlabeledQT -> ConsensusMapNormalizer -> output
#consensus_txt = "/home/vitalv/cyano_dataset_20180112/consensusMapNormalizer/consensusMap_co2.csv"
consensus_txt = "/home/vitalv/cyano_dataset_20180112/consensusMapNormalizer/consensusMap_YFP.csv"

#eic_txt is a file with XICs obtained from EICExtractor: 
#XICs are extracted from the RT-corrected mzML files (trafoXML from MapAlignerPose) 
#at the consensus features locations defined from consensusMapNormalizer -> consensus2edta 
#eic_txt = "/home/vitalv/cyano_dataset_20180112/eicExtractor/consensusMap_co2.csv"
eic_txt = "/home/vitalv/cyano_dataset_20180112/eicExtractor/consensusMap_YFP.csv"


num_samples = 5
num_replica = 4



def read_consensus(fn):
	'''
	Read text table from consensusXML exported by OpenMS TextExporter
	'''
	cons_header = []
	pept_header = []
	runs_name = []

	# fetch the headers for consensus features, unassigned peptides and experiments' names. 
	for row in csv.reader(open(fn), delimiter='\t'):
		if row[0] == '#CONSENSUS':
			cons_header = row
		elif row[0] == '#UNASSIGNEDPEPTIDE':
			pept_header = row
		elif row[0] == 'MAP':
			#runs_name.append(row[2])
			runs_name.append(row[2].split("/")[-1].split(".")[0])

	# read consensus features
	s = StringIO()
	with open(fn) as fh:
		for line in fh:
			if line.startswith("CONSENSUS"):
				s.write(line)
	s.seek(0)
	cons = pandas.read_csv(s, sep='\t', header=None, names=cons_header)
	co_peps = []
	with open(fn) as fh:
		for line in fh:
			if line.startswith("CONSENSUS"):
				co_peps.append('')
			elif line.startswith('PEPTIDE') and co_peps[-1] == '':
				# choose the first recorded peptide sequence as consensus sequence
				co_peps[-1] = line.split("\t")[5]
	cons['peptide_0'] = co_peps
	
	# read uassigned peptides as consensus features
	s = StringIO()
	with open(fn) as fh:
		for line in fh:
			if line.startswith("UNASSIGNEDPEPTIDE"):
				s.write(line)
	s.seek(0)
	ua_peps = pandas.read_csv(s, sep='\t', header=None, names=pept_header)
	ua_peps = ua_peps.groupby(['sequence', 'charge']).mean()
	
	return cons, ua_peps, runs_name



def read_eic(fn):
	'''
	Read the detailed output file from EIC extraction in OpenMS.
	Two isotopic peaks (M and M+1) are extracted for each consensus feature and unassigned peptide
	Calculate geometric average of intensities from the two isotopic peaks, as well as diviations of RT and mass.
	'''
	
	with open(fn) as fh:
		eic_str = [StringIO() for _ in range(4)]
		sample_header = [i for i in fh.next().rstrip().split(',') if i] # sample row
		fh.next() # empty row
		cols = fh.next().rstrip().split(',') # quantity headers
		for ix in range(2, len(cols)):
			i = int((ix-2)/5)
			cols[ix] = '_'.join([sample_header[i], cols[ix]]) # rename columns according to sample names
		
		ix = 0
		for line in fh:
			eic_str[ix % 4].write(line)
			ix += 1

		[sio.seek(0) for sio in eic_str]

	# obtain quantities from M and M+1 of target and decoy features, separately.
	eic = pandas.read_csv(eic_str[0], header=None, names=cols) # Monoisotopic    
	eic_iso = pandas.read_csv(eic_str[1], header=None, names=cols) # 13C isotope
	eic_decoy0 = pandas.read_csv(eic_str[2], header=None, names=cols)
	eic_decoy1 = pandas.read_csv(eic_str[3], header=None, names=cols)
	
	for samp in sample_header:
		int_ix = samp + '_intensity'
		rt_ix = samp + '_dRT'
		ppm_ix = samp + '_dppm'

		# values of target features
		eic[samp + '_int_1'] = (eic[int_ix] * eic_iso[int_ix]) ** 0.5
		eic[samp + '_dRT_1'] = eic[rt_ix] - eic_iso[rt_ix]
		eic[samp + '_ppm_1'] = eic[ppm_ix] - eic_iso[ppm_ix]

		# values of decoy features
		eic[samp + '_int_d1'] = (eic_decoy0[int_ix] * eic_decoy1[int_ix]) ** 0.5
		eic[samp + '_dRT_d0'] = eic_decoy0[rt_ix] 
		eic[samp + '_dRT_d1'] = eic_decoy0[rt_ix] - eic_decoy1[rt_ix]
		eic[samp + '_ppm_d0'] = eic_decoy0[ppm_ix]
		eic[samp + '_ppm_d1'] = eic_decoy0[ppm_ix] - eic_decoy1[ppm_ix]
		
	return eic




eic = read_eic(eic_txt)
cons, uapep, ic0 = read_consensus(consensus_txt)


print "Consensus features: %s"% len(cons)
print "Consensus features with an assigned peptide: %s" % len(cons[cons.peptide_0 != ''])
print "Proportion of features with an assigned peptide: %s" % str(round(len(cons[cons.peptide_0 != ''])*100/len(cons),2))
print "Unassigned peptides: %s" % len(uapep)
print "len(eic) %s " % len(eic) + " = len(cons) %s " % len(cons) + " + len(uapep) %s" % len(uapep)

#ic0 = sorted(ic0)
#Note this sorting does not always work the way I want. T.ex: Cyano_1000_R1 will be sorted before Cyano_100_R1, and Cyano_200_ will be sorted before Cyano_60
#It will sort Cyano_1000_RX, Cyano_100_RX, Cyano_200_RX, Cyano_300_RX, Cyano_60_RX. Just change positions between the 60 and the 1000 replicates: 
#ic0 = ic0[-4:]+ic0[4:-4]+ic0[:4]

# replace Intensity Column names
icols = [i for i in cons.columns if i.startswith('intensity_')]
cons = cons[icols + ['quality_cf', 'peptide_0', 'charge_cf', 'rt_cf', 'mz_cf']]
cons.rename(columns=dict(zip(icols[1:], ic0)), inplace=True)

#The YFP consensus includes the 4 1%CO2 replicas, remove them?
#[cons.drop(c, axis=1, inplace=True) for c in cons.columns if "co2" in c]
#ic0 = [c for c in ic0 if "co2" not in c]

# please pardon me for the confusing variable namings here.
# ic1x => feature Intensity Column 1 :XIC (geometric average of M and M+1) 
ic1x = sorted([i for i in eic.columns if i.endswith('int_1')])

# column names of consensus (targets)
#crt0 => rt difference between consensus f and found eic
crt0 = sorted([i for i in eic.columns if i.endswith('dRT')])
#crt1 => rt difference between monoisotopic found eic and C13 found eic
crt1 = sorted([i for i in eic.columns if i.endswith('dRT_1')])
#cppm0 => mz difference (in ppm) between consensus f and found eic
cppm0 = sorted([i for i in eic.columns if i.endswith('ppm')])
#cppm1 => mz difference (in ppm) between monoisotopic found eic and C13 found eic
cppm1 = sorted([i for i in eic.columns if i.endswith('ppm_1')])

# column names of decoys
dic1x = sorted([i for i in eic.columns if i.endswith('int_d1')])
drt0 = sorted([i for i in eic.columns if i.endswith('dRT_d0')])
drt1 = sorted([i for i in eic.columns if i.endswith('dRT_d1')])
dppm0 = sorted([i for i in eic.columns if i.endswith('ppm_d0')])
dppm1 = sorted([i for i in eic.columns if i.endswith('ppm_d1')])



############################################################################################################################################################



# combine Feature and XIC into dataframe dx
dx = pandas.concat([eic, cons], axis=1)
dx['peptide_0'] = [i for i in cons.peptide_0]  + [ i[0] for i in uapep.index.tolist()]
dx['charge_cf'] = [i for i in cons.charge_cf]  + [ i[1] for i in uapep.index.tolist()]
dx['mz_cf'] = [i for i in cons.mz_cf] + uapep.mz.tolist()
dx['rt_cf'] = [i for i in cons.rt_cf] + uapep.rt.tolist()
dx = dx[dx.rt_cf > 0]

dx['peptide'] = [ re.sub('C\(Carbamidomethyl\)', 'C', str(i)) for i in dx.peptide_0]
dx['baseseq'] = [ re.sub('\(.+?\)', '', str(i)) for i in dx.peptide_0]

dx['mods'] = [ sorted(re.findall('\(.+?\)', str(i))) for i in dx.peptide ]
dx['uniq'] = [ "%s%d%s" % (x.baseseq, x.charge_cf, ''.join(x.mods)) for i, x in dx.iterrows()] #iterrows() takes veeery long

# cross-run quantifications from feature-based linking.
dx['f_overlap'] = [ numpy.count_nonzero(numpy.nan_to_num(i)) for i in dx[ic0].values ]
# cross-run quantifications from ion-based linking. 
dx['e_overlap'] = [ numpy.count_nonzero(i) for i in dx[ic1x].values ]
# median intensity in each individual run
dx['medianEIC'] = dx[ic1x].median(axis=1)

print "number of features:\t", len(dx)
print "number of runs:\t", len(ic1x)



# XIC extraction rate for consensus features
e = numpy.array(dx[dx.intensity_cf.isnull() == False].e_overlap).tolist() 
print e.count(len(ic1x)) * 100 / len(e), "% of ", len(e) , " consensus features are extracted across all runs." #(all, identified or not)
#this is important to compare with the proportion of IDENTIFIED UNIQUE peptides across runs plot.
#len(ic1x) : total number of samples where both an ion and its corresponding 13C isotope could be found for each feature 
#counts the  number of times the total number of samples was found in the list of all consensus features

#co2 dataset
#71.2357509307 % of  78163  consensus features are extracted across all runs.

#yfp dataset (+4 CO2 1% replicas, 19 runs)
#67.3573640271 % of  76541  consensus features are extracted across all runs.


############################################################################################################################################################



# Reference set: consensus features with no missing values in both feature linking and XIC extraction
refdx = dx[(dx.f_overlap ==len(ic0)) & (dx.e_overlap == len(ic0))]


fcol = numpy.concatenate([crt0, crt1, cppm0, cppm1])
dcol = numpy.concatenate([drt0, drt1, dppm0, dppm1])

# unit-less transformation. From the paper:
#Since these factors have different units and intervals of changes, are normalized by their own standard deviations and thus converted to 
# unitless quantities that can be simply combined

refz = (refdx[fcol] - refdx[fcol].mean()) / refdx[fcol].std()
decoyz = (refdx[dcol] - refdx[fcol].mean().values) / refdx[fcol].std().values
testz = (dx[fcol] - refdx[fcol].mean()) / refdx[fcol].std()


#Note, if ic1x columns were not sorted (when creating them from eic columns ending with int_1) sort them now. otherwise samples of the same group (replicas) will not be reshaped (grouped) together!:
ic1x = sorted(ic1x)

#for YFP dataset (+4 CO2 1% replicas, 19 runs):
ic1x_reshape = \
[i for i in numpy.array([i for i in dic1x if "300_ja8" in i]).reshape(1,3)] + \
[i for i in numpy.array([i for i in dic1x if not "300_ja8" in i]).reshape(4,4)]

#for i in numpy.array(ic1x).reshape(num_samples, num_replica):
for i in ic1x_reshape: #for YFP dataset
	cv = refdx[i].std(axis=1) / refdx[i].mean(axis=1)
	cv = numpy.sqrt(cv)
	for run_id in i:
		refz[run_id + '_cv'] = cv


#dic1x = sorted(dic1x)

#for YFP dataset (+4 CO2 1% replicas, 19 runs)
dic1x_reshape = \
[i for i in numpy.array([i for i in dic1x if "300_ja8" in i]).reshape(1,3)] + \
[i for i in numpy.array([i for i in dic1x if not "300_ja8" in i]).reshape(4,4)]

#for i in numpy.array(dic1x).reshape(num_samples, num_replica):
for i in dic1x_reshape: #YFP dataset misses 1 repl
	cv = refdx[i].std(axis=1) / refdx[i].mean(axis=1)
	cv = numpy.sqrt(cv)
	for run_id in i:
		decoyz[run_id + '_cv'] = cv

#for i in numpy.array(ic1x).reshape(num_samples, num_replica):
for i in ic1x_reshape:
	cv = dx[i].std(axis=1) / dx[i].mean(axis=1)
	cv = numpy.sqrt(cv)
	for run_id in i:
		testz[run_id + '_cv'] = cv
		
ccv = [i for i in refz.columns if i.endswith('cv')]
dcv = [i for i in decoyz.columns if i.endswith('cv')]


fil_cols = zip(ic1x, crt0, crt1, cppm0, cppm1, ccv)
fil_cols_decoy = zip(dic1x, drt0, drt1, dppm0, dppm1, dcv)


zscores = []
for f in fil_cols:
	e = refdx[f[0]] > 0
	rz = refz[e]
	z = [numpy.sum(v*v) for v in rz[list(f[1:])].values]
	zscores = zscores + z

dzscores = []
for f in fil_cols_decoy:
	e = refdx[f[0]] > 0
	dd = decoyz[e]
	z = [numpy.sum(v*v) for v in dd[list(f[1:])].values]
	dzscores = dzscores + z


#####PLOT FDR ION-PIP Score ##############
bins = numpy.arange(-6.6, 4.5, 0.1)
pandas.Series(-numpy.log(zscores)).hist(bins=bins, alpha = 0.5, color ='b', lw=0, label='Target')
pandas.Series(-numpy.log(dzscores)).hist(bins=bins, alpha = 0.5, color='r', lw=0, label='Decoy')
# pandas.Series(-numpy.log(tzscores)).hist(bins=bins, alpha = 0.3, color='y', lw=0, label='All')
plt.legend()
plt.xlabel('Score')
plt.ylabel('# Features')
plt.show()
##########################################




score_list = sorted(zip(zscores + dzscores, [0 for i in zscores] + [1 for i in zscores]), key=lambda x: x[0])
score_cutoff = 0
hit_count = 0
decoy_count = 0

# set FDR threshold 
fdr = 0.05
for s in score_list:
	hit_count += 1
	score_cutoff = -numpy.log(s[0])
	if s[1] > 0:
		decoy_count += 1
	if decoy_count / hit_count >= fdr:
		print "cutoff score:", score_cutoff
		break




#FILTER: Set to 0 XICs for those consensus that have scores lower than score_cutoff:

#fil_cols = zip(ic1x, crt0, crt1, cppm0, cppm1, ccv)
# ic1x => EIC Intensity Column 1 :XIC (geometric average of M and M+1) 
for f in fil_cols:
	score = numpy.array([numpy.sum(v*v) for v in testz[list(f[1:])].values])
	score = -numpy.log(score)
	print list(score > score_cutoff).count(True) / len(score)
	dx[f[0]][list(score <= score_cutoff)] = 0 # XIC filter: "The features that failed to pass the threshold were considered missing (zero-intensity)"


	
dx['medianEIC'] = dx[dx[ic1x] > 0][ic1x].median(axis=1)
dx['e_overlap'] = [ numpy.count_nonzero(i) for i in dx[ic1x].values ]


print "Features before filtering:\t", len(dx)




#"After filling missing values, the consensus map was further filtered by removing features that failed to be quantified in at least one replicate in each sample"
# quantified at least in one run in each sample
#ix = numpy.min([numpy.mean(dx[i].values, axis=1) for i in numpy.array(ic1x).reshape(num_samples, num_replica)], axis=0) > 0 
ix = numpy.min([numpy.mean(dx[i].values, axis=1) for i in ic1x_reshape], axis=0) > 0 

dx = dx[ix]


print "Features after filtering:\t", len(dx)
print "Unique peptide sequences:\t", len(dx.baseseq.unique())










icols = zip(sorted(ic0), sorted(ic1x))
for a in icols[:]:
	do = dx[pandas.notnull(dx[a[0]])]
	do = do[do[a[1]] > 0]
	xy = do[list(a)].apply(numpy.log2).values
	r2 = pearsonr(xy[:,0], xy[:, 1])[0] ** 2
	print a[1], len(xy), r2
	
	if r2 < 0.5:
		# Remove problematic runs with low correlation between feature abundance and XIC intensity.
		icols.remove(a)
		print a, " has been removed due to low correlation."

# 20180112_co2-0p15_r1_GC2_01_5489.mzML_int_1 6377 0.909839953335
# 20180112_co2-0p15_r2_GC3_01_5490.mzML_int_1 6586 0.908063933816
# 20180112_co2-0p15_r3_GC4_01_5491.mzML_int_1 6205 0.907578554596
# 20180112_co2-0p15_r4_GC5_01_5492.mzML_int_1 6448 0.889573870027
# 20180112_co2-0p20_r1_GB6_01_5484.mzML_int_1 6517 0.910221269159
# 20180112_co2-0p20_r2_GB7_01_5485.mzML_int_1 6627 0.914645232267
# 20180112_co2-0p20_r3_GB8_01_5486.mzML_int_1 6779 0.918745155007
# 20180112_co2-0p20_r4_GC1_01_5487.mzML_int_1 6490 0.923326631875
# 20180112_co2-0p30_r1_GB2_01_5479.mzML_int_1 6911 0.935146576714
# 20180112_co2-0p30_r2_GB3_01_5480.mzML_int_1 6869 0.919658603912
# 20180112_co2-0p30_r3_GB4_01_5481.mzML_int_1 6545 0.911873260304
# 20180112_co2-0p30_r4_GB5_01_5482.mzML_int_1 6970 0.911969045027
# 20180112_co2-0p50_r1_GA6_01_5474.mzML_int_1 6712 0.936189248628
# 20180112_co2-0p50_r2_GA7_01_5475.mzML_int_1 6780 0.927628713612
# 20180112_co2-0p50_r3_GA8_01_5476.mzML_int_1 6407 0.934017471431
# 20180112_co2-0p50_r4_GB1_01_5477.mzML_int_1 6793 0.923298374275
# 20180112_co2-1p00_r1_GA2_01_5469.mzML_int_1 6221 0.93340142286
# 20180112_co2-1p00_r2_GA3_01_5470.mzML_int_1 6059 0.923595020186
# 20180112_co2-1p00_r3_GA4_01_5471.mzML_int_1 6000 0.920571656523
# 20180112_co2-1p00_r4_GA5_01_5472.mzML_int_1 6000 0.918828785897

# 20180112_300_ja2_r1_GC6_01_5494.mzML_int_1 5727 0.835055452615
# 20180112_300_ja2_r2_GC7_01_5495.mzML_int_1 5775 0.84736047898
# 20180112_300_ja2_r3_GC8_01_5496.mzML_int_1 5630 0.841888429302
# 20180112_300_ja2_r4_GD1_01_5497.mzML_int_1 5638 0.836229390161
# 20180112_300_ja8_r2_GD3_01_5500.mzML_int_1 6574 0.725071993082
# 20180112_300_ja8_r3_GD4_01_5501.mzML_int_1 5511 0.851943966971
# 20180112_300_ja8_r4_GD5_01_5502.mzML_int_1 5411 0.824187576703
# 20180112_60_ja2_r1_GD6_01_5504.mzML_int_1 6201 0.839603223823
# 20180112_60_ja2_r2_GD7_01_5505.mzML_int_1 6491 0.797720854365
# 20180112_60_ja2_r3_GD8_01_5506.mzML_int_1 6415 0.79768483836
# 20180112_60_ja2_r4_GE1_01_5507.mzML_int_1 6645 0.788358444378
# 20180112_60_ja8_r1_GE2_01_5509.mzML_int_1 6425 0.833516422644
# 20180112_60_ja8_r2_GE3_01_5510.mzML_int_1 6347 0.833226157013
# 20180112_60_ja8_r3_GE4_01_5511.mzML_int_1 6428 0.812601501637
# 20180112_60_ja8_r4_GE5_01_5512.mzML_int_1 6279 0.831171250139
# 20180112_co2-1p00_r1_GA2_01_5469.mzML_int_1 6291 0.840450867054
# 20180112_co2-1p00_r2_GA3_01_5470.mzML_int_1 5757 0.84994171058
# 20180112_co2-1p00_r3_GA4_01_5471.mzML_int_1 5477 0.840896977678
# 20180112_co2-1p00_r4_GA5_01_5472.mzML_int_1 5853 0.835984804572




############################################################################################################################################################


# log2 transform of selected intensity colums from the dx
data = dx[ic0 + ic1x].apply(numpy.log2)


# impute remaining missing values as the minimum detectable intensity (int_1):
# NOTE: This imputation is for -inf, which are the ic1x (int_1) cols (ion extracted intensities), NOT consensus(ic0) features intensities!

data[data==-numpy.inf] = data[data!=-numpy.inf][ic1x].min().min()


# NOTE: Don't do imputation YFP peptides for the 1% co2 replicas that go with the YFP dataset!
from Bio import SeqIO
yfp_seqs = SeqIO.parse(open('/home/vitalv/database/pJA_YFP_vectors.fasta'),'fasta')
yfp_seqs_str = ""
for y in yfp_seqs:
	yfp_seqs_str += y.seq.tostring()

data['baseseq'] = dx.baseseq
ixs = []
for i, r in data.iterrows():
	if r.baseseq in yfp_seqs_str:
		ixs.append(i)
		print r.baseseq
		data.loc[i, ic1x] = [numpy.log2(v) for v in dx.loc[i, ic1x].values]



# or impute missing values as zero intensity
#data[data==-numpy.inf] = 0

data['rt'] = dx['rt_cf']
data['mass'] = dx['mz_cf'] * dx['charge_cf']
data['charge'] = dx['charge_cf']
data['peptide'] = dx.peptide
data['baseseq'] = dx.baseseq
data['uniq'] = dx.uniq 

for a, b in icols:
	dn = data[(data[a].notnull()) & (data[b] > 0)] #If impute int_1 with min intensity: all data[b] will be > 0 since I just imputed these (0 in dx, -inf in data, after log2) 
	mx = numpy.matrix(dn[[a,b,'mass']])
	# KNN regression (k=5 by default) for predicting feature abundance (log2) based on ion intensity and precosor mass
	regr = KNeighborsRegressor().fit(mx[:,1:], mx[:,0])
	a_ = regr.predict(numpy.matrix(data[[b, 'mass']]))[:,0]
	a_[data[b].values==0] = 0 # keep missing values from extraction
	data[a + '_'] = a_

ic2 = [a+'_' for a in ic0]
tmp = dx[ic0].apply(numpy.log2)
tmp.values[numpy.isnan(tmp.values)] = data[ic2].values[numpy.isnan(tmp.values)]
data[ic0] = tmp.apply(numpy.exp2) # converting back to linear scale



print data[ic0].values.flatten().tolist().count(1) * 1. / data[ic0].values.flatten().__len__()







'''
# quantified vs. identified #
#print len(data[data.peptide != '']) *100./ data[ic0].__len__(), "% features have been assigned with peptide sequence."
print len(data[data.peptide.notnull()]) *100./ data[ic0].__len__(), "% features have been assigned with peptide sequence."
data[ic0].apply(numpy.log10).median(axis=1).hist(bins=numpy.arange(4,10,0.1), lw=0, alpha=0.9)
#data[data.peptide != ''][ic0].apply(numpy.log10).median(axis=1).hist(bins=numpy.arange(4,10,0.1), lw=0, alpha=0.9)
data[data.peptide.notnull()][ic0].apply(numpy.log10).median(axis=1).hist(bins=numpy.arange(4,10,0.1), lw=0, alpha=0.9)
plt.xlabel("Abundance (log10)")
plt.ylabel("Frequency")
plt.title("Quantified VS Identified")
'''




pepdata = data[data.peptide.notnull()]

# remove duplications
pepdata.drop_duplicates(subset=['uniq'], inplace=True)

# remove modified peptides?
# pepdata = pepdata[[len(i) == 0 for i in pepdata.mods]]


# report matrix
mx = pepdata.groupby(['baseseq'])[ic0].sum()



# save peptide-level quantification result to a CSV file
#mx.to_csv(r'peptide_quant_co2.csv')
mx.to_csv(r'peptide_quant_yfp.csv')



#/usr/bin/python3.6 /home/vitalv/diffacto_old/src/run_diffacto.py -peptide='/home/vitalv/cyano_dataset_20180112/peptide_quant_co2.csv' -db='/home/vitalv/database/Synechocystis_PCC6803_protein_sequences_YFP.fasta' -samples='/home/vitalv/cyano_dataset_20180112/diffacto/samples.csv' -ref=False -out='/home/vitalv/cyano_dataset_20180112/diffacto/diffacto_out_peptides.DeMixQ.tsv' -out_type='peptide'

#/usr/bin/python3.6 /home/vitalv/diffacto_old/src/run_diffacto.py -peptide='/home/vitalv/cyano_dataset_20180112/peptide_quant_yfp.csv' -db='/home/vitalv/database/Synechocystis_PCC6803_protein_sequences_YFP.fasta' -samples='/home/vitalv/cyano_dataset_20180112/diffacto/samples_YFP.csv' -ref=False -out='/home/vitalv/cyano_dataset_20180112/diffacto/diffacto_out_peptides_YFP.tsv' -out_type='peptide'