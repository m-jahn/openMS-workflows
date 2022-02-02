#!/usr/bin/python

import os 
import process_handler
import argparse
import sys
from multiprocessing import Pool


'''
Author: Vital Vialas (vital.vialas@scilifelab.se)


LABEL FREE QUANTIFICATION WORKFLOW:

Runs OpenMS tools for Label Free Quantification on LC-MS/MS data
The key steps are:

-FeatureFinderIdentification. Finds MS1 features based on positions determined from MS2 (high confidence identification) spectra
	ftp://ftp.mi.fu-berlin.de/pub/OpenMS/release-documentation/html/TOPP_FeatureFinderIdentification.html

-MapAlignerPoseClustering. Map alignment between runs. Correct retention time distortions
	http://ftp.mi.fu-berlin.de/pub/OpenMS/release-documentation/html/TOPP_MapAlignerPoseClustering.html

-FeatureLinkerUnlabeledQT. Creates a consensus map of equivalent, aligned features
	http://ftp.mi.fu-berlin.de/pub/OpenMS/release-documentation/html/TOPP_FeatureLinkerUnlabeledQT.html

-EICextractor. Obtains XICs at the m/z and RT coordinates of the consensus features
	http://ftp.mi.fu-berlin.de/pub/OpenMS/release1.9-documentation/html/TOPP_EICExtractor.html

-DeMixQ processing. Fills in missing feature abundance values estimated from a KNN (k=5) averaging abundances from other quantified features having the most similar *FDR-filtered* XIC patterns in the same run
	


DEPENDENCIES:

OpenMS >= 2.0 
Python >= 2.7
-numpy
-pandas
-lxml
-pymzml http://pymzml.github.io/
-pyteomics http://pythonhosted.org//pyteomics/


REQUIRES:

-Folder with mzML files 
-Folder with search identification results in idXML (OpenMS format). One per mzML file
	(identification result files in idXML format can be obtained running msgfplus.py and qvality.py)

- DeMixQ consensus2edta.py: https://github.com/userbz/DeMix-Q/blob/master/scripts/consensus2edta.py 
	*will be downloaded to current (__file__) folder if not present 
- DeMixQ post-processing script: https://github.com/userbz/DeMix-Q/blob/master/scripts/DeMixQ_data_processing.py
	*will be downloaded to current (__file__) folder if not present 

- process_handler.py, bundled together with this script. Makes it possible to run locally or in a cluster with the Slurm queueing system


RUN:

python lfq_workflow.py 
	-project_folder=/absolute/or/relative/path/to/folder/where/results/will/be/created 
	-mzML_folder=/absolute/or/relative/path/to/mzML/ 
	-id_results_folder=/absolute/or/relative/path/to/id_results/ 
	-n_samples 5 
	-n_replica 4

'''    


print("\n--lfq_worflow--\n")
ph = process_handler.ProcessHandler()()
print("ProcessHandler object created on pipeline.py\n")


#Read arguments , basic checkings and create folders for results 
#-------------------------------------------------------------------------------------------------------------------------------------

apars = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
apars.add_argument('-project_folder', help = 'folder where all results will be created', required=True)
apars.add_argument('-mzML_folder', default= 'mzML/', help = 'dataset folder with mzML files', required=True)
apars.add_argument('-id_results_folder', help='folder containing FDR filtered idXML ', required=True)
apars.add_argument('-n_samples',  help = 'number of samples', required=True)
apars.add_argument('-n_replica', help='number of replicate runs per sample ', required=True)
args = apars.parse_args()

project_folder = args.project_folder
mzML_folder = args.mzML_folder
id_results_folder =  args.id_results_folder
n_samples = args.n_samples
n_replica = args.n_replica

if not os.path.isdir(mzML_folder):
	sys.exit('\nmzML_folder not found\n')
if not os.path.isdir(id_results_folder):
	sys.exit('\nid_results_folder not found\n')

#Create paths for storing results of each of the steps in the workflow:
if not os.path.isdir(project_folder): os.mkdir(project_folder)

feature_finder_folder = os.path.join(project_folder, "featureFinderId/")
id_mapper_folder = os.path.join(project_folder, "idMapper/")
map_aligner_folder = os.path.join(project_folder, "mapAlignerPose/")
feature_linker_folder = os.path.join(project_folder, "featureLinker/")
consensus_map_normalizer_folder = os.path.join(project_folder, "consensusMapNormalizer/")
text_exporter_folder = os.path.join(project_folder, "textExporter/")
consensus2edta_folder = os.path.join(project_folder, "consensus2edta/")
map_rt_transformer_folder = os.path.join(project_folder, "mapRTTransformer/")
eic_extractor_folder = os.path.join(project_folder, "eicExtractor/")
demixq_processing_folder = os.path.join(project_folder, "DeMixQ")



#-------------------------------------------------------------------------------------------------------------------------------------
#FeatureFinderIdentification

def batchFeatureFinderIdentification(dataset_dir, id_results_folder, feature_finder_folder=feature_finder_folder):
	if not os.path.isdir(feature_finder_folder): os.mkdir(feature_finder_folder)
	pool = ph.getPool()
	for mzML_file in sorted(os.listdir(dataset_dir)):
		id_file = os.path.join(id_results_folder, mzML_file.replace(".mzML", ".idXML"))
		pool.applyAsync(processFeatureFinderId, [os.path.join(dataset_dir, mzML_file), id_file, feature_finder_folder])
	results = pool.checkPool()
	featureFinder_job_ids = [r.get() for r in results]
	return featureFinder_job_ids

def processFeatureFinderId(mzML_file, id_file, feature_finder_folder):
	#demixq paper parameters:
	#mz_window = 10.0 # default 10
	#n_isotopes = 2 # default 2
	#peak_width = 60.0 # default 60 
	#mapping_tolerance = 0.0 default 0
	threads = 10
	featureFinderCmd = "FeatureFinderIdentification " +\
									 " -in %s"%  mzML_file + \
									 " -id %s"%  id_file + \
									 " -threads %s"%threads +\
									 " -out %s"%feature_finder_folder+mzML_file.split("/")[-1].replace(".mzML", ".featureXML")
	featureFinder_job_id = ph.executeCmd(featureFinderCmd, jobType="featureFinder")
	return featureFinder_job_id




#-------------------------------------------------------------------------------------------------------------------------------------
#IdMapper

def batchIdMapper(dependent_on_job_ids, dataset_dir, id_results_folder, feature_finder_folder=feature_finder_folder, id_mapper_folder=id_mapper_folder):
	if not os.path.isdir(id_mapper_folder):	os.mkdir(id_mapper_folder)
	pool = ph.getPool()
	for idXML_file in sorted([i for i in os.listdir(id_results_folder) if i.endswith(".idXML")]):
		featureXML_file = os.path.join(feature_finder_folder, idXML_file.replace(".idXML", ".featureXML"))
		idXML_file = os.path.join(id_results_folder,idXML_file)		
		pool.applyAsync(processIdMapper, [idXML_file, featureXML_file, id_mapper_folder, dependent_on_job_ids])
	results = pool.checkPool()
	idMapper_job_ids = [r.get() for r in results]
	return idMapper_job_ids

def processIdMapper(idXML_file, featureXML_file, id_mapper_folder, dependent_on_job_ids):
	#demixq paper parameters rt_tol = 15, mz_tol = 5
	rt_tol = 30
	mz_tol = 10
	threads = 10
	out_featureXML_file = os.path.join(id_mapper_folder, os.path.basename(featureXML_file))
	idMapperCmd = "IDMapper -id %s"% idXML_file + \
								" -rt_tolerance %s" % rt_tol + \
								" -mz_tolerance %s" % mz_tol + \
								" -in %s"% featureXML_file + \
								" -threads %s"%threads + \
								" -out %s"% out_featureXML_file
	idMapper_job_id = ph.executeCmd(idMapperCmd, jobType="idMapper", dependent_on_job_ids=dependent_on_job_ids)
	return idMapper_job_id




#-------------------------------------------------------------------------------------------------------------------------------------
#MapAlignerPoseClustering


def batchMapAlignerPose(dependent_on_job_ids, dataset_dir, id_mapper_folder=id_mapper_folder, map_aligner_folder=map_aligner_folder):

	if not os.path.isdir(map_aligner_folder): os.mkdir(map_aligner_folder)
	pool = ph.getPool()
	#remember: id_files_folder might be empty but I still want to processMapAligner and executeCmd so I can make it dependent on id_mapper_job_ids
	#but luckily the names of these idXML files is the name of the mzML files in dataset dir but located under id_mapper_folder and ending in .idXML
	in_list = sorted([os.path.join(id_mapper_folder, mzML_file.replace(".mzML", ".featureXML")) for mzML_file in os.listdir(dataset_dir) if mzML_file.endswith('.mzML')])
	out_list = sorted([os.path.join(map_aligner_folder, os.path.basename(idXML_file)) for idXML_file in in_list])
	trafo_out_list = sorted([idXML_file.replace(".featureXML", ".trafoXML") for idXML_file in out_list ])
	pool.applyAsync(processMapAlignerPose, [dependent_on_job_ids, in_list, out_list, trafo_out_list])
	results = pool.checkPool()
	mapAligner_job_ids = [r.get() for r in results]
	return mapAligner_job_ids 

def processMapAlignerPose(dependent_on_job_ids, in_list, out_list, trafo_out_list):
	threads = 10
	#untouched default values from the OpenMS MapAlignerPoseClustering
	#http://ftp.mi.fu-berlin.de/pub/OpenMS/release-documentation/html/TOPP_MapAlignerPoseClustering.html
	mapAlignerCmd = "MapAlignerPoseClustering -in %s"%' '.join(in_list) + \
								" -threads %s" % threads + \
								" -trafo_out %s"%' '.join(trafo_out_list) + \
								" -out %s"%' '.join(out_list)
	mapAligner_job_id = ph.executeCmd(mapAlignerCmd, jobType="mapAligner", dependent_on_job_ids=dependent_on_job_ids)
	return mapAligner_job_id




#-------------------------------------------------------------------------------------------------------------------------------------
#FeatureLinkerUnlabeledQT

def batchFeatureLinkerUnlabeledQT(dependent_on_job_ids, dataset_dir, map_aligner_folder=map_aligner_folder, feature_linker_folder=feature_linker_folder):
	if not os.path.isdir(feature_linker_folder): os.mkdir(feature_linker_folder)
	pool = ph.getPool()
	in_list = sorted([os.path.join(map_aligner_folder, mzML_file.replace(".mzML", ".featureXML")) for mzML_file in os.listdir(dataset_dir)if mzML_file.endswith('.mzML')])
	pool.applyAsync(processFeatureLinkerUnlabeledQT, [in_list, feature_linker_folder, dependent_on_job_ids])
	results = pool.checkPool()
	featureLinker_job_id  = [r.get() for r in results]
	return featureLinker_job_id

def processFeatureLinkerUnlabeledQT(in_list, feature_linker_folder, dependent_on_job_ids):
	rt_diff = 100 #default(paper): 180 (s)
	mass_diff = 0.3 #default(paper): 5 (ppm)
	threads = 10
	featureLinkerCmd = "FeatureLinkerUnlabeledQT -in %s"%' '.join(in_list) +\
										 " -algorithm:distance_RT:max_difference %s"% rt_diff +\
										 " -algorithm:distance_MZ:max_difference %s"% mass_diff +\
										 " -algorithm:distance_MZ:unit Da" +\
										 " -threads %s"% threads +\
										 " -out %s"%os.path.join(feature_linker_folder,"consensusMap.consensusXML") 
	featureLinker_job_id = ph.executeCmd(featureLinkerCmd, jobType="featureLinker", dependent_on_job_ids=dependent_on_job_ids)
	return featureLinker_job_id




#-------------------------------------------------------------------------------------------------------------------------------------
#ConsensusMapNormalizer

def batchConsensusMapNormalizer(dependent_on_job_ids, dataset_dir, consensus_map_normalizer_folder=consensus_map_normalizer_folder):
	if not os.path.isdir(consensus_map_normalizer_folder): os.mkdir(consensus_map_normalizer_folder)
	pool = ph.getPool()
	in_cons = os.path.join(feature_linker_folder,"consensusMap.consensusXML")
	pool.applyAsync(processConsensusMapNormalizer, [in_cons, consensus_map_normalizer_folder, dependent_on_job_ids])
	results = pool.checkPool()
	consensusMapNormalizer_job_id  = [r.get() for r in results]
	return consensusMapNormalizer_job_id

def processConsensusMapNormalizer(in_cons, consensus_map_normalizer_folder, dependent_on_job_ids):
	algorithm = 'median'
	threads = 10
	consensusMapNormalizerCmd = "ConsensusMapNormalizer -in %s" % in_cons + \
								 " -algorithm_type %s" % algorithm + \
								 " -out  %s " % os.path.join(consensus_map_normalizer_folder, "consensusMap.consensusXML")
	consensusMapNormalizer_job_id = ph.executeCmd(consensusMapNormalizerCmd, jobType="consensusMapNormalizer", dependent_on_job_ids=dependent_on_job_ids)
	return consensusMapNormalizer_job_id




#-------------------------------------------------------------------------------------------------------------------------------------
#TextExporter

def textExporter(dependent_on_job_ids):
	if not os.path.isdir(text_exporter_folder):	os.mkdir(text_exporter_folder)
	consensusXML_file = os.path.join(consensus_map_normalizer_folder, "consensusMap.consensusXML")
	textExporterCmd = "TextExporter -in %s"% consensusXML_file + \
										" -out %s"% os.path.join(text_exporter_folder, os.path.basename(consensusXML_file.replace(".consensusXML", ".csv")))
	textExporter_job_id = ph.executeCmd(textExporterCmd, jobType="textExporter", dependent_on_job_ids=dependent_on_job_ids)
	return textExporter_job_id




#-------------------------------------------------------------------------------------------------------------------------------------
#consensus2edta.py

def consensus2edta(dependent_on_job_ids):
	if not os.path.isdir(consensus2edta_folder): os.mkdir(consensus2edta_folder)
	python_script =  "consensus2edta.py" #current folder. same as __file__
	consensus_csv = os.path.join(text_exporter_folder, "consensusMap.csv")
	decoy_rt = 300
	decoy_ppm = 50
	out_dir = consensus2edta_folder
	consensus2edtaCmd = "python %s %s %d %d %s" % (python_script, consensus_csv, decoy_rt, decoy_ppm, out_dir)
	consensus2edta_job_id = ph.executeCmd(consensus2edtaCmd, jobType="consensus2edta", dependent_on_job_ids=dependent_on_job_ids)
	return consensus2edta_job_id
 




#-------------------------------------------------------------------------------------------------------------------------------------
#MapRTTransformer

def batchMapRTTransformer(dependent_on_job_ids, dataset_dir, map_aligner_folder=map_aligner_folder, map_rt_transformer_folder=map_rt_transformer_folder):
	if not os.path.isdir(map_rt_transformer_folder): os.mkdir(map_rt_transformer_folder)
	pool = ph.getPool()
	for mzML_file in sorted([os.path.join(dataset_dir, mzML_file) for mzML_file in os.listdir(dataset_dir) if mzML_file.endswith('.mzML')] ):
		trafoXML_file =  os.path.join(map_aligner_folder, os.path.basename(mzML_file).replace(".mzML", ".trafoXML"))
		pool.applyAsync(processMapRTTransformer, [mzML_file, trafoXML_file, map_aligner_folder, map_rt_transformer_folder, dependent_on_job_ids])
	pool.checkPool()
	results = pool.checkPool()
	mapRTTransformer_job_ids = [r.get() for r in results]
	return mapRTTransformer_job_ids 

def processMapRTTransformer(mzML_file, trafoXML_file, map_aligner_folder, map_rt_transformer_folder, dependent_on_job_ids):
	threads = 10
	out_mzML_file = os.path.join(map_rt_transformer_folder, os.path.basename(mzML_file))
	mapRTTransformerCmd = "MapRTTransformer -in %s"% mzML_file + \
												" -trafo_in %s"% trafoXML_file + \
												" -model:type b_spline" + \
												" -model:linear:symmetric_regression " + \
												" -threads %s" % threads +\
												" -out %s" % out_mzML_file
	mapRTTransformer_job_id = ph.executeCmd(mapRTTransformerCmd, jobType="mapRTTransformer", dependent_on_job_ids=dependent_on_job_ids)
	return mapRTTransformer_job_id





#-------------------------------------------------------------------------------------------------------------------------------------
#EICExtractor

def batchEICExtractor(dependent_on_job_ids, dataset_dir, map_rt_transformer_folder=map_rt_transformer_folder, consensus2edta_folder=consensus2edta_folder, eic_extractor_folder=eic_extractor_folder):
	if not os.path.isdir(eic_extractor_folder):	os.mkdir(eic_extractor_folder)
	edta_file = consensus2edta_folder + "consensusMap.edta"
	mzML_file_list = sorted([os.path.join(map_rt_transformer_folder, mzML_file) for mzML_file in os.listdir(dataset_dir) if mzML_file.endswith(".mzML")])
	processEICExtractorlist(dependent_on_job_ids, edta_file, mzML_file_list, eic_extractor_folder)

def processEICExtractorlist(dependent_on_job_ids, edta_file, mzML_file_list, eic_extractor_folder):
	rt_tol = 120 #DeMixQ default
	mz_tol = 10 #DeMixQ default
	threads = 10
	eicExtractorCmd = "EICExtractor -in %s "% ' '.join(mzML_file_list) + \
										" -pos %s "% edta_file + \
										" -rt_tol %d " % rt_tol + \
										" -mz_tol %d " % mz_tol + \
										" -threads %d " % threads + \
										" -out %s " % os.path.join(eic_extractor_folder, edta_file.split("/")[-1].replace(".edta", ".csv") )
	eicExtractor_job_id = ph.executeCmd(eicExtractorCmd, jobType="eicExtractor", dependent_on_job_ids=dependent_on_job_ids)
	return eicExtractor_job_id




#-------------------------------------------------------------------------------------------------------------------------------------
#DeMixQ_processing

def deMixQ_processing(dependent_on_job_ids, n_samples, n_replica):

	#DeMixQ_data_processing
	if not os.path.isdir(demixq_processing_folder): os.mkdir(demixq_processing_folder)
	demixq_processing_script = 'DeMixQ_data_processing.py' #current folder. same as __file__
	eic_file = os.path.join(eic_extractor_folder, "consensusMap.csv")
	cons_file = os.path.join(text_exporter_folder, "consensusMap.csv")
	if os.path.exists(eic_file) and os.path.exists(cons_file):
		demixq_processing_cmd = "python %s " % demixq_processing_script +\
								" -eic %s " % eic_file +\
								" -cons %s " % cons_file +\
								" -n_samples %s " % n_samples +\
								" -n_repeats %s " % n_replica +\
								" -knn_k 0 " +\
								" -out %s  " % os.path.join(demixq_processing_folder, 'peptides.DeMixQ.csv') 
		demixq_processing_job_id = ph.executeCmd(demixq_processing_cmd, jobType="demixq_processing", dependent_on_job_ids=dependent_on_job_ids)
	else:
		print("make sure consensusMap.csv exists in %s and consensusMap.csv in %s " %(eic_extractor_folder, text_exporter_folder))








def main():

	featureFinder_job_ids = batchFeatureFinderIdentification(dataset_dir=mzML_folder, id_results_folder=id_results_folder, feature_finder_folder=feature_finder_folder)
	idMapper_job_ids = batchIdMapper(dataset_dir=mzML_folder, id_results_folder=id_results_folder, dependent_on_job_ids=featureFinder_job_ids)
	mapAlignerPose_job_ids = batchMapAlignerPose(dataset_dir=mzML_folder, dependent_on_job_ids=idMapper_job_ids)
	featureLinker_job_ids = batchFeatureLinkerUnlabeledQT(dataset_dir=mzML_folder, dependent_on_job_ids=mapAlignerPose_job_ids)
	consensusMapNormalizer_job_ids = batchConsensusMapNormalizer(dataset_dir=mzML_folder, dependent_on_job_ids=featureLinker_job_ids)
	textExporter_job_id = textExporter(dependent_on_job_ids=consensusMapNormalizer_job_ids)
	if not os.path.exists('consensus2edta.py'):
		os.system('wget https://raw.githubusercontent.com/userbz/DeMix-Q/master/scripts/consensus2edta.py')
	consensus2edta_job_id = consensus2edta(dependent_on_job_ids=textExporter_job_id)
	mapRTTransformer_job_ids = batchMapRTTransformer(dataset_dir=mzML_folder, dependent_on_job_ids=mapAlignerPose_job_ids)
	eicExtractor_job_id = batchEICExtractor(dataset_dir=mzML_folder, dependent_on_job_ids=mapRTTransformer_job_ids+[consensus2edta_job_id])

	print("\n finished LFQ workflow. Starting DeMixQ_data_processing\n")
	if not os.path.exists('DeMixQ_data_processing.py'):
		print("DeMixQ_data_processing.py not found in current dir. Downloading from https://raw.githubusercontent.com/userbz/DeMix-Q/master/scripts/DeMixQ_data_processing.py")
		os.system('wget https://raw.githubusercontent.com/userbz/DeMix-Q/master/scripts/DeMixQ_data_processing.py')
	#remember I have made some changes to DeMixQ_data_processing.py (ic0 -column headers- is base name not whole path, and missing value imputation is min, not 0)
	deMixQ_processing(dependent_on_job_ids=[consensus2edta_job_id]+[eicExtractor_job_id], n_samples=n_samples, n_replica=n_replica)

	'''
	#Next step would be to run diffacto on the generated peptides.DeMixQ.csv
	/usr/bin/python3.6 /home/vitalv/cyano/label-free-demixq-diffacto-pipeline/run_diffacto.py  -peptide='/home/vitalv/cyano_dataset_20180112_YFP/DeMixQ/peptides.DeMixQ.csv' -db='/home/vitalv/database/Synechocystis_PCC6803_protein_sequences_YFP.fasta' -samples='/home/vitalv/cyano_dataset_20180112_YFP/diffacto/samples.csv' -ref=False -out='/home/vitalv/cyano_dataset_20180112_YFP/diffacto/peptides.DeMixQ.Diffacto.tsv' -out_type='peptide'
	'''


if __name__ == '__main__':

	main()
