
import os
import sys
import subprocess as sp
from multiprocessing import Pool
import datetime
import re

pipeline_folder = os.getcwd()

jobParams_ = dict()
jobParams_["allocation"] = "snic2016-7-94"
#This run times work well for lfq_workflow.py using featureFinderId instead of FeatureFinderCentroided (original DeMixQ)
#FeatureFInderCentroided generates many more features, lasts much longer, generates much larger files -> Pump those numbers up 
jobParams_["featureFinder"] = {"numCores": 10, "minutesRuntime": 60}
jobParams_["idMapper"] = {"numCores": 1, "minutesRuntime": 10}
jobParams_["featureLinker"] = {"numCores": 10, "minutesRuntime": 60}
jobParams_["mapAligner"] = {"numCores": 6, "minutesRuntime": 30}
jobParams_["textExporter"] = {"numCores": 1, "minutesRuntime": 10}
jobParams_["consensus2edta"] = {"numCores": 1, "minutesRuntime": 10}
jobParams_["eicExtractor"] = {"numCores": 6, "minutesRuntime": 60}
jobParams_["mapRTTransformer"] = {"numCores": 6, "minutesRuntime": 60}


'''
Make ProcessHandler a singleton class:
A class that is only instatiated once. 
executeCmd is then a method of this class
but then create LocalProcessHandler and SlurmProcessHandler as subclasses that inherit from ProcessHandler
and these two subclasses deal with the appropiate cmd
NOte: 
Singleton does not mean that it can't be called as many times as you like,
but everytime it will return the same object!
Create a Metaclass Singleton.
Then make ProcessHandler an instance of that metaclass
http://stackoverflow.com/questions/6760685/creating-a-singleton-in-python

'''


class Singleton(type):
	
	'''metaclass for ProcessHandler'''

	_instances = {}

	def __call__(cls, *args, **kwargs):
		#when I call ProcessHandler(), cls is: <class 'process_handler.ProcessHandler'>
		#print "cls is: %s" % cls    #cls is: <class '__main__.ProcessHandler'> when I call ProcessHandler()
		if cls not in cls._instances:
			#print "First time ProcessHandler() is called. cls._instances is %s" % cls._instances
			cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
			#print "super(Singleton, cls) is %s" % super(Singleton, cls) #<super: <class 'Singleton'>, <Singleton object>>
			#print "super(Singleton, cls).__call__(*args, **kwargs) is %s: " % super(Singleton, cls).__call__(*args, **kwargs)
		#else: #need this else if I want to do smthn on __init__ when calling ProcessHandler, otherwise delete
		#	#print "ProcessHandler() has already been called. cls._instances is %s" % cls._instances
		#	cls._instances[cls].__init__(*args, **kwargs) #In case I want to run __init__ everytime the class is called (so I can check if it's type Local or type Slurm)
		#	#print "cls._instances[cls].__init__(*args, **kwargs) is %s" % cls._instances[cls].__init__(*args, **kwargs)
		return cls._instances[cls]


class ProcessHandler(object):

	'''
	There will only be one instance of ProcesHandler
	Inherits from Singleton metaclass: ProcessHandler class is an instance of Singleton clss
	'''

	__metaclass__ = Singleton

	#def __init__(self): #__init__  is used when the class is called to initialize the instance: ph = ProcessHandler()
	#	print "ProcessHandler() intialized"
	#	self.getType() #wont work here __init__ can only return self or None

	def __call__(self): #__call__ method is called when the instance is called: ph()
		return self.getType()# = self.getType

	def getType(self):
		if sp.call("sinfo -o %v", shell=True) == 0:
			return SlurmProcessHandler()
		else: 
			print("slurm not found in system")
			return LocalProcessHandler()


	def getPool(self): #Note: should this also be in SlurmProcessHandler?
		return MyPool();



class MyPool:

	def __init__(self, processes=1):
		self.pool = Pool(processes)
		self.results = []

	def applyAsync(self, f, args):
		r = self.pool.apply_async(f, args)
		self.results.append(r)

	def checkPool(self):
		try:
			for res in self.results:
				res.get(0xFFFFFFFF)
			self.pool.close()
			self.pool.join()
			return self.results#try get submitted job id in msgfplus.py:pool.applyAsync(processWithMSGFPlus
		except KeyboardInterrupt:
			print "Caught KeyboardInterrupt, terminating workers"
			self.pool.terminate()
			self.pool.join()
			sys.exit()



class LocalProcessHandler(ProcessHandler):

	#__metaclass__ = Singleton

	def executeCmd(self, cmd, jobType, dependent_on_job_ids=[]):
		print ""
		print cmd 
		print ""
		#proc = sp.Popen(cmd, stdout=sp.PIPE, shell=True)
		proc = sp.Popen(cmd, shell=True)
		#print proc
		out,err = proc.communicate()
		return out, err



class SlurmProcessHandler(ProcessHandler):

	#def __init__(self):
	#	#super(SlurmProcessHandler, self).__init__()
	#	#self.slurm_version = self.get_slurm_version()
	#	print "detected slurm installation"


	#def get_slurm_version(self):
	#	out,err = sp.Popen("sinfo -o %v", shell=False, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
	#	#return "Slurm installation detected: \n" + out


	def executeCmd(self, cmd, jobType, dependent_on_job_ids=[]):
		jobType, scriptName = self.writeScript(cmd, jobType)
		submitted_job_id = self.submitJob(scriptName, jobType, dependent_on_job_ids)
		return submitted_job_id


	def writeScript(self, cmd, jobType):

		if re.search("-in \/.+", cmd): 
			file_name = jobType + "_" + re.search("-in \/.+", cmd).group().split("/")[-1].split(".")[0] + ".sh"
		else:
			file_name = jobType + ".sh"

		job_file = open(os.path.join(os.path.join(pipeline_folder,jobType), file_name), 'w+')

		runTimeString = str(datetime.timedelta(minutes=jobParams_[jobType]["minutesRuntime"])).replace(" days, ","-")  
		
		job_file.write("#!/bin/bash -l\n")
		job_file.write("#SBATCH -A %s\n" % jobParams_["allocation"])
		if jobParams_[jobType]["numCores"] > 15:
			job_file.write("#SBATCH -p node\n")
		else:
			job_file.write("#SBATCH -p core\n")
		job_file.write("#SBATCH -n %d\n" % (jobParams_[jobType]["numCores"]))
		job_file.write("#SBATCH -t %s\n" % (runTimeString))
		job_file.write("#SBATCH -J %s\n" % (jobType))
		job_file.write("#SBATCH -o %s" % os.path.join(pipeline_folder,jobType)+"/slurm-%j.out" ) #Note: folder and jobType must be same!

		job_file.write("\n%s\n" % (cmd))

		job_file.close()
		
		os.chmod(job_file.name,0755)
		
		return jobType, job_file.name #scriptName


	def submitJob(self, scriptName, jobType, dependent_on_job_ids):
		
		if not dependent_on_job_ids:
			cmd = "sbatch %s" % (scriptName)
		else:
			if isinstance(dependent_on_job_ids, list) and len(dependent_on_job_ids)> 1:
				dependent_on_job_ids = (":").join(dependent_on_job_ids)
			cmd = "sbatch --depend=afterok:%s %s"%(dependent_on_job_ids, scriptName)
		proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
		out,err = proc.communicate()
		rc = proc.returncode

		print "\n%s\n"%cmd

		submitted_job_id = re.search("Submitted batch job \d+", out).group().split("job ")[1]

		return submitted_job_id

		sys.stdout.flush()
		#if dryRun_:
		#rc = 0
		#else:
		#rc = subprocess.call(cmd, shell=True)

		if rc == 1:
			print("Error while submitting job %s" % scriptName)
			return 0


