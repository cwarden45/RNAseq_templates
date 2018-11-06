import sys
import re
import os
import commands

jobPrefix="CW"
max_concurrent = 12

parameterFile = "parameters.txt"
finishedSamples = ()

threads = ""
ref = ""
alignmentFolder = ""
readsFolder = ""
strandType = ""
email = ""

inHandle = open(parameterFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	param = lineInfo[0]
	value = lineInfo[1]

	if param == "Reads_Folder_MAC":
		readsFolder = value
	
	if param == "Alignment_Folder_MAC":
		alignmentFolder = value
		
	if param == "Bowtie2_Ref":
		ref = value
		
	if param == "Threads":
		threads = value

	if param == "strand":
		strandType = value
		
	if param == "Cluster_Email":
		email = value
		
if (threads == "") or (threads == "[required]"):
	print "Need to enter a value for 'Threads'!"
	sys.exit()

if (strandType == "") or (strandType == "[required]"):
	print "Need to enter a value for 'strand'!"
	sys.exit()

if (email == "") or (email == "[required]"):
	print "Need to enter a value for 'Cluster_Email'!"
	sys.exit()
	
if (ref == "") or (ref == "[required]"):
	print "Need to enter a value for 'Bowtie2_Ref'!"
	sys.exit()

if (readsFolder == "") or (readsFolder == "[required]"):
	print "Need to enter a value for 'Reads_Folder_MAC'!"
	sys.exit()
	
if (alignmentFolder == "") or (alignmentFolder == "[required]"):
	print "Need to enter a value for 'Alignment_Folder_MAC'!"
	sys.exit()
	
fileResults = os.listdir(readsFolder)

submitAll = "master_tophat_queue.sh"
masterHandle = open(submitAll,"w")
text = "#!/bin/bash\n"
masterHandle.write(text)

jobCount = 0
jobHash={}

for file in fileResults:
	result = re.search("(.*)_S\d+_L\d{3}_R1_001.fastq$",file)
	
	if result:
		sample = result.group(1)
		jobCount += 1
		jobName = jobPrefix+str(jobCount)
		print jobName
			
		if (sample not in finishedSamples):
			print sample
			shellScript = sample + ".sh"
			text = "sbatch " + shellScript + "\n"
			masterHandle.write(text)

			outHandle = open(shellScript, "w")
			text = "#!/bin/bash\n"
			text = text + "#SBATCH -J "+jobName+"\n"
			text = text + "#SBATCH --mail-type=ALL\n"
			text = text + "#SBATCH --mail-user="+email+"\n"
			text = text + "#SBATCH -n "+threads+"\n"#one thread
			text = text + "#SBATCH -N 1\n"
			text = text + "#SBATCH --mem=8g\n"
			text = text + "#SBATCH --time=48:00:00\n"
			text = text + "#SBATCH --output="+jobName+".log\n\n"

			text = text + "module load samtools/1.6\n"
			#text = text + "export LD_LIBRARY_PATH=/net/isi-dcnl/ifs/user_data/Seq/cwarden/bowtie2-2.3.2/tbb2017_20170604oss/lib/intel64/gcc4.7:$LD_LIBRARY_PATH\n"
			text = text + "export PATH=/net/isi-dcnl/ifs/user_data/Seq/cwarden/bowtie2-2.3.2:$PATH\n\n"			
			outHandle.write(text)

				
			outputSubfolder = alignmentFolder +"/" + sample
			text = "mkdir " + outputSubfolder + "\n"
			outHandle.write(text)
										
			read1 = readsFolder + "/" + file

			tophatStrand = ""
			if strandType == "no":
				tophatStrand = "fr-unstranded"
			elif strandType == "reverse":
				tophatStrand = "fr-firststrand"
			elif strandType == "yes":
				tophatStrand = "fr-secondstrand"
			else:
				print "Need to provide TopHat mapping for strand: " + strandType
				sys.exit()
			
			text = "/net/isi-dcnl/ifs/user_data/Seq/cwarden/tophat-2.1.1.Linux_x86_64/tophat2 -o " + outputSubfolder + " -p " + threads + " --no-coverage-search --library-type " + tophatStrand + " " + ref + " " + read1 + "\n" 
			outHandle.write(text)
									
			topHatBam = outputSubfolder + "/accepted_hits.bam"																			
			userBam = alignmentFolder + "/" + sample + ".bam"
			
			text = "mv " + topHatBam + " " + userBam +"\n"
			outHandle.write(text)

			text = "/opt/SAMtools/1.6/bin/samtools index " + userBam + "\n"
			outHandle.write(text)

			text = "gzip " + read1
			outHandle.write(text)
			
			outHandle.close()
			
			if jobCount > max_concurrent:
				#test code from https://hpc.nih.gov/docs/job_dependencies.html
				
				depJobName = jobPrefix+str(jobCount-max_concurrent)
				depJobID = ""
				if depJobName in jobHash:
					depJobID=jobHash[depJobName]
				else:
					print "Error mapping ID for " + depJobName
					sys.exit()
			
				cmd = "sbatch --depend=afterany:"+ depJobID + " " + shellScript
				status, outtext = commands.getstatusoutput(str(cmd))
				
				numResult = re.search("(\d+)",outtext)
				if numResult:
					jobnum = numResult.group(1)
				else:
					print "Modify code to parse output: " + outtext
					sys.exit()
				
				jobHash[jobName] = jobnum
				print jobName + "-->"+jobnum
			else:
				cmd = "sbatch " + shellScript
				status, outtext = commands.getstatusoutput(str(cmd))

				numResult = re.search("(\d+)",outtext)
				if numResult:
					jobnum = numResult.group(1)
				else:
					print "Modify code to parse output: " + outtext
					sys.exit()
				
				jobHash[jobName] = str(jobnum)
				print jobName + "-->"+jobnum			
