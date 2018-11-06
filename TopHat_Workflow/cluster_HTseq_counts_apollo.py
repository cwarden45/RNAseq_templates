import sys
import re
import os
import commands

jobPrefix="cwHT"
max_concurrent = 12

parameterFile = "parameters.txt"
finishedSamples = ()

alignmentFolder = ""
gtf_file = ""
lncRNAgtf = ""
strandType = ""
email = ""
libType = ""

inHandle = open(parameterFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	param = lineInfo[0]
	value = lineInfo[1]
	
	if param == "Alignment_Folder_MAC":
		alignmentFolder = value
		
	if param == "txGTF_MAC":
		gtf_file = value

	if param == "lncRNA_GTF_MAC":
		lncRNAgtf = value		

	if param == "strand":
		strandType = value
		
	if param == "Cluster_Email":
		email = value

	if param == "Read_Pairing":
		libType = value	
		
if (strandType != "yes") and (strandType != "no") and (strandType != "reverse"):
	print "Need to provide HT-Seq mapping for strand: " + strandType
	sys.exit()
		
fileResults = os.listdir(alignmentFolder)

jobCount = 0
jobHash={}

for file in fileResults:
	result = re.search(".bam$",file)
	fullPath = os.path.join(alignmentFolder, file)
	
	if result:
		sample = re.sub(".bam$","",file)
		sortResult = re.search(".name.sort.bam",file)
		
		if not sortResult:
			jobCount += 1
			jobName = jobPrefix+str(jobCount)
			print jobName
		
		if (sample not in finishedSamples) and (not sortResult):
			print sample
	
			shellScript = "htseq_" + sample + ".sh"
			
			outHandle = open(shellScript, "w")
			text = "#!/bin/bash\n"
			text = text + "#SBATCH -J "+jobName+"\n"
			text = text + "#SBATCH --mail-type=ALL\n"
			text = text + "#SBATCH --mail-user="+email+"\n"
			text = text + "#SBATCH -n 1\n"#one thread
			text = text + "#SBATCH -N 1\n"
			text = text + "#SBATCH --mem=4g\n"
			text = text + "#SBATCH --time=12:00:00\n"
			text = text + "#SBATCH --output="+jobName+".log\n\n"
			
			##module load htslib/1.6
			
			#text = text + "module load Python/2.7.14-foss-2017a\n\n"
			#code above is supposed to be sufficient for htseq-count, but I need to use alternative set-up below
			text = text + "module load Python/3.6.1-foss-2017a\n"
			#installed with "python setup.py install --user", with git code
			#with executable /net/isi-dcnl/ifs/user_data/Seq/cwarden/More_Programs/htseq/build/scripts-3.6/htseq-count
			text = text + "module load samtools/1.6\n\n"
			outHandle.write(text)

			#when I've tested single-end RNA-Seq, I've gotten the same results for name and position sorted .bam (and manual says parameter is ignored for single-end data)
			#However, for paired-end RNA-Seq, there can be differences, and the name-sorted .bam seems to be a little more accurate
			if libType == "SE":
				countsFile = sample + "_gene_counts.txt"
				text = "/net/isi-dcnl/ifs/user_data/Seq/cwarden/More_Programs/htseq/build/scripts-3.6/htseq-count -f bam -r pos -s " + strandType + " " + fullPath + " " + gtf_file + " > " + countsFile + "\n"
				outHandle.write(text)

				countsFile = sample + "_lncRNA_counts.txt"
				#switch to " -i gene_id " for MiTranscriptome
				text = "/net/isi-dcnl/ifs/user_data/Seq/cwarden/More_Programs/htseq/build/scripts-3.6/htseq-count -f bam -r pos -s " + strandType + " -i gene_name " + fullPath + " " + lncRNAgtf + " > " + countsFile + "\n"
				#outHandle.write(text)	
			elif libType == "PE":
				nameSortedBam = sample + ".name.sort.bam"
				sortPrefix = re.sub(".bam$","",nameSortedBam)
				text = "/opt/SAMtools/1.6/bin/samtools sort -n " + fullPath + " " + sortPrefix + "\n"
				outHandle.write(text)

				countsFile = sample + "_gene_counts.txt"
				text = "/net/isi-dcnl/ifs/user_data/Seq/cwarden/More_Programs/htseq/build/scripts-3.6/htseq-count -f bam -s " + strandType + " " + nameSortedBam + " " + gtf_file + " > " + countsFile + "\n"
				outHandle.write(text)

				countsFile = sample + "_lncRNA_counts.txt"
				#switch to " -i gene_id " for MiTranscriptome
				text = "/net/isi-dcnl/ifs/user_data/Seq/cwarden/More_Programs/htseq/build/scripts-3.6/htseq-count -f bam -s " + strandType + " -i gene_name " + nameSortedBam + " " + lncRNAgtf + " > " + countsFile + "\n"
				#outHandle.write(text)

				text = "rm " + nameSortedBam + "\n"
				outHandle.write(text)
			else:
				print "'Read_Pairing' must be single-end('SE') or paired-end ('PE')"
				sys.exit()
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
				status, outtext = commands.getstatusoutput(cmd)

				numResult = re.search("(\d+)",outtext)
				if numResult:
					jobnum = numResult.group(1)
				else:
					print "Modify code to parse output: " + outtext
					sys.exit()
				
				jobHash[jobName] = str(jobnum)
				print jobName + "-->"+jobnum		
