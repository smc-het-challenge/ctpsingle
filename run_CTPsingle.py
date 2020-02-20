import sys
import os
import getopt

MAX_MUTATIONS = 10000000 # FOR DEBUGGING PURPOSES. MAKE SURE EITHER TO REMOVE IT OR TO SET TO SOME LARGE VALUE IN THE FINAL VERSION

gender		= "unknown"
multiplier	= 2
minimumCoverage = 20


class CNA:
	def __init__(self):
		self.chrom = ""
		self.start = 0
		self.end   = 0
		self.major_allele = 0
		self.minor_allele = 0
		self.frac  = 0.0



class SNV:
	def __init__(self, chromosome, position, refNucleotide, altNucleotide):
		self.chrom  = chromosome
		self.pos    = position
		self.refNuc = refNucleotide
		self.altNuc = altNucleotide
		self.cluster = "NA"

	def setCluster(self, cluster):
		assert self.cluster == "NA", "ERROR. Cluster of mutation already set. Can not be changed."
		self.cluster = cluster

	def __eq__(self, other):
		return self.chrom == other.chrom and self.pos == other.pos and self.refNuc == other.refNuc and self.altNuc == other.altNuc
	
	def __str__(self):
		return "SNV: chrom=" + self.chrom + ", pos=" + str(self.pos) + ", refNuc=" + self.refNuc + ", altNuc=" + self.altNuc


class Cluster:
	def __init__(self, clusterID, cellFraction):
		self.ID = clusterID
		self.cellFraction = cellFraction
		self.numMutations = 0
	
	def increaseNumMutations(self):
		self.numMutations += 1
	
	def getNumMutations(self):
		return self.numMutations

	def getCellFraction(self):
		return self.cellFraction
		
	def __eq__(self, other):
		return self.ID == other.ID
	
		


'''
    Use neutral + unknown for CTPsingle input.
    chromosome and position arguments are chromosome and position of SNV
'''
def SNV_copy_number_status(chromosome, position, CNA_array):

	for CNA_event in CNA_array: 
		if position >= CNA_event.start and position <= CNA_event.end and chromosome == CNA_event.chrom:
			#print(chromosome + "\t" + str(position) + ":\t" + CNA_event.chrom + " " + str(CNA_event.start) + "," + str(CNA_event.end))
			if CNA_event.frac != 1:
				return "subclonal"
			elif CNA_event.major_allele == 1 and CNA_event.minor_allele == 0: # clonal single copy deletion
				return "deletion"
			elif CNA_event.major_allele == 1 and CNA_event.minor_allele == 1:
				return "neutral"
			else:
				return "clonal"
	
	return "unknown" # this means that we did not find any CNA event in the array CNV_array spanning given SNV




############################################ Reading arguments ####################################
pathBattenbergFile = ""
pathSNVFile	   = ""
outputFilesPrefix  = "" # program will first prepare .frq file and then run CTPsingle and report several output files

helpMessage =  "\nArguments are the following:\n"
helpMessage += "\t-c battenberg CNV file\n"
helpMessage += "\t-s SNV file\n"
helpMessage += "\t-o [optional argument] prefix of the output files (if not set, prefix will be obtained from path to Battenberg file by trimming everything after .battenberg\n"
helpMessage += "\n-------- EXAMPLE:\t"
helpMessage += "python run_CTPsingle.py -c example-data/P1-noXY.battenberg.txt -s example-data/P1-noXY.mutect.vcf -o example-data/subchallenge"
helpMessage += "\n"
       
if sys.argv[1] == "-h":
	print(helpMessage)
	sys.exit(2)

try:
	myopts, args = getopt.getopt(sys.argv[1:],"c:s:o:")

except getopt.GetoptError as e:
	print (str(e))
	print(helpMessage)
	sys.exit(2)

for option, argument in myopts:
	if option == '-c':
        	pathBattenbergFile = argument.strip()
	elif option == '-s':
		pathSNVFile = argument.strip()
	elif option == '-o':
		outputFilesPrefix = argument.strip()


assert os.path.exists(pathBattenbergFile), "ERROR. Path to Battenberg file set to " + pathBattenbergFile + ". File does not exist."
assert os.path.exists(pathSNVFile), "ERROR. Path to SNV file set to " + pathSNVFile + ". File does not exist."


if outputFilesPrefix == "":
	outputFilesPrefix = os.sep.join(pathBattenbergFile.split(os.sep)[:-1]) + os.sep

	for component in pathBattenbergFile.strip().split(os.sep)[-1].split("."): # take Battenberg file name and split it by "."
		if component.upper() == "BATTENBERG":
			break
		else:
			outputFilesPrefix += component + "."

	outputFilesPrefix = outputFilesPrefix.rstrip(".") + ".CTPsingle" # consider removing + ".CTPsingle" if necessary


############################################ Reading BATTENBERG FILE ##############################
battenberg_CNA_segments	= []
diploidGenomeLength	= 0  # length of diploid portion of the genome
nonDiploidGenomeLength  = 0  # length of non-diploid portion of the genome

battenbergFile = open(pathBattenbergFile, 'r')
battenbergFile.readline() # read header line

for line in battenbergFile:
	lineColumns = line.strip().split()

	CNAevent = CNA()
	CNAevent.chrom = lineColumns[0]
	CNAevent.start = int(lineColumns[1])
	CNAevent.end   = int(lineColumns[2])
	CNAevent.major_allele = int(lineColumns[7])
	CNAevent.minor_allele = int(lineColumns[8])
	CNAevent.frac  = float(lineColumns[9])

	if CNAevent.chrom == "X" or CNAevent.chrom == "Y":
		continue

	battenberg_CNA_segments.append(CNAevent)

	if CNAevent.minor_allele == 1 and CNAevent.major_allele == 1 and CNAevent.frac == 1:
		diploidGenomeLength += CNAevent.end  - CNAevent.start + 1
	else:
		nonDiploidGenomeLength += CNAevent.end - CNAevent.start + 1

battenbergFile.close()

fracNonDiploidGenome = (float(nonDiploidGenomeLength))/(diploidGenomeLength + nonDiploidGenomeLength)



################################################# Reading SNV file #########################################################
all_SNVs  = []
used_SNVs = []
copy_neutral_SNVs = 0

pathFrqFile = outputFilesPrefix + ".frq"
Frq_file = open(pathFrqFile, "w")
Frq_file.write("Chromosome Position Mutant Reference Mcount Rcount Multiplier Gender\n")
SNV_file = open(pathSNVFile, "r")

for line in SNV_file:
	if line[0] == "#":
		continue

	lineColumns = line.split("\t")
	SNV_chr	    = lineColumns[0]
	SNV_pos	    = int(lineColumns[1])
	SNV_refNuc  = lineColumns[3]
	SNV_altNuc  = lineColumns[4]
	all_SNVs.append(SNV(SNV_chr, SNV_pos, SNV_refNuc, SNV_altNuc))

    	if SNV_chr.upper() == "X" or SNV_chr.upper() == "Y":
		continue


        SNV_batt_CN_status = SNV_copy_number_status(SNV_chr, SNV_pos, battenberg_CNA_segments)
	if SNV_batt_CN_status in ["subclonal", "deletion", "clonal"]:
		#print(SNV_chr + "\t" + str(SNV_pos) + "\t" + SNV_batt_CN_status)
		continue
	if SNV_batt_CN_status not in ["neutral", "unknown"]:
		print("ERROR. Unknown SNV_batt_CN_status " + SNV_batt_CN_status)
		continue

	copy_neutral_SNVs += 1
	

       	FILTER = lineColumns[6];
	if FILTER != "PASS":
		continue
	
	INFO = lineColumns[7]
	if "SOMATIC" not in INFO:
		continue


	tumor_reference_coverage = -1
	tumor_altered_coverage   = -1

	FORMAT   = lineColumns[8];
	fieldIDs = FORMAT.split(":")
	
	for i in range(len(fieldIDs)):
		if fieldIDs[i] == "AD":
			tumorColumn = lineColumns[10].strip()
			tumor_AD = tumorColumn.split(":")[i]
			tumor_reference_coverage = int(tumor_AD.split(",")[0])
			tumor_altered_coverage   = int(tumor_AD.split(",")[1])
			break

	if tumor_altered_coverage + tumor_reference_coverage < minimumCoverage:
		continue


	Frq_file.write("{} {} {} {} {} {} {} {}\n".format(SNV_chr, SNV_pos, SNV_altNuc, SNV_refNuc, tumor_altered_coverage, tumor_reference_coverage, multiplier, gender))
	used_SNVs.append(SNV(SNV_chr, SNV_pos, SNV_refNuc, SNV_altNuc))
	if len(used_SNVs) > MAX_MUTATIONS:
		break

Frq_file.close()




# First make sure that no SNV appears multiple times in the input file (mutations at the same position with different variant nucleotide are ok)
repeatingMutations = False
for i in range(len(all_SNVs)):
	for j in range(i):
		if all_SNVs[i] == all_SNVs[j]:
			repeatingMutations = True
			break

if repeatingMutations == False:
	os.system("Rscript CTPsingle_pancancer.R -f " + pathFrqFile + " -o " + outputFilesPrefix + " -m " + "GammaAdjMatrices")


os.remove(pathFrqFile)

# If CTPsingle didn't report any output just exit without doing anything further
if os.path.exists(outputFilesPrefix + ".CTPsingle.purity") == False:

	for subchallenge in ["1A", "1B", "1C", "2A", "3A"]:
		subchallengeOutputFile = open(outputFilesPrefix + "." + subchallenge + ".txt", "w")
		subchallengeOutputFile.write("NA")
		subchallengeOutputFile.close()

	sys.exit(2)


# read purity value from CTPsingle output file, remove that file and create new file
purityFile = open(outputFilesPrefix + ".CTPsingle.purity", "r")
purity = float(purityFile.readline().rstrip().split()[1])
purityFile.close()
os.remove(outputFilesPrefix + ".CTPsingle.purity")
file_1A = open(outputFilesPrefix + ".1A.txt", "w")
file_1A.write(str(purity) + "\n")
file_1A.close()


clusters = {}
clusterAssignmentsFile = open(outputFilesPrefix + ".CTPsingle.clusters", "r")
clusterAssignmentsFile.readline()
lines = clusterAssignmentsFile.readlines()
assert len(lines) == len(used_SNVs), "ERROR. Cluster assignments file reported by CTPsingle has different number of lines than input .frq file."

for i in range(len(lines)):
	line = lines[i]
	lineColumns = line.strip().split()
	assert used_SNVs[i].chrom  == lineColumns[0], "ERROR. Position of mutation in .clusters file differs from the position in the corresponding .frq file! " + str(used_SNVs[i])
	assert used_SNVs[i].pos    == int(lineColumns[1]), "ERROR. Position of mutation in .clusters file differs from the position in the corresponding .frq file! " + str(used_SNVs[i])
	assert used_SNVs[i].altNuc == lineColumns[2], "ERROR. Variant nucleotide in .clusters file differs from variant nucleotide in the corresponding .frq file! "  + str(used_SNVs[i])
	assert used_SNVs[i].refNuc == lineColumns[3], "ERROR. Reference nucleotide in .clusters file differs from reference nucleotide in the corresponding .frq file! " + str(used_SNVs[i]) 

	clusterID = lineColumns[4]
	if clusterID == "0":
		clusterID = "GHOST_CLUSTER"
	used_SNVs[i].setCluster(clusterID)
	if clusters.has_key(clusterID) == False:
		clusters[clusterID] = Cluster(clusterID, float(purity * float(lineColumns[5])))

	clusters[clusterID].increaseNumMutations()


sizeGhostCluster = 0
if clusters.has_key("GHOST_CLUSTER"):
	sizeGhostCluster = clusters["GHOST_CLUSTER"].getNumMutations()
	del clusters["GHOST_CLUSTER"]

numClusters = len(clusters)
file_1B = open(outputFilesPrefix + ".1B.txt", "w")
file_1B.write(str(numClusters) + "\n")
file_1B.close()


clusterAssignmentsFile.close()
os.remove(outputFilesPrefix + ".CTPsingle.clusters")

'''
DESCRIPTION:
In its output cluster assignment file, CTPsingle associates each mutation with single cluster. There might exist cluster with ID 0 which represents ghost cluster.
When reporting output we associate numbers from 1 onwards to clusters. Ghost cluster is removed from consideration by this timepoint.
So, internally all the time we work with cluster IDs as in CTPsingle output and only when we report final output files, instead of writing particular clusterID we write
outputClusterID[clusterID] which transforms ID of this cluster into ID required in the output (as specified by DREAM CHALLENGE).

Note that, in tree files reported by CTPsingl, 0 is unrelated to ghost cluster and is used to denote the root node (which is population of healthy cells).
'''

clusterIDs = sorted([int(x) for x in clusters.keys()]) # sorting is not essential, but we did it to make debugging easier
for i in range(numClusters):
	clusterIDs[i] = str(clusterIDs[i])


outputClusterID = {} # outputClusterID[x] = y means: if CTPsingle used x as ID for cluster, this ID will be printed as y in the output
currentIndex    = 1  # output cluster IDs should be 1, 2, ... so that tree representation works properly
for clusterID in clusterIDs:
	outputClusterID[clusterID] = currentIndex 
	currentIndex += 1
outputClusterID["GHOST_CLUSTER"] = currentIndex


file_1C = open(outputFilesPrefix + ".1C.txt", "w")
for i in range(numClusters):
	clusterID    = clusterIDs[i]	
	numMutations = clusters[clusterID].getNumMutations()
	cellFraction = clusters[clusterID].getCellFraction()
	file_1C.write("\t".join([str(outputClusterID[clusterID]), str(numMutations),  str(cellFraction)]) + "\n")

numUnclusteredMutations = len(all_SNVs) - len(used_SNVs) + sizeGhostCluster
if numUnclusteredMutations > 0:
	file_1C.write("\t".join([str(outputClusterID["GHOST_CLUSTER"]), str(numUnclusteredMutations), "0"]) + "\n")

file_1C.close()



file_2A = open(outputFilesPrefix + ".2A.txt", "w")
index_usedSNVs = 0  # iterate through all SNVs. When we meet SNV used for clustering take and print its cluster, otherwise write NA. 

for index_allSNVs in range(len(all_SNVs)):
	if all_SNVs[index_allSNVs] == used_SNVs[index_usedSNVs]:
		clusterID = used_SNVs[index_usedSNVs].cluster
		if clusterID != "GHOST_CLUSTER":
			file_2A.write(str(outputClusterID[clusterID]) + "\n")
		else:
			file_2A.write(str(outputClusterID["GHOST_CLUSTER"]) + "\n")
		index_usedSNVs += 1
	else:
		file_2A.write(str(outputClusterID["GHOST_CLUSTER"]) + "\n")

file_2A.close()



currentTreeIndex   = 1 # index of the tree that is currently under consideration (among output trees). It is expected that output trees are indexed from 1 to #trees for tree of given size.
printedTreeIndex   = 1 # some of the trees only will be reported (depending on the score of the tree). This keeps index of the next tree to be reported. Index is used in filenames.
bestBranchingScore = -1
while True:
	pathCurrentTree = outputFilesPrefix + ".treeSize_" + str(numClusters) + "-treeIndex_" + str(currentTreeIndex) + ".CTPsingle.tree"
	if os.path.exists(pathCurrentTree) == False:
		break
	else:
		treeFile = open(pathCurrentTree, "r")
		treeFileLines = treeFile.readlines()
		treeFile.close()
		os.remove(pathCurrentTree)
		score = float(treeFileLines[0].strip().split()[-1])
		
		if score < 0.00001: # we know that optimal score is 0 so if score is not 0 just continue
			parent = {}
			numChildren = {} # numChildren[clusterID] with obvious meaning
			
			# cluster IDs in reported tree files are not consistent with these reported in clustering file. Resolve the ambiguity			
			treeClusterToTrueClusterID = {}
			treeClusterToTrueClusterID["0"] = "ROOT"
			for line in treeFileLines:
				lineColumns = line.strip().split()
				treeClusterID = lineColumns[-3]
				frequency = float(lineColumns[-2])*purity # need to adjust for purity to get cellular fraction
				for clusterID in clusterIDs:
					if abs(clusters[clusterID].getCellFraction() - frequency) < 0.01:	
						treeClusterToTrueClusterID[treeClusterID] = clusterID
				assert treeClusterToTrueClusterID.has_key(treeClusterID), "ERROR. Did not find true cluster ID matching to tree cluster ID " + treeClusterID

			for line in treeFileLines:
				if len(line.strip()) == 0:
					continue
				lineColumns     = line.strip().split()
				parentClusterID = lineColumns[-4]
			
				parentClusterID = treeClusterToTrueClusterID[lineColumns[-4]]
				childClusterID  = treeClusterToTrueClusterID[lineColumns[-3]]

				parent[childClusterID] = parentClusterID
				if numChildren.has_key(parentClusterID) == False:
					numChildren[parentClusterID] =  1
				else:
					numChildren[parentClusterID] += 1
			
			assert numChildren["ROOT"] == 1, "ERROR. Root degree greater than 1 in " + pathCurrentTree

			branchingScore = sum([(numChildren[x]-1) for x in list(numChildren.keys())])
			if branchingScore > bestBranchingScore:
				bestBranchingScore = branchingScore
				treeSize = sum([numChildren[x] for x in list(numChildren.keys())])
				assert treeSize == numClusters, "ERROR. Tree size (without 0 node which is population of healthy cells) differs from the number of clusters!"
				if os.path.exists(outputFilesPrefix + ".3A.txt"):
					os.remove(outputFilesPrefix + ".3A.txt")

				file_3A = open(outputFilesPrefix + ".3A.txt", "w")
				for row in range(1, treeSize+1):
					for clusterID in clusterIDs:
						if outputClusterID[clusterID] == row:
							if parent[clusterID] == "ROOT":
								file_3A.write(str(outputClusterID[clusterID]) + "\t" + "0")
								file_3A.write("\n")
							else:
								file_3A.write(str(outputClusterID[clusterID]) + "\t" + str(outputClusterID[parent[clusterID]]))
								file_3A.write("\n")
				file_3A.close()
				printedTreeIndex += 1		
		
	currentTreeIndex += 1	


if printedTreeIndex == 1: # means no tree was printed above
	file_3A = open(outputFilesPrefix + ".3A.txt", "w")
	file_3A.write("NA")
	file_3A.close()

print("All done")
