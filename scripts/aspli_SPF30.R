###############################################################################
#                                                                             #
# This R script analyzes differential gene expresion and alternative splicing #
# between Col-0 and spf30-1                                                   #
# It was used to analyse the data described in Romanowski A. et al., 2020     #
#                                                                             #
###############################################################################

#########################################################################
# Author: Andrew Romanowski (aromanowski@leloir.org) - Yanovsky Lab     #
#########################################################################

### Notes
# Many lines of codes are commented. Uncomment them as needed for your usage

#########################################################################
#                                                                       #
#                     Requires                                          #
#                                                                       #
#########################################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install()}
if(!require(ASpli)) BiocManager::install("ASpli")
if(!require(GenomicFeatures)) BiocManager::install("GenomicFeatures")
if(!require(GenomicRanges)) BiocManager::install("GenomicRanges")
if(!require(AnnotationDbi)) BiocManager::install("AnnotationDbi")
if(!require(BiocParallel)) BiocManager::install("BiocParallel")

#########################################################################
#                                                                       #
#                     Includes                                          #
#                                                                       #
#########################################################################

library(ASpli)
library(BiocParallel)
library(GenomicFeatures)
library(GenomicRanges)
library(AnnotationDbi)

###############################################
#             Begin Analysis                  #
###############################################

# set base directory
basedir <- "c:/SPF30/comp"
setwd(basedir)

# Make the file genome sqlite from the transcriptome, using AtRTD2 transcriptome genes.gtf
# This GTF was edited so that it had the correct chromosome naming
#genome<-makeTxDbFromGFF("AtRTD2_edited.gtf", format = "gtf")
#saveDb(genome, file="genomeAtRTD2.sqlite")
#features<-binGenome(genome)
#save(features, file="featuresAtRTD2.sqlite")

# Load the genome
genome <- loadDb(file="genomeAtRTD2.sqlite")
# saveDb(genome, file = "genomeAtRTD2.sqlite")
load(file="featuresAtRTD2.RData")
# save(features, file = "featuresAtRTD2.RData")

###########################################################
# Summary of AtRTD2 features:                             #
#                                                         #
# * Number of extracted Genes = 34212                     #
# * Number of extracted Exon Bins = 238662                #
# * Number of extracted intron bins = 178027              #
# * Number of extracted trascripts = 82190                #
# * Number of extracted junctions = 151944                #
# * Number of AS bins (not include external) = 41863      #
# * Number of AS bins (include external) = 41941          #
# * Classified as:                                        #
#   ES bins = 1686	(4%)                                  #
#   IR bins = 13033	(31%)                                 #
#   Alt5'ss bins = 4244	(10%)                             #
#   Alt3'ss bins = 7683	(18%)                             #
#   Multiple AS bins = 15217	(36%)                       #
#   classified as:                                        #
#                 ES bins = 1627	(11%)                   #
#                 IR bins = 5060	(33%)                   #
#                 Alt5'ss bins = 2941	(19%)               #
#           			Alt3'ss bins = 5001	(33%)               #
###########################################################

############################################################
# Create a target file with description of the experiment  #
#                                                          #
# sample	bam	condition                                    #
# Col_LL_A	Col_LL_A.bam	ctrl	                           #
# Col_LL_B	Col_LL_B.bam	ctrl	                           #
# Col_LL_C	Col_LL_C.bam	ctrl	                           #
# SMN_LL_A	SMN_LL_A.bam	spf30-1	                         #
# SMN_LL_B	SMN_LL_B.bam	spf30-1	                         #
# SMN_LL_C	SMN_LL_C.bam	spf30-1	                         #
############################################################
targets <- read.table("targets.txt", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
bam <- loadBAM(targets,cores=8)
save(bam, file="spf30-1_BAM.RData")

##############################################
# OUTPUT module                              #
# Count tables, PSI, PIR (AS with junctions) #
##############################################
# Get raw counts from the BAM data
# l is read length from the experiment. Here it is set to 100bp
counts <- readCounts(features, bam, cores=7, readLength = 100, targets=targets, maxISize = 5000)
save(counts, file="counts_spf30-1.Rdata") # Save RAW count data

# next set the conditions, the group and pair as in the target file (column condition)
condition <- c(rep("ctrl",3), rep("spf30-1",3))
group <- factor(c(rep("ctrl",3), rep("spf30-1",3)))
pair <- c("ctrl","spf30-1")

# Get splicing events counts with ASpli
as <- AsDiscover(counts, targets, features, bam, readLength = 100, threshold = 5, cores = 1)
save(as, file="as_spf30-1.Rdata")

#######################################################
# OUTPUT module                                       #
# Differential expression and differential bin usage  #
#######################################################
du <-DUreport(counts = counts, targets = targets, forceGLM = TRUE)
save(du, file="du_spf30-1.RData")

######################
#    Print results   #
######################
writeAll(counts = counts, du = du, as = as, output.dir = "spf30-1_Results")
