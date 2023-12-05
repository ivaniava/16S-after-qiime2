# Libraries ----
library(microDecon) # for decontamination using blanc and controls
library(dplyr)
library(tibble)
library(qiime2R)
library(tidyverse)
library(phyloseq)
library(DECIPHER)
library(phangorn)
library(MicEco)
library(openxlsx)
library(readxl)
library(phylosmith)
library(MicEco)
library(miaViz)
library(microbiome)
library(eulerr)
library(qpcR)
library(microbiomeutilities)

set.seed(16818)

# Samples with read count per ASV without normalization ----
ASV_table <- read.qza(dada2table_All.qza")
ASV_matrix <- ASV_qza$data
ASV_matrix <- as.data.frame(ASV_matrix)

samples_metadata <- read.table(".../files/FFPE_All.tsv", sep = '\t', header = T)
names(ASV_matrix) <- samples_metadata$feature.type.number[match(names(ASV_matrix), samples_metadata$sample.id)]

# Reorder file for microdecon: ASV_id, blanks/controls, Samples(ordered per group), taxa. ----
ASV_matrix <- select(ASV_matrix, "Blanc", "Control", 
                     "AD differentiated 1", "AD differentiated 2", "AD differentiated 3", ...)

#File for microdecon, for put back the rownames as first col ----
ASV_matrix <- rownames_to_column(ASV_matrix, var = "ASV_id") 

ASV_decon <- decon(ASV_matrix, numb.blanks = 2, numb.ind = c(30, 16, 36, 20, 3, 16, 17, 9, 8), taxa = F, runs = 2, 
                   thresh = 0.7, prop.thresh = 5e-05, regression = 0, low.threshold = 40, up.threshold = 400)

#Table of ASV count per sample decontaminated ----
ASV_clean <- ASV_decon$decon.table
ASV_clean$Mean.blank<-NULL
names(ASV_clean) <- samples_metadata$sample.id[match(names(ASV_clean), samples_metadata$feature.type.number)]
colnames(ASV_clean)[1] = "Feature.ID"
write.table(ASV_clean, file = "PQ00193_FFPE/FFPE_All/ASV/decontamination/ASV_clean.tsv", row.names = F, quote = F, sep = '\t')

#Tables of ASV removed and summarized results ----
ASV_removed <- ASV_decon$OTUs.removed
colnames(ASV_removed)[1] = "Feature.ID"
ASV_removed_reads <- ASV_decon$reads.removed
colnames(ASV_removed_reads)[1] = "Feature.ID"
ASV_removed_mean <- ASV_decon$mean.per.group 
colnames(ASV_removed_mean)[1] = "Feature.ID"
ASV_removed_sum <- ASV_decon$sum.per.group 
colnames(ASV_removed_sum)[1] = "Feature.ID"

#ASV with taxonomy to observe which bacterias were eliminated and remained in the analysis ----
ASV_qza_tax <- read_qza("PQ00193_FFPE/FFPE_All/taxonomy/taxonomy_silva_All.qza")
ASV_tax <- ASV_qza_tax$data
ASV_tax_removed <- merge(ASV_removed, ASV_tax[, c("Feature.ID", "Taxon")], by="Feature.ID") 
ASV_tax_removed <- select(ASV_tax_removed, "Feature.ID", "Taxon", everything())
ASV_tax_remained <- merge(ASV_clean, ASV_tax[, c("Feature.ID", "Taxon")], by="Feature.ID")
ASV_tax_remained2 <- ASV_tax_remained[, c(1, 157)]
write.table(ASV_tax_removed, file = "PQ00193_FFPE/FFPE_All/ASV/decontamination/ASV_tax_removed.tsv", row.names = F, quote = F, sep = '\t')
write.table(ASV_tax_remained2, file = "PQ00193_FFPE/FFPE_All/ASV/decontamination/ASV_tax_remained.tsv", row.names = F, quote = F, sep = '\t')

# Remove ASV ----
# For some diversity analysis, we can add a phylogenetic tree, this can be made aligning DNA sequences to a db (fragment insertion sepp) or 
# made a de novo tree. For this purpose, the file containing the sequences dada2seqs_AD_SQ_NAT.qza is uploaded and substract from it,
# the 18 ASV removed using the microdecon package in R (ASV_tax_removed_allgroups). 

ASV_tax_removed_allgroups <- ASV_tax_removed
ASV_tax_removed_allgroups <- ASV_tax_removed_allgroups[-c(3,4,6,12,15),]
ASV_dna_seqs <- read.table("PQ00193_FFPE/FFPE_All/ASV/dada2seqs_All.tsv", sep = '\t', header = T)
ASV_dna_seqs_clean <- anti_join(ASV_dna_seqs, ASV_tax_removed_allgroups, by = "Feature.ID")
write.table(ASV_dna_seqs_clean, file = "PQ00193_FFPE/FFPE_All/ASV/decontamination/ASV_dna_seqs_clean.tsv", row.names = F, quote = F, sep = '\t')

ASV_decontam <- read.table("PQ00193_FFPE/FFPE_AD_SQ_NAT/ASV/decontamination/ASV_clean.tsv", sep = '\t', header = T)
ASV_decontam <- tibble::column_to_rownames(ASV_decontam, var = "Feature.ID")
ASV =   otu_table(data.frame(ASV_decontam), taxa_are_rows = TRUE)

ASV_tax_decontam <- read.table("PQ00193_FFPE/FFPE_AD_SQ_NAT/ASV/decontamination/ASV_tax_remained.tsv", sep = '\t', header = T)
taxonomy <-  ASV_tax_decontam %>% separate(Taxon, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species"), sep=";") %>% column_to_rownames(var = 'Feature.ID') # specify that the "Feature ID" column of data are rownames
TAX =   tax_table(as.matrix(taxonomy))

metadata_All <- read.table("PQ00193_FFPE/manifest_files/FFPE_All.tsv", sep = '\t', header = T)
META = sample_data(data.frame(metadata_All, row.names = metadata_All$sample.id))

head(taxonomy)
head(taxa_names(TAX))
head(taxa_names(ASV))
head(sample_names(ASV))
head(sample_names(META))

#Create phyloseq object ----
ps <- phyloseq(ASV, TAX, META)

#Create and tree with sequences new ASV names  ----
sequences <- DNAStringSet(seqs)
alignment <- AlignSeqs(sequences, anchor=NA, verbose=T) # Create a phylo tree using de novo approach. The first step uses DECHIPER package to align sequences (for ~7000 ASV sequeneces, it took 11 minutes)
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA") # The second step, uses phagorns to create the phylogenetic tree with the DECHIPER output. Reading the alingment file 
dm <- dist.ml(phangAlign) #It took almost 10 mim, the tree is made with a distance based method, so a distance object is created with dist.ml
treeNJ <- NJ(dm) # Note, tip order != sequence order. It took for almost 20 min, construction of an unrooted tree usind NJ (neighbor joining)
class(treeNJ)

#### this step are for take a NJ tree as input and generate another tree using ML method (not for me at the moment, it took longer than 3 hours, it was stopped)
#fit = pml(treeNJ, data=phangAlign)
#fitGTR <- update(fit, k=4, inv=0.2) # fit a GTR+G+I (Generalized time-reversible with Gamma rate variation) maximum likelihood tree using the neighbor-joining tree as a starting point
#fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))
#fitGTR <- pml_bb(treeNJ, model="GTR+G(4)+I")

#Add to phyloseq object, tree and sequences ----
taxa_names(treeNJ) <- taxa_names(ASV)
setequal(taxa_names(ASV), treeNJ$tip.label)
ps <- merge_phyloseq(ps, treeNJ, sequences)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps))) #Change ASV names from dada2 format, to own names "ASV1 to ASVn"
head(taxa_names(ps))
saveRDS(ps, 'PQ00193_FFPE/phyloseq_obj/ps_decontam.rds') #For further use, to only load the phyloseq object
ps <- readRDS(".../phyloseq_obj/ps_decontam.rds") #For further use, to load de rds file created containing the previously phyloseq object created


# Filtering of samples and ASV regarding: reads, Bacteria, prevalence ----

# Phyloseq for filtering data considering previously table decontaminated with "Microdecon". The next step is used for 
# filtering out samples with low number of reads, ASV from mitochondria and chloroplast, and ASV prevalence of 1 ASV 
# with 1 read (low prevalence)

#Add reads per sample column to metadata
sample_data(ps)$total_reads <- sample_sums(ps) 
all_reads <- sample_data(ps) #metadata with total read number per sample

#Filter 1 out ASV unassigned and archaea at kingdom, to only remain bacteria ----
unique(tax_table(ps)[, "Kingdom"])
table(tax_table(ps)[, "Kingdom"], exclude = NULL)
ps_1 <- subset_taxa(ps, !is.na(Kingdom) & !Kingdom %in% c("Unassigned", "d__Archaea"))
table(tax_table(ps_1)[, "Kingdom"], exclude = NULL)

#Filter 2/3 out ASV in Phylum corresponding to chloroplast (order) and mitochondria (family) ----
sort(unique(tax_table(ps_1)[, "Order"]))
table(tax_table(ps_1)[, "Order"], exclude = NULL)
sort(unique(tax_table(ps_1)[, "Family"]))
table(tax_table(ps_1)[, "Family"], exclude = NULL)
head(tax_table(ps_1))
ps_2 <- subset_taxa(ps_1, (tax_table(ps_1)[,"Order"]!=" o__Chloroplast") | is.na(tax_table(ps_1)[,"Order"]))
ps_3 <- subset_taxa(ps_2, (tax_table(ps_2)[,"Family"]!=" f__Mitochondria") | is.na(tax_table(ps_2)[,"Family"]))

#Filter 4 out considering prevalence, if 1 ASV is present with only one read. In the table created, column 1 is the average of read counts
#and column 2 represents the sum of the reads. At phylum level, the ones presented with an average of 1 and 1 read is removed. 
#prev_ps = apply(X = otu_table(ps_3), MARGIN = ifelse(taxa_are_rows(ps_3), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
#prev_ps2 = data.frame(Prevalence = prev_ps, TotalAbundance = taxa_sums(ps_3), tax_table(ps_3))
#dfprev <- plyr::ddply(prev_ps2, "Phylum", function(prev_ps2){cbind(mean(prev_ps2$Prevalence),sum(prev_ps2$Prevalence))}) 
#prev_phyla <- c(" p__Latescibacterota", " p__Margulisbacteria",  " p__Methylomirabilota", " p__SAR324_clade(Marine_group_B)", 
                #" p__Thermotogota", " p__WPS-2") #phylum of median of reads = 1 and 1 read
#ps_4 <- subset_taxa(ps_3, !Phylum %in% prev_phyla)

#Filter 5 exclude phylum NA and change names  ----
unique(tax_table(ps_3)[, "Phylum"])
table(tax_table(ps_3)[, "Phylum"], exclude = NULL) #There is 39 NA for Phylum, we decide not remove them "Unknown_phylum"
ps_4 <- subset_taxa(ps_3, !is.na(Phylum))

tax_table(ps_4) <- gsub("d__", "", tax_table(ps_4)) 
tax_table(ps_4) <- gsub(" p__", "", tax_table(ps_4)) 
tax_table(ps_4) <- gsub(" c__", "", tax_table(ps_4))
tax_table(ps_4) <- gsub(" o__", "", tax_table(ps_4))
tax_table(ps_4) <- gsub(" f__", "", tax_table(ps_4))
tax_table(ps_4) <- gsub(" g__", "", tax_table(ps_4))
tax_table(ps_4) <- gsub(" s__", "", tax_table(ps_4))   
tax_table(ps_4)

#Remove samples from the same  patient and group, retain only 1 (the one with higher number of reads) ----
ps_5 <- subset_samples(ps_4, !sample.id %in% c("..", "..."))

saveRDS(ps_5, '.../phyloseq_obj/ps_decont_tax.rds') #For further use, to only load the phyloseq object
ps_5 <- readRDS(".../phyloseq_obj/ps_decont_tax.rds") #For further use, to load de rds file created containing the previously phyloseq object created
sample_data(ps_5)

#Filter 6/7, remove low frequency ASV. Remain ASV seen 1 time in at least 10% of the samples per pairwise or three group comparison. ----
#First, we subset data and then, remove lower frequency ASV and remove samples with less than 1000 reads. 

ps_ADSQ <- subset_samples(ps_5, subtype %in% c("AD", "SQ")) 
ps_fADSQ <- prune_samples(sample_sums(ps_ADSQ) > 1000, ps_ADSQ)
