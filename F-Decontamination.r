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
ASV_qza <- read_qza("PQ00193_FFPE/FFPE_All/ASV/dada2table_All.qza")
ASV_matrix <- ASV_qza$data
ASV_matrix <- as.data.frame(ASV_matrix)

samples_metadata <- read.table("PQ00193_FFPE/manifest_files/FFPE_All.tsv", sep = '\t', header = T)
names(ASV_matrix) <- samples_metadata$feature.type.number[match(names(ASV_matrix), samples_metadata$sample.id)]

# Reorder file for microdecon: ASV_id, blanks/controls, Samples(ordered per group), taxa. ----
ASV_matrix <- select(ASV_matrix, "Blanc", "Control", 
                     "AD differentiated 1", "AD differentiated 2", "AD differentiated 3", "AD differentiated 4", "AD differentiated 5", "AD differentiated 6", "AD differentiated 7", "AD differentiated 8", "AD differentiated 9", "AD differentiated 10",
                     "AD differentiated 11", "AD differentiated 12", "AD differentiated 13", "AD differentiated 14", "AD differentiated 15","AD differentiated 16", "AD differentiated 17", "AD differentiated 18", "AD differentiated 19", "AD differentiated 20",
                     "AD differentiated 21", "AD differentiated 22", "AD differentiated 23", "AD differentiated 24", "AD differentiated 25","AD differentiated 26", "AD differentiated 27", "AD differentiated 28", "AD differentiated 29", "AD differentiated 30",
                     "AD undifferentiated 1", "AD undifferentiated 2", "AD undifferentiated 3","AD undifferentiated 4", "AD undifferentiated 5", "AD undifferentiated 6", "AD undifferentiated 7","AD undifferentiated 8", "AD undifferentiated 9", "AD undifferentiated 10", 
                     "AD undifferentiated 11","AD undifferentiated 12", "AD undifferentiated 13", "AD undifferentiated 14", "AD undifferentiated 15","AD undifferentiated 16",
                     "SQ differentiated 1", "SQ differentiated 2", "SQ differentiated 3", "SQ differentiated 4", "SQ differentiated 5", "SQ differentiated 6", "SQ differentiated 7", "SQ differentiated 8", "SQ differentiated 9", "SQ differentiated 10",
                     "SQ differentiated 11", "SQ differentiated 12", "SQ differentiated 13", "SQ differentiated 14", "SQ differentiated 15","SQ differentiated 16", "SQ differentiated 17", "SQ differentiated 18", "SQ differentiated 19", "SQ differentiated 20",
                     "SQ differentiated 21", "SQ differentiated 22", "SQ differentiated 23", "SQ differentiated 24", "SQ differentiated 25","SQ differentiated 26", "SQ differentiated 27", "SQ differentiated 28", "SQ differentiated 29", "SQ differentiated 30",
                     "SQ differentiated 31", "SQ differentiated 32", "SQ differentiated 33", "SQ differentiated 34", "SQ differentiated 35","SQ differentiated 36",
                     "SQ undifferentiated 1", "SQ undifferentiated 2", "SQ undifferentiated 3","SQ undifferentiated 4", "SQ undifferentiated 5", "SQ undifferentiated 6", "SQ undifferentiated 7","SQ undifferentiated 8", "SQ undifferentiated 9", "SQ undifferentiated 10", 
                     "SQ undifferentiated 11","SQ undifferentiated 12", "SQ undifferentiated 13", "SQ undifferentiated 14", "SQ undifferentiated 15","SQ undifferentiated 16", "SQ undifferentiated 17","SQ undifferentiated 18", "SQ undifferentiated 19","SQ undifferentiated 20",
                     "LCLC 1", "LCLC 2", "LCLC 3", 
                     "SCLC 1", "SCLC 2", "SCLC 3", "SCLC 4", "SCLC 5", "SCLC 6", "SCLC 7", "SCLC 8", "SCLC 9", "SCLC 10", "SCLC 11", "SCLC 12", "SCLC 13", "SCLC 14", "SCLC 15", "SCLC 16",
                     "MET 1", "MET 2", "MET 3", "MET 4", "MET 5", "MET 6", "MET 7", "MET 8", "MET 9", "MET 10", "MET 11", "MET 12", "MET 13", "MET 14", "MET 15", "MET 16", "MET 17",  
                     "Others 1", "Others 2", "Others 3", "Others 4", "Others 5", "Others 6", "Others 7", "Others 8", "Others 9",
                     "NORMAL 1", "NORMAL 2", "NORMAL 3", "NORMAL 4", "NORMAL 5", "NORMAL 6","NORMAL 7", "NORMAL 8")

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
ps <- readRDS("PQ00193_FFPE_LC/phyloseq_obj/ps_decontam.rds") #For further use, to load de rds file created containing the previously phyloseq object created


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
ps_5 <- subset_samples(ps_4, !sample.id %in% c("McAllister_19_25_A_A2", "McAllister_19_39_A", "McAllister_19_39_B", 
                                               "McAllister_19_8_A_A1", "McAllister_19_8_B_A5", "McAllister_20_1_A_A"))

saveRDS(ps_5, 'PQ00193_FFPE_LC/phyloseq_obj/ps_decont_tax.rds') #For further use, to only load the phyloseq object
ps_5 <- readRDS("PQ00193_FFPE_LC/phyloseq_obj/ps_decont_tax.rds") #For further use, to load de rds file created containing the previously phyloseq object created
sample_data(ps_5)

#Filter 6/7, remove low frequency ASV. Remain ASV seen 1 time in at least 10% of the samples per pairwise or three group comparison. ----
#First, we subset data and then, remove lower frequency ASV and remove samples with less than 1000 reads. 

ps_ADSQ <- subset_samples(ps_5, subtype %in% c("AD", "SQ")) 
ps_fADSQ <- prune_samples(sample_sums(ps_ADSQ) > 1000, ps_ADSQ) 
saveRDS(ps_fADSQ, 'PQ00193_FFPE_LC/phyloseq_obj/ps_fADSQ.rds')
with(sample_data(ps_fADSQ), table(subtype))

ps_ADsmall <- subset_samples(ps_5, subtype %in% c("AD", "SCLC")) 
ps_fADsmall <- prune_samples(sample_sums(ps_ADsmall) > 1000, ps_ADsmall) 
saveRDS(ps_fADsmall, 'PQ00193_FFPE_LC/phyloseq_obj/ps_fADsmall.rds')
with(sample_data(ps_fADsmall), table(subtype))

ps_SQsmall <- subset_samples(ps_5, subtype %in% c("SQ", "SCLC")) 
ps_fSQsmall <- prune_samples(sample_sums(ps_SQsmall) > 1000, ps_SQsmall) 
saveRDS(ps_fSQsmall, 'PQ00193_FFPE_LC/phyloseq_obj/ps_fSQsmall.rds')
with(sample_data(ps_fSQsmall), table(subtype))

ps_ADSQsmall <- subset_samples(ps_5, subtype %in% c("AD", "SQ", "SCLC")) 
ps_fADSQsmall <- prune_samples(sample_sums(ps_ADSQsmall) > 1000, ps_ADSQsmall) 
saveRDS(ps_fADSQsmall, 'PQ00193_FFPE_LC/phyloseq_obj/ps_fADSQsmall.rds')
with(sample_data(ps_fADSQsmall), table(subtype))

ps_ADd <- subset_samples(ps_5, feature.type %in% c("AD differentiated", "AD undifferentiated"))
ps_fADd <- prune_samples(sample_sums(ps_ADd) > 1000, ps_ADd) 
saveRDS(ps_fADd, 'PQ00193_FFPE_LC/phyloseq_obj/ps_ADd.rds')
with(sample_data(ps_fADd), table(feature.type))

ps_ADnor <- subset_samples(ps_5, subtype %in% c("AD", "NORMAL")) 
ps_fADnor <- prune_samples(sample_sums(ps_ADnor) > 1000, ps_ADnor) 
saveRDS(ps_fADnor, 'PQ00193_FFPE_LC/phyloseq_obj/ps_fADnor.rds')
with(sample_data(ps_fADnor), table(subtype))

ps_SQd <- subset_samples(ps_5, feature.type %in% c("SQ differentiated", "SQ undifferentiated"))
ps_fSQd <- prune_samples(sample_sums(ps_SQd) > 1000, ps_SQd) 
saveRDS(ps_fSQd, 'PQ00193_FFPE_LC/phyloseq_obj/ps_SQd.rds')
with(sample_data(ps_fSQd), table(feature.type))

# Venn diagrams ----

ps_ADSQ_small_nor <- subset_samples(ps_5, subtype %in% c("AD", "SQ", "SCLC", "NORMAL")) 
ps_fADSQ_small_nor <- prune_samples(sample_sums(ps_ADSQ_small_nor) > 1000, ps_ADSQ_small_nor) 
saveRDS(ps_fADSQ_small_nor, 'PQ00193_FFPE_LC/phyloseq_obj/ps_fADSQ_small_nor.rds')
ps_fADSQ_small_nor <- readRDS("PQ00193_FFPE_LC/phyloseq_obj/ps_fADSQ_small_nor.rds")
with(sample_data(ps_fADSQ_small_nor), table(subtype))
with(sample_data(ps_fADSQ_small_nor), table(feature.type))

pseq.rel <- microbiome::transform(ps_fADSQ_small_nor, "compositional")
group_subtype <- unique(as.character(meta(pseq.rel)$subtype))
group_feature.type <- unique(as.character(meta(pseq.rel)$feature.type))
print(group_subtype)
print(group_feature.type)

pseq.rel.f <- format_to_besthit(pseq.rel)
taxa_names(pseq.rel.f)[1:5]

## ASV more than 0 of abundance (presence/absence) and present in at leats 1 sample ----

list_core_subtype <- c() # an empty object to store information

for (n in group_subtype){ # for each variable n in DiseaseState
  ps.sub <- subset_samples(pseq.rel.f, subtype == n) # Choose sample from DiseaseState by n
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0, 
                         prevalence = 0) 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core_subtype[[n]] <- core_m
  } # add to a list core taxa for each group.

print(list_core_subtype)

list_core_feature.type <- c() # an empty object to store information

for (n in group_feature.type){ # for each variable n in DiseaseState 
  ps.sub <- subset_samples(pseq.rel.f, feature.type == n) # Choose sample from DiseaseState by n
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0, 
                         prevalence = 0)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core_feature.type[[n]] <- core_m} # add to a list core taxa for each group.

print(list_core_feature.type)

mycols_subtype <- c("AD" ="#FF0000", "SQ" ="#00008B", "SCLC" = "#008B8B", "NORMAL" = "#32CD32")
mycols_feature.type <- c("AD differentiated" ="#fc8d59", "AD undifferentiated" = "#fee090",
                         "SQ differentiated" = "#74add1", "SQ undifferentiated" = "#abd9e9", 
                         "SCLC" = "#008B8B", "NORMAL" = "#32CD32")

list_ASV_subtype_all <- do.call(qpcR:::cbind.na, list_core_subtype)
xlsx::write.xlsx(list_ASV_subtype_all, 'PQ00193_FFPE_LC/ASV_uniques.xlsx', sheetName = "subtype_allshared")

list_ASV_feature.type_all <- do.call(qpcR:::cbind.na, list_core_feature.type)
xlsx::write.xlsx(list_ASV_feature.type_all, 'PQ00193_FFPE_LC/ASV_uniques.xlsx', sheetName = "feature.type_allshared", append = TRUE)

unique_taxas_subtype <- unique_taxa(pseq.rel.f, treatment = "subtype")
list_ASV_subtype_uniques <- data.frame(lapply(unique_taxas_subtype, "length<-", max(lengths(unique_taxas_subtype))))
xlsx::write.xlsx(list_ASV_subtype_uniques, 'PQ00193_FFPE_LC/ASV_uniques.xlsx', sheetName = "subtype_uniques", append = TRUE)

unique_taxas_feature.type <- unique_taxa(pseq.rel.f, treatment = "feature.type")
list_ASV_feature.type_uniques <- data.frame(lapply(unique_taxas_feature.type, "length<-", max(lengths(unique_taxas_feature.type))))
xlsx::write.xlsx(list_ASV_feature.type_uniques, 'PQ00193_FFPE_LC/ASV_uniques.xlsx', sheetName = "feature.type_uniques", append = TRUE)

list_common <- common_taxa(pseq.rel.f, treatment = "subtype", subset = NULL, n = 'all')
xlsx::write.xlsx(list_common, 'PQ00193_FFPE_LC/ASV_uniques.xlsx', sheetName = "common_ASV", append = TRUE)

plot(venn(list_core_subtype), fills = mycols_subtype)
plot(venn(list_core_feature.type), fills = mycols_feature.type) #cant compute more than 5 groups

## ASV more than 5% of abundance and present in at leats a 50% of samples ----

list_core_subtype_restriction <- c() # an empty object to store information

for (n in group_subtype){ # for each variable n in DiseaseState
  ps.sub <- subset_samples(pseq.rel.f, subtype == n) # Choose sample from DiseaseState by n
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.05, # this could include the minimal abundance of the sample, if compositional, from 0 to 1 
                         prevalence = 0.1) #the presence in at least a 10% of samples 10/100
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core_subtype_restriction[[n]] <- core_m
} # add to a list core taxa for each group.

print(list_core_subtype_restriction)

list_ASV_subtype_restriction <- do.call(qpcR:::cbind.na, list_core_subtype_restriction)
xlsx::write.xlsx(list_ASV_subtype_restriction, 'PQ00193_FFPE_LC/ASV_uniques.xlsx', sheetName = "subtype_restriction", append = TRUE)

plot(venn(list_core_subtype_restriction), fills = mycols_subtype)

## Easier way to make venn diagram, however, cant extract ASV info for each group, thats why the previouly was performed ----
venn <- ps_venn(pseq.rel, "subtype", fill =c("AD" ="#FF0000", "NORMAL" = "#32CD32", "SCLC" = "#008B8B", "SQ" ="#00008B"))

## pie chart, with percentages of ASV shared with ----

ps_upset <- read.delim("PQ00193_FFPE_LC/ASV_percent.txt", sep = "\t", header = T)
ps_upset$Group <- factor(ps_upset$Group, levels = c("Normal", "SCLC", "AD", "SQ"))
ps_upset$Feature <- factor (ps_upset$Feature, levels = c("Unique AD", "Unique SQ", "Unique SCLC", "Shared in Lung cancer", "Normal in Lung cancer"))

detach(package:plyr)    
library(dplyr)

ps_upset2 <- ps_upset %>%                                    # Calculate percentage by group
  group_by(Group) %>%
  mutate(perc = Count / sum(Count)) %>% 
  as.data.frame()

colors_features <- c("Unique AD" = "#FF0000",
                     "Unique SQ" = "#00008B",
                     "Unique SCLC" = "#008B8B",
                     "Normal in Lung cancer" = "#32CD32",
                     "Shared in Lung cancer" = "gold")

ASV_plot_nor_SCLC_AD_SQ <- ggplot(ps_upset2, aes(x = 1, y = Count, fill = Feature)) +  
  geom_bar(stat="identity", position="fill", alpha = 0.5) +
  geom_text(aes(label = ifelse(percents >= 0, paste0(scales::percent(percents)), "")),
        stat = "identity", position = position_fill(vjust=0.5), size = 3.5) +
  theme_bw() +
  theme(aspect.ratio = 3.5, 
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.text.x = element_text(size = 12), 
        legend.position = "right",
        legend.key = element_blank(),
        panel.grid  = element_blank()) +
  facet_wrap(~Group, ncol = 4) +
  ylab("") +
  xlab("") +
  scale_fill_manual(values = colors_features)
  
ASV_plot_nor_SCLC_AD_SQ

ggsave(filename = "ASV_plot_nor_SCLC_AD_SQ.svg", path = "PQ00193_FFPE_LC/ASV_plot/", plot = ASV_plot_nor_SCLC_AD_SQ, device = "svg", height = 6, width = 6, dpi = 300)

ASV_plot_nor_AD <- ggplot(subset(ps_upset2, Group %in% c("Normal", "AD")), aes(x = 1, y = Count, fill = Feature)) +  
  geom_bar(stat="identity", position="fill", alpha = 0.5) +
  geom_text(aes(label = ifelse(perc >= 0, paste0(scales::percent(perc)), "")),
            stat = "identity", position = position_fill(vjust=0.5), size = 3) +
  theme_bw() +
  theme(aspect.ratio = 2.5, 
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.text.x = element_text(size = 9), 
        legend.text = element_text(size = 8),
        legend.position = "right",
        legend.key = element_blank(),
        panel.grid  = element_blank()) +
  facet_wrap(~Group, ncol = 2) +
  ylab("") +
  xlab("") +
  scale_fill_manual(values = colors_features)

ASV_plot_nor_AD

ggsave(filename = "ASV_plot_nor_AD.svg", path = "PQ00193_FFPE_LC/ASV_plot/", plot = ASV_plot_nor_AD, device = svg, height = 3.5, width = 3.5, dpi = 300)

ASV_plot_nor_SQ <- ggplot(subset(ps_upset2, Group %in% c("Normal", "SQ")), aes(x = 1, y = Count, fill = Feature)) +  
  geom_bar(stat="identity", position="fill", alpha = 0.5) +
  geom_text(aes(label = ifelse(percents >= 0, paste0(scales::percent(percents)), "")),
            stat = "identity", position = position_fill(vjust=0.5), size = 3.5) +
  theme_bw() +
  theme(aspect.ratio = 3.5, 
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.text.x = element_text(size = 12), 
        legend.position = "right",
        legend.key = element_blank(),
        panel.grid  = element_blank()) +
  facet_wrap(~Group, ncol = 4) +
  ylab("") +
  xlab("") +
  scale_fill_manual(values = colors_features)

ASV_plot_nor_SQ

ggsave(filename = "ASV_plot_nor_SQ.svg", path = "PQ00193_FFPE_LC/ASV_plot/", plot = ASV_plot_nor_SQ, device = "svg", height = 6, width = 6, dpi = 300)

ASV_plot_nor_SCLC <- ggplot(subset(ps_upset2, Group %in% c("Normal", "SCLC")), aes(x = 1, y = Count, fill = Feature)) +  
  geom_bar(stat="identity", position="fill", alpha = 0.5) +
  geom_text(aes(label = ifelse(percents >= 0, paste0(scales::percent(percents)), "")),
            stat = "identity", position = position_fill(vjust=0.5), size = 3.5) +
  theme_bw() +
  theme(aspect.ratio = 3.5, 
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.text.x = element_text(size = 12), 
        legend.position = "right",
        legend.key = element_blank(),
        panel.grid  = element_blank()) +
  facet_wrap(~Group, ncol = 4) +
  ylab("") +
  xlab("") +
  scale_fill_manual(values = colors_features)

ASV_plot_nor_SCLC

ggsave(filename = "ASV_plot_nor_SCLC.svg", path = "PQ00193_FFPE_LC/ASV_plot/", plot = ASV_plot_nor_SCLC, device = "svg", height = 6, width = 6, dpi = 300)

ASV_plot_nor_ADSQ <- ggplot(subset(ps_upset2, Group %in% c("Normal", "AD", "SQ")), aes(x = 1, y = Count, fill = Feature)) +  
  geom_bar(stat="identity", position="fill", alpha = 0.5) +
  geom_text(aes(label = ifelse(percents >= 0, paste0(scales::percent(percents)), "")),
            stat = "identity", position = position_fill(vjust=0.5), size = 3.5) +
  theme_bw() +
  theme(aspect.ratio = 3.5, 
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.text.x = element_text(size = 12), 
        legend.position = "right",
        legend.key = element_blank(),
        panel.grid  = element_blank()) +
  facet_wrap(~Group, ncol = 4) +
  ylab("") +
  xlab("") +
  scale_fill_manual(values = colors_features)

ASV_plot_nor_ADSQ

ggsave(filename = "ASV_plot_nor_ADSQ.svg", path = "PQ00193_FFPE_LC/ASV_plot/", plot = ASV_plot_nor_ADSQ, device = "svg", height = 6, width = 6, dpi = 300)


# Venn diagrams ----

ps_ADnor <- readRDS("PQ00193_FFPE_LC/phyloseq_obj/ps_fADnor.rds") 
ps_fADnor <- prune_samples(sample_sums(ps_ADnor) > 1000, ps_ADnor) 
with(sample_data(ps_fADnor), table(subtype))

pseq.rel <- microbiome::transform(ps_fADnor, "compositional")
group_subtype <- unique(as.character(meta(pseq.rel)$subtype))
print(group_subtype)

pseq.rel.f <- format_to_besthit(pseq.rel)
taxa_names(pseq.rel.f)[1:5]

## ASV more than 0 of abundance (presence/absence) and present in at leats 1 sample ----

list_core_subtype <- c() # an empty object to store information

for (n in group_subtype){ # for each variable n in DiseaseState
  ps.sub <- subset_samples(pseq.rel.f, subtype == n) # Choose sample from DiseaseState by n
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0, 
                         prevalence = 0) 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core_subtype[[n]] <- core_m
} # add to a list core taxa for each group.

print(list_core_subtype)

mycols_subtype <- c("AD" ="#FF0000", "NORMAL" = "#32CD32")

list_ASV_subtype_all <- do.call(qpcR:::cbind.na, list_core_subtype)
xlsx::write.xlsx(list_ASV_subtype_all, 'PQ00193_FFPE_LC/ASV_uniques.xlsx', sheetName = "subtype_allshared")

list_ASV_feature.type_all <- do.call(qpcR:::cbind.na, list_core_feature.type)
xlsx::write.xlsx(list_ASV_feature.type_all, 'PQ00193_FFPE_LC/ASV_uniques.xlsx', sheetName = "feature.type_allshared", append = TRUE)

unique_taxas_subtype <- unique_taxa(pseq.rel.f, treatment = "subtype")
list_ASV_subtype_uniques <- data.frame(lapply(unique_taxas_subtype, "length<-", max(lengths(unique_taxas_subtype))))
xlsx::write.xlsx(list_ASV_subtype_uniques, 'PQ00193_FFPE_LC/ASV_uniques.xlsx', sheetName = "subtype_uniques", append = TRUE)

unique_taxas_feature.type <- unique_taxa(pseq.rel.f, treatment = "feature.type")
list_ASV_feature.type_uniques <- data.frame(lapply(unique_taxas_feature.type, "length<-", max(lengths(unique_taxas_feature.type))))
xlsx::write.xlsx(list_ASV_feature.type_uniques, 'PQ00193_FFPE_LC/ASV_uniques.xlsx', sheetName = "feature.type_uniques", append = TRUE)

list_common <- common_taxa(pseq.rel.f, treatment = "subtype", subset = NULL, n = 'all')
xlsx::write.xlsx(list_common, 'PQ00193_FFPE_LC/ASV_uniques.xlsx', sheetName = "common_ASV", append = TRUE)

plot(venn(list_core_subtype), fills = mycols_subtype)


## ASV more than 5% of abundance and present in at leats a 10% of samples ----

list_core_subtype_restriction <- c() # an empty object to store information

for (n in group_subtype){ # for each variable n in DiseaseState
  ps.sub <- subset_samples(pseq.rel.f, subtype == n) # Choose sample from DiseaseState by n
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.05, # this could include the minimal abundance of the sample, if compositional, from 0 to 1 
                         prevalence = 0.1) #the presence in at least a 10% of samples 10/100
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core_subtype_restriction[[n]] <- core_m
} # add to a list core taxa for each group.

print(list_core_subtype_restriction)

list_ASV_subtype_restriction <- do.call(qpcR:::cbind.na, list_core_subtype_restriction)
xlsx::write.xlsx(list_ASV_subtype_restriction, 'PQ00193_FFPE_LC/ASV_uniques.xlsx', sheetName = "subtype_restriction", append = TRUE)

plot(venn(list_core_subtype_restriction), fills = mycols_subtype) 
