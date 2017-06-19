#diale3i lagani
#http://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/rnaSeq_DE.pdf

#Our source of label vector:
#https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP018008&go=go
#https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP018008

#Define a function for loading packages
load_ = function(pkg, bioC=T) {
	#character.only has to be set to True in order for require() or library() to realize it's dealing with a variable
	if(!require(pkg, character.only=T, quietly = T)) {
		if(bioC){
			source(file = "http://bioconductor.org/biocLite.R")
			biocLite(pkg, dependencies=T)
		} else {
			install.packages(pkg)
		}
	}
	library(pkg, character.only=T)
}

#Installing and loading the required packages
load_("codetools", bioC=F)
load_("glue",bioC=F)
load_("doRNG", bioC=F)
load_("recount")
load_("DESeq2")
load_("edgeR")
load_("AnnotationDbi")
load_("org.Hs.eg.db")
load_("clusterProfiler")
load_('biomaRt')
load_('plotly', bioC=F)
load_('corrplot', bioC=F)

#fores pou 8a tre3ei to bootstrap
n<-1

#====================================================
#===================  Obtaining the Dataset   ==================
#====================================================

#project of interest
project_id <- 'SRP018008';

# Download the gene-level RangedSummarizedExperiment data
download_study(project_id)

# Load the object rse_gene
load(file.path(project_id, 'rse_gene.Rdata'))

#finding the labels we are missing
sra_run_table<-read.table("SraRunTable.csv",header = TRUE,sep="\t")
# We get the columns for sample id and sample classes (cancer / normal)
sra_run_subset <- sra_run_table[,c(11,12)]

# Using a regular expression we remove the prefix before the cancer/normal
# The first condition for the data are to be B123_Normal or B123_Cancer so we just erase the prefix

sra_run_subset[2]<-apply(sra_run_subset[2],2, function(x) sub(".*_","",x))

# On the other hand e.g B123-12N (for normal) or B123-12T (for tumor), we replace the whole value
sra_run_subset[2]<-apply(sra_run_subset[2],2, function(x) sub(".*-.*N","Normal",x))
sra_run_subset[2]<-apply(sra_run_subset[2],2, function(x) sub(".*-.*T","Cancer",x))

# Get the samples from the summary experiment and put them in a dataframe
transcripts<- as.data.frame(sort(rse_gene$sample))

sra_run_subset<- sra_run_subset[order(sra_run_subset$SRA_Sample_s),]

# Rename the column name
colnames(transcripts)[1]<-"RSA_Samples"

#extracting the count data
count_data <- assay(rse_gene)

#we have the label vector, we know need to sort it corresponding to the recount data
characteristics_vec<-c()
labels<-c()
for (i in 1:length(rse_gene$sample)){
	#for each column downloaded from recount
	name<-rse_gene$sample[i]
	for (j in 1:length(sra_run_subset$SRA_Sample_s)){
		#match the one in SRA table and keep the label
		name2<-sra_run_subset$SRA_Sample_s[j]
		if(name==name2){
			labels[i]<-sra_run_subset$Sample_Name_s[j]
		}
	}
}

labels

#====================================================
#====================   Quality Check  ======================
#====================================================

color = labels
color[color=='Cancer'] = 'red'
color[color=='Normal'] = 'blue'
png("BoxPlot1.png")
boxplot(count_data, xlab = "observations", ylab = "counts", col = color, outline = FALSE)
dev.off()

#Many counts are close to zero and provide skewed boxplots
#remove the genes that have very small counts

toKeep <- apply(count_data, 1, sum) > 50 * dim(count_data)[2];
count_data <- count_data[toKeep, ];
dim(count_data)

png("BoxPlot2.png")
boxplot(count_data, xlab = "observations", ylab = "counts", col = color, outline = FALSE)
dev.off()

#The boxplots in the far right end seem to give a distorted picture compared to the others
#these probably correspond to the samples that weren't in the original paper
count_data<-count_data[,1:75]
labels<-labels[1:75]

png("BoxPlot3.png")
boxplot(count_data, xlab = "observations", ylab = "counts", col = color, outline = FALSE)
dev.off()


#Normalization
Normalize <- function(x){(x-min(x))/(max(x)-min(x))}
count_data = Normalize(count_data)

#Visualize through Principal Component Analysis
PCs = as.data.frame(prcomp(count_data)$rotation)

library(plotly)

plot_pca <- plot_ly(data = PCs, x = ~PC1, y = ~PC2, color = labels, colors = 'Spectral')
plot_pca <- layout(plot_pca, title = "Bladder Cancer PCs",
       xaxis = list(title = "PC 1"),
       yaxis = list(title = "PC 2"))

plot_pca

png("CorrPlot.png")
Corre = cor(count_data)
corrplot(Corre, method='color')
dev.off()

#====================================================
#====================================================
#====================================================

#antwnis 19/6/17
categories_vec <- c()
for(i in 1:n){

#des posoi einai oi cancer kai oi normal
cancer_patients<-length(which(labels=="Cancer")) #53
normal_patients<-length(which(labels=="Normal")) #40

#krata ta indexes ston count data pou einai to ka8e group
cancer_col_indexes <- which(labels=="Cancer")
normal_col_indexes <- which(labels=="Normal")

#xwrise ta group meta3i tous se 2 mikrotera
normal_data<- count_data[,normal_col_indexes]
dim(normal_data)

cancer_data<- count_data[,cancer_col_indexes]
dim(cancer_data)

#kataskeuase nea indexes pou 8a xrisimopoih8oun me replacement
cancer_new_indexes<-sample(cancer_patients, cancer_patients, replace = TRUE, prob = NULL)
normal_new_indexes<-sample(normal_patients, normal_patients, replace = TRUE, prob = NULL)

#kataskeauase ta permutation tou ka8e group analoga me ta indexes
final_normal_data<-normal_data[,normal_new_indexes]
final_cancer_data<-cancer_data[,cancer_new_indexes]

#kane bind ta nea group se ena megalo
count_data<-cbind(final_normal_data,final_cancer_data)

#to charactericsd_vec2 einai pleon me tin seira ola ta normal mazi ola ta cancer mazi giati
#etsi ta kaname sto cbind
labels<-c(rep("Normal",normal_patients),rep("Cancer",cancer_patients))

#sbise ta colnames giati twra exoume epanotopo8etisi kai den sineragazetai to DGElist
colnames(count_data) <- NULL

#antwnis 19/6/17

#prwto link
#https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html

#object for edgeR

dgList <- DGEList(counts=count_data, group=factor(labels)) #genes=rownames(count_data),


png("test.png")
plotMDS(dgList, method="bcv", col=as.numeric(dgList$samples$group))
legend("bottomleft", as.character(unique(dgList$samples$group)), col=1:3, pch=20)
dev.off()

#design matrix: (one hot encoding)
design.mat <- model.matrix(~ 0 + dgList$samples$group)
colnames(design.mat) <- levels(dgList$samples$group)

#estimating dispersion #autes oi 3 entoles kanoun to idio pragma? need to explore further
d2 <- estimateGLMCommonDisp(dgList,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="auto")
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
png("test2.png")
plotBCV(d2)
dev.off()

#trito link
#Compute genewise exact tests for differences in the means between two groups of negative-binomially distributed counts.
et <- exactTest(d2)
results_edgeR <- topTags(et,adjust.method="BH", n = nrow(count_data), sort.by = "PValue")
results_edgeR <- topTags(et, n = nrow(count_data), sort.by = "PValue")

#krata mono ta data
edger_table<-results_edgeR$table
#apo auta mono ta differentially expressed genes
#correction p-value
diff_genes<-edger_table[edger_table$PValue<0.05/nrow(edger_table),]
#mono ta onomata twn genes
#sbise kati malakies pou exei me teleies kai tetoia
#it ruins the entrez id, so we get rid of the dots
diff_genes_names_<-diff_genes$genes
#sbise kati malakies pou exei me teleies kai tetoia
diff_genes_names<-gsub("\\..*","",diff_genes_names_)

universe_genes<-gsub("\\..*","",edger_table$genes)

Corre = cor(count_data[diff_genes_names_,])
corrplot(Corre)

###########################################################################

#Mapping from ensembl ID to entrez ID
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

G_list1 <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
"entrezgene", "description"),values=diff_genes_names,mart= mart)

G_list2 <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
"entrezgene", "description"),values=universe_genes,mart= mart)

entrez_diff_genes_names<-G_list1[,2]
entrez_diff_genes_names<-as.character(entrez_diff_genes_names[!is.na(entrez_diff_genes_names)])

universe_genes<-G_list2[,2]
universe_genes<-as.character(universe_genes[!is.na(universe_genes)])


ego <- enrichGO(gene          = entrez_diff_genes_names,
                universe      = universe_genes,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "fdr",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                minGSSize     = 15,
                maxGSSize     = 500,
                readable      = TRUE)
#

head(ego)[, 1:7]


filteredEgo <- gofilter(ego, level = 4);
dim(filteredEgo)

categories_vec<-c(categories_vec,filteredEgo$ID)
#print(categories_vec)
}


categories_freq<-table(categories_vec)

png("test4.png")
hist(categories_freq, breaks = 12, col = "lightblue", border = "pink")
dev.off()
