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
load_('ggplot2', bioC=F)
load_('corrplot', bioC=F)
load_('vioplot', bioC=F)

#epeidi o lagani gamietai kai prepei na paroume ta xaraktiristika apo allou
my_data<-read.table("SraRunTable.txt",header = TRUE,sep="\t")

#project of interest
project_id <- 'SRP018008';

# Download the gene-level RangedSummarizedExperiment data
download_study(project_id)

# Load the object rse_gene
load(file.path(project_id, 'rse_gene.Rdata'))
class(rse_gene)

#extracting the count data
count_data <- assay(rse_gene)

#sto telos to characteristic vec2 8a exei "cancer" "no-cancer"
#stis katalliles 8eseis pou 8a eixe to kanoniko dataset an den
#itan malakismeno
characteristics_vec<-c()
for (i in 1:length(rse_gene$sample)){
	#for each column of rse gene
	name<-rse_gene$sample[i]
	for (j in 1:length(my_data$SRA_Sample_s)){
		#match the one in SRA table and keep the label
		name2<-my_data$SRA_Sample_s[j]
		if(name==name2){
			characteristics_vec[i]<-my_data$Sample_Name_s[j]
			break
		}
	}
}

#maybe change
characteristics_vec2<-c()
for(i in 1:length(characteristics_vec)){
	characteristics_vec2[i]<-as.character(my_data$Sample_Name_s[characteristics_vec[i]])
}
characteristics_vec2

#make binary
for(i in 1:length(characteristics_vec2)){
	if(grepl("Cancer",characteristics_vec2[i])){
		characteristics_vec2[i]<-"Cancer"
	}
	else{
		characteristics_vec2[i]<-"Normal"
	}
}
characteristics_vec2

#keep data with some significance
#remove the genes that have very small counts

#do normalization
#filter
toKeep <- apply(count_data, 1, sum) > 50 * dim(count_data)[2];
count_data <- count_data[toKeep, ];
dim(count_data)

#prwto link
#https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html

#object for edgeR
dgList <- DGEList(counts=count_data, genes=rownames(count_data),group=factor(characteristics_vec2))

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

png("test3.png")
dotplot(filteredEgo)
dev.off()
