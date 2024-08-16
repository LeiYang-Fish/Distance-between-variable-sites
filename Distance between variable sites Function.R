######################################## If you want to show the results by gene
library(seqinr)
library(Biostrings)

# Define the function DIST_Var, which can return a list of distances among neighbor variable sites.
DIST_Var <- function(x){
  Both_paralogs <- read.fasta(x,seqtype = "DNA", as.string = TRUE,forceDNAtolower = FALSE,set.attributes = FALSE)
  S1 <- unlist(strsplit(as.character(Both_paralogs[1]), ""))
  S2 <- unlist(strsplit(as.character(Both_paralogs[2]), ""))
  Paralog_diff <- which(S1!=S2)
  dist_var_sites <- sapply(2:length(Paralog_diff), function(x) Paralog_diff[x]-Paralog_diff[x-1]); dist_var_sites
}

# Go to the folder for RAG1. Each species has one .txt file that contains its two gene copy sequences.
List_compile<- c()
for (filename in list.files(pattern=".txt")) {
  List_compile <- c(List_compile, DIST_Var(filename))
}

write.csv(List_compile,file = "dist_var_lists_RAG1.csv", row.names = FALSE)
# Repeat above for EGR2B, EGR3, IRBP2, and RAG2 and combine the results into a single .csv file. One column for each gene.

# Plot the results
Paralog_distances <- read.csv("dist_var_lists_by_gene_All_5_genes.csv")
stripchart(Paralog_distances,xlab="Gene",ylab="Distance between neighboring variable sites (bp)",
           group.names=c("RAG1","EGR2B","EGR3","IRBP2","RAG2"),col="blue",vertical=TRUE,pch=1)




######################################## If you want to show the results by gene and species
library(seqinr)
library(Biostrings)

# Define the function DIST_Var
DIST_Var <- function(x){
  Both_paralogs <- read.fasta(x,seqtype = "DNA", as.string = TRUE,forceDNAtolower = FALSE,set.attributes = FALSE)
    S1 <- unlist(strsplit(as.character(Both_paralogs[1]), ""))
    S2 <- unlist(strsplit(as.character(Both_paralogs[2]), ""))
  Paralog_diff <- which(S1!=S2)
  dist_var_sites <- sapply(2:length(Paralog_diff), function(x) Paralog_diff[x]-Paralog_diff[x-1]); dist_var_sites
  }

# Define the function cbind.fill to fill NA to shorter columns
cbind.fill<-function(...){
  nm <- list(...) 
  nm<-lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x)))))
}

# Bind the results for all species
List_compile<- list()
file_list <- list()
for (filename in list.files(pattern=".txt")) {
  file_list <- append(file_list, filename)
  List_compile <- cbind.fill(List_compile, DIST_Var(filename))
  }

# To get the abbreviations of species names
file_list
Genus_list <- gsub("_.*","", file_list)
Species_list <- gsub("[:.:].*","",gsub(".*_","", file_list))
Name_abrev <- lapply(1:length(file_list), function (x) 
  paste0(substring(as.character(Genus_list[x]),1,2),".",substring(as.character(Species_list[x]),1,2)))

List_compile <- List_compile[,-1] # remove the first column which contains only NAs

write.csv(List_compile,file = "dist_var_lists.csv", row.names = FALSE) # write as CSV to check and also remove the row numbers

# Plot the results
Paralog_distances <- read.csv("dist_var_lists.csv")
stripchart(Paralog_distances,xlab="Species with cloned RAG1 copy identified two alleles",ylab="Distance between neighboring variable sites (bp)",
           group.names=Name_abrev,col="blue",vertical=TRUE,pch=1)



