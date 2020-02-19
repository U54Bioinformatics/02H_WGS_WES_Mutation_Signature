ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(data.table)
library("gridExtra")
library("MutationalPatterns")
library('NMF')
library(BSgenome.Hsapiens.UCSC.hg19)
library(RColorBrewer)
library(GenomicRanges)
library(GenomicFeatures)
library(BSgenome)
library(ref_genome, character.only = TRUE)
library(BSgenome.Celegans.UCSC.ce2)


Mutation_file<-fread("File.txt")
head(Mutation_file)
alfa = with(Mutation_file, GRanges(chr, IRanges(start=pos, end=pos)))
values(alfa) <- Mutation_file[,c("Sample","ref", "alt")]



types = mut_type(alfa)
context = mut_context(alfa, ref_genome)
type_context = type_context(alfa, ref_genome)


#### create GRange file in a list for each sample 
g<- list()
list<- unique(Mutation_file$Sample)

for(i in (1:length(list)))
{
  Mutation_file_single<- Mutation_file[Mutation_file$Sample== list[i],]
  alfa_single<- with(Mutation_file_single, GRanges(chr, IRanges(start=Mutation_file_single$pos, end=Mutation_file_single$pos), REF= ref, ALT=alt))
  names(alfa_single)<- Mutation_file_single$name
  genome(alfa_single) <- "hg19"
  g[[i]]<- (alfa_single)
}


names(g) <- list
type_occurrences <- mut_type_occurrences(g, ref_genome) ##### 



##plot 6 classes and CpGprevalence

p1 = plot_spectrum(type_occurrences, CT = TRUE, by= names(g))


mut_mat <- mut_matrix(vcf_list = g, ref_genome = ref_genome)

##plot 96 classes 
plot_96_profile(mut_mat, ymax = 0.06) 



##Run NMF for denovo signatures extraction

mut_mat = mut_mat + 0.0001
estimate = nmf(mut_mat, rank=2:6, method="brunet", nrun=1000, seed=123456)
plot(estimate)


## select best solution

nmf_res <- extract_signatures(mut_mat, rank = 4, nrun = 10) 



##Assign signature names

colnames(nmf_res$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D")
head(nmf_res$signatures)


##Plot the 96-profile of the signatures

plot_96_profile(nmf_res$signatures, ymax = 0.05)


rownames(nmf_res$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D")
pc1 <- plot_contribution(nmf_res$contribution, nmf_res$signature, mode = "relative",
                         palette =c(brewer.pal(8,"Dark2")), coord_flip = FALSE)



##plot the contribution of the signatures
pch1 <- plot_contribution_heatmap(nmf_res$contribution,sig_order = c("Signature A", "Signature B", "Signature C", "Signature D"), cluster_samples=FALSE)


#Upload 30 Signature COSMIC catalogue for fitting part
sp_url <- paste("http://cancer.sanger.ac.uk/cancergenome/assets/", "signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
cancer_signatures = cancer_signatures[as.vector(new_order),]
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
cancer_signatures = as.matrix(cancer_signatures[,4:33]) ###### all COSMIC Signatures(columns signatures, rows 96 classes)


df = cos_sim_matrix(nmf_res$signatures, cancer_signatures)
df1 <-plot_cosine_heatmap(df)



cos_sim_samples_signatures = cos_sim_matrix(mut_mat, cancer_signatures)
# plot heatmap
plot_cos_sim_samples_signatures1 = plot_cosine_heatmap(cos_sim_samples_signatures)
plot_cos_sim_samples_signatures1


cancer_signatures <- as.data.frame.matrix(cancer_signatures)
### include new signature in the fitting catalogue
cancer_signatures$MM1<- as.numeric(nmf_res$signatures[,4])/sum(nmf_res$signatures[,4])
mut_mat[,1:ncol(mut_mat)] = apply(mut_mat[,1:ncol(mut_mat)], 2, function(x) as.numeric(as.character(x)))
cancer_signatures[,1:ncol(cancer_signatures)] = apply(cancer_signatures[,1:ncol(cancer_signatures)], 2,function(x) as.numeric(as.character(x)))
fit_res <- fit_to_signatures(mut_mat, as.matrix(cancer_signatures))



# get relative contribution
contribution = t(fit_res$contribution)
# relative contribution
contribution_norm = contribution / rowSums(contribution)
# get maximum contribution over all samples per signature
max_contribution_norm = apply(contribution_norm, 2, function(x) max(x)) 
# get signatures with at least 10% contribution in at least 1 sample
select = which(max_contribution_norm > 0.1)


cancer_signatures_selected <- as.data.frame.matrix(cancer_signatures[,c('list of signatures with at least 10% contribution in at least 1 sample')])


cancer_signatures_selected$MM1<- as.numeric(nmf_res$signatures[,4])/sum(nmf_res$signatures[,4])
mut_mat[,1:ncol(mut_mat)] = apply(mut_mat[,1:ncol(mut_mat)], 2, function(x) as.numeric(as.character(x)))
cancer_signatures_selected[,1:ncol(cancer_signatures_selected)] = apply(cancer_signatures_selected[,1:ncol(cancer_signatures_selected)], 2,function(x) as.numeric(as.character(x)))
fit_res <- fit_to_signatures(mut_mat, as.matrix(cancer_signatures_selected))


plot_contribution(fit_res$contribution,cancer_signatures[,1:ncol(cancer_signatures)],coord_flip = FALSE, mode = "relative", palette = c(brewer.pal(12,"Paired"), brewer.pal(8,"Dark2"),brewer.pal(11,"Spectral")))

plot_contribution_heatmap(fit_res$contribution,cluster_samples = FALSE, method = "complete")





