## run-sleuth.R by Rohan Maddamsetti.

## TODO: ask whether genes with bad start codons have little to no expression.

library(tidyverse)
library(sleuth)
library(cowplot)
library(RColorBrewer)
library(VennDiagram)
library(ggrepel)

## q-value threshold as a global value:
THRESHOLD <- 0.01

## helper function for running sleuth analysis.
make.transcript.table <- function(proj.dir, analysis.dir, s) {
    data <- read.table(
        file.path(proj.dir, "results", analysis.dir, s, "abundance.tsv"),
        header = TRUE, sep="\t", stringsAsFactors=FALSE) %>%
        select(target_id, tpm) %>%
        mutate(sample=s)
    return(data)
}

## helper to measure the overlap between the significant mRNA sets.
dice.similarity <- function(x,y) {
    x <- unique(x)
    y <- unique(y)
    return(2*length(intersect(x,y))/(length(x)+length(y)))
}

## helper to measure the overlap between the significant mRNA sets.
jaccard.index <- function(x,y) {
    x <- unique(x)
    y <- unique(y)
    return(length(intersect(x,y))/length(union(x,y)))
}

## Make a PCA plot of either single strains, or all three strains together.
## I hacked sleuth::plot_pca internals to make the plot look the way I want.
my_plot_pca <- function(obj, pc_x = 1L, pc_y = 2L, use_filtered = TRUE, units = "est_counts", point_size = 3, point_alpha = 0.8, my.color=NA, ...) {
    stopifnot(is(obj, "sleuth"))
    stopifnot(sleuth:::check_norm_status(obj))
    units <- sleuth:::check_quant_mode(obj, units)
    mat <- NULL
    if (use_filtered) {
        mat <- sleuth:::spread_abundance_by(obj$obs_norm_filt, units, 
            obj$sample_to_covariates$sample)
    } else {
        mat <- sleuth:::spread_abundance_by(obj$obs_norm, units, obj$sample_to_covariates$sample)
    }
    pca_res <- prcomp(t(mat))
    pcs <- sleuth:::as_df(pca_res$x[, c(pc_x, pc_y)])
    pcs$sample <- rownames(pcs)
    rownames(pcs) <- NULL
    pc_x <- paste0("PC", pc_x)
    pc_y <- paste0("PC", pc_y)
    pcs <- dplyr::left_join(pcs, obj$sample_to_covariates, by = "sample")
    
    p <- ggplot(pcs,aes_string(pc_x,pc_y, colour="Clone",shape="Treatment")) +
        geom_point(size = point_size, alpha = point_alpha) +
        xlab("principal component 1") +
        ylab("principal component 2")

    if (my.color %in% c('blue','purple','red')) {
        p <- p + scale_colour_manual(values=c(my.color))
    }
    else {
        p <- p + scale_colour_manual(values=c('blue','purple','red'))
    }
    return(p)
}

## project directory is the parent of the src directory.
proj.dir <- file.path(dirname(getwd()))

## get annotation parsed from the headers of the transcriptome fasta
## file by parse_S288c_headers.py. Note that separation char is '|'.
annot.file <- file.path(proj.dir,"results","transcript_annotation.csv")
rna.annotation <- read.csv(annot.file, sep='|',
                           header=TRUE, as.is = TRUE) %>% rename(target_id=mRNA) %>%
    mutate(gene_symbol=ifelse(is.na(gene_symbol),gene,gene_symbol))

## metadata showing experimental design (Genotype + Treatment), with all three strains.
metadata <- read.csv(
    file.path(proj.dir, "data", "sleuth-metadata.csv"),
    header=TRUE, stringsAsFactors=FALSE)

#############################################################
## analysis of single genotypes, and their response to differing conditions.

run.clone.sleuth <- function(proj.dir, clone) {
    stopifnot(clone %in% c('HMY12', 'HMY127', 'HMY362'))
    
    ## generate analysis.dir from the clone name.
    analysis.dir <- paste0(clone,"-ref-kallisto-output")
    
    ## assert that the analysis.dir is valid.
    stopifnot(analysis.dir %in% c("S288c-ref-kallisto-output",
                                  "HMY12-ref-kallisto-output",
                                  "HMY127-ref-kallisto-output",
                                  "HMY362-ref-kallisto-output"))
    if (analysis.dir == "S288c-ref-kallisto-output") {
        print("S288c reference transcriptome-based analysis")
    } else {
        print(paste(clone, "reference transcriptome-based analysis"))
    }
    
    ## filter metadata on just the clone of interest.
    metadata <- metadata %>% filter(Clone == clone)
    
    ## add path names of kallisto output directories to the metadata table.
    metadata <- metadata %>%
        mutate(path = file.path(proj.dir,"results", analysis.dir, sample, "abundance.h5"))
    
    ## use partial function application on make.transcript.table.
    make.transcript.tb <- partial(make.transcript.table,proj.dir, analysis.dir)
    
    ## make a df of inferred transcript abundances for each sample.  
    transcript.df <- map(metadata$sample, make.transcript.tb) %>%
        reduce(left_join)
    
    my.annotation <- transcript.df %>% select(target_id) %>% distinct()
    
    so <- sleuth_prep(metadata, target_mapping = my.annotation,
                      extra_bootstrap_summary = TRUE) %>%
        sleuth_fit(~ Treatment, 'E') %>%
        sleuth_fit( ~1, 'intercept') %>%
        ## test contribution of environment in predicting expression.
        sleuth_lrt('intercept','E')
    return(so)
}

HMY12.clone.so <- run.clone.sleuth(proj.dir, "HMY12")
HMY127.clone.so <- run.clone.sleuth(proj.dir, "HMY127")
HMY362.clone.so <- run.clone.sleuth(proj.dir, "HMY362")

### Figure 2.
### A) PCA plots
### B) variance per principal component.
### C) Venn Diagram.
### Put panel C together with AB in 'post' with Illustrator.

HMY12.pca <- my_plot_pca(HMY12.clone.so, my.color='blue') + theme_classic() + ggtitle("HMY12") + xlim(-60000,60000) + ylim(-60000,60000) + guides(color=FALSE, shape=FALSE)
HMY127.pca <- my_plot_pca(HMY127.clone.so, my.color='purple') + theme_classic() + ggtitle("HMY127") + xlim(-60000,60000) + ylim(-60000,60000) + guides(color=FALSE, shape=FALSE) 
HMY362.pca <- my_plot_pca(HMY362.clone.so, my.color='red') + theme_classic() + ggtitle("HMY362") + xlim(-60000,60000) + ylim(-60000,60000) + guides(color=FALSE, shape=FALSE) 
Fig2A <- plot_grid(HMY12.pca, HMY127.pca, HMY362.pca, nrow=3)
## plot variance per principal component.
HMY12.pca.var.plot <- plot_pc_variance(HMY12.clone.so) + theme_classic() + ggtitle("HMY12")
HMY127.pca.var.plot <- plot_pc_variance(HMY127.clone.so) + theme_classic() + ggtitle("HMY127")
HMY362.pca.var.plot <- plot_pc_variance(HMY362.clone.so) + theme_classic() + ggtitle("HMY362")
Fig2B <- plot_grid(HMY12.pca.var.plot, HMY127.pca.var.plot, HMY362.pca.var.plot, nrow=3)

Fig2AB <- plot_grid(Fig2A, Fig2B,labels=c('A','B'),nrow=1)
ggsave("../results/figures/Fig2AB.pdf",Fig2AB,width=7,height=7)

### Significant genes.
HMY12.test.results <- sleuth_results(
  HMY12.clone.so, 'intercept:E', 'lrt', show_all=FALSE) %>%
dplyr::filter(qval <= THRESHOLD)

## quick check to see where FLO11 occurs.
all.HMY12.test.results <- sleuth_results(
  HMY12.clone.so, 'intercept:E', 'lrt', show_all=TRUE)
FLO11.HMY12 <- all.HMY12.test.results %>% left_join(rna.annotation) %>%
    filter(gene=='YIR019C')

HMY127.test.results <- sleuth_results(
  HMY127.clone.so, 'intercept:E', 'lrt', show_all=FALSE) %>%
dplyr::filter(qval <= THRESHOLD)

HMY362.test.results <- sleuth_results(
  HMY362.clone.so, 'intercept:E', 'lrt', show_all=FALSE) %>%
dplyr::filter(qval <= THRESHOLD)

annotated.HMY12.test.results <- left_join(HMY12.test.results, rna.annotation)
annotated.HMY127.test.results <- left_join(HMY127.test.results, rna.annotation)
annotated.HMY362.test.results <- left_join(HMY362.test.results, rna.annotation)

## write results to file.
write.csv(annotated.HMY12.test.results,file=file.path(proj.dir,"results","HMY12_LRT.csv"))
write.csv(annotated.HMY127.test.results,file=file.path(proj.dir,"results","HMY127_LRT.csv"))
write.csv(annotated.HMY362.test.results,file=file.path(proj.dir,"results","HMY362_LRT.csv"))

## look at what is in common in these gene sets.
HMY12.genes <- annotated.HMY12.test.results$gene
HMY127.genes <- annotated.HMY127.test.results$gene
HMY362.genes <- annotated.HMY362.test.results$gene

## find genes that are DE in all three strains.
gene.intrsctn <- Reduce(intersect,list(HMY12.genes,HMY127.genes,HMY362.genes))

dice.similarity(HMY12.genes,HMY127.genes)
jaccard.index(HMY12.genes,HMY127.genes)

dice.similarity(HMY12.genes,HMY362.genes)
jaccard.index(HMY12.genes,HMY362.genes)

dice.similarity(HMY127.genes,HMY362.genes)
jaccard.index(HMY127.genes,HMY362.genes)

## let's do a randomization test to see if the cardinality of the tf.intrsctn is
## significant or not.
all.genes <- rna.annotation$gene

calc.intersection.randomization.pval <- function(gene_vec, x, c1,c2,c3,reps=10000) {
    calc.intersection.size <- function(i) {
        v1 <- sample(gene_vec,size=c1)
        v2 <- sample(gene_vec,size=c2)
        v3 <- sample(gene_vec,size=c3)
        return(length(Reduce(intersect,list(v1,v2,v3))))
    }
    greater.than.threshold <- function(y) ifelse(y > x, 1, 0)
    
    rando.intersection.sizes <- map(seq_len(10000),.f=calc.intersection.size)
    in.tail <- map(rando.intersection.sizes,.f=greater.than.threshold)
    return(Reduce(sum,in.tail)/length(in.tail))
}

## overlap in gene sets is highly significant:
## p < 0.00001
calc.intersection.randomization.pval(all.genes,length(gene.intrsctn),
                                     length(HMY12.genes),
                                     length(HMY127.genes),
                                     length(HMY362.genes))

## Figure 2C: Venn Diagram of these numbers above.
Fig2C <- venn.diagram(x=list(HMY12.genes, HMY127.genes, HMY362.genes),
                      category.names = c(paste0("HMY12 ",
                                                '(',length(HMY12.genes),')'),
                                         paste0("HMY127 "
                                               ,'(',length(HMY127.genes),')'),
                                         paste0("HMY362 ",
                                                '(',length(HMY362.genes),')')),
                      filename = NULL,
                      output = TRUE ,
                      imagetype="svg" ,
                      height = 480 , 
                      width = 480 , 
                      resolution = 300,
                      compression = "lzw",
                      lwd = 1,
                      col=c("blue", 'purple', 'red'),
                      fill = c(alpha("blue",0.3), alpha('purple',0.3), alpha('red',0.3)),
                      cex=2,
                      fontfamily = "sans",
                      cat.cex=1.5,
                      cat.default.pos = "outer",
                      cat.pos = c(-27, 27, 135),
                      cat.dist = c(0.055, 0.055, 0.085),
                      cat.fontfamily = "sans",
                      cat.col = c("blue", 'purple', 'red'),
                      rotation = 1
                      )

## Put panels AB and C together using Illustrator.
ggsave(Fig2C,file=file.path(proj.dir,"results","figures","Fig2C.pdf"))

## do a rank correlation test between strains.
## Need equal length vectors for these tests:
## filter on the intersection of significantly DE genes.
HMY12.vec <- HMY12.genes[HMY12.genes %in% gene.intrsctn]
HMY127.vec <- HMY127.genes[HMY127.genes %in% gene.intrsctn]
HMY362.vec <- HMY362.genes[HMY362.genes %in% gene.intrsctn]

## get the rank of each item in each vec in gene.intrsctn,
## in order to do a rank correlation test.
HMY12.rank.vec <- match(HMY12.vec, gene.intrsctn)
HMY127.rank.vec <- match(HMY127.vec, gene.intrsctn)
HMY362.rank.vec <- match(HMY362.vec, gene.intrsctn)

## highly similar response in all pairs of three clones,
## although HMY127 and HMY362 signal is weaker (but still highly significant).
## z = 14.9, p < 10^-15.
cor.test(HMY12.rank.vec, HMY127.rank.vec, method="kendall")
## z = 17.7, p < 10^-15.
cor.test(HMY12.rank.vec, HMY362.rank.vec, method="kendall")
## z = 5.3, p < 10^-7.
cor.test(HMY127.rank.vec, HMY362.rank.vec, method="kendall")


####################################################################
## Make a spreadsheet of all 7 disjoint subsets of genes that are differentially
## expressed. Are these subsets enriched for any particular networks?

## encode these subsets using bit encoding.
## first bit == 1 if in HMY12
## second bit == 1 if in HMY127
## third bit == 1 if in HMY362

## now let's calculate the 7 disjoint sets in the Venn Diagram
## for what's DE in the three strains.

## Let's examine the intersection of all three.
##annotated.HMY12.111.results <- annotated.HMY12.test.results %>%
##    filter(gene %in% gene.intrsctn)
##write.csv(annotated.HMY12.111.results,file=file.path(proj.dir,"results","HMY12_111_LRT.csv"))

## 110
only.HMY12.and.127 <- HMY12.genes[(HMY12.genes %in% HMY127.genes) &
                                  !(HMY12.genes %in% gene.intrsctn)]
## 101
only.HMY12.and.362 <- HMY12.genes[(HMY12.genes %in% HMY362.genes) &
                                 !(HMY12.genes %in% gene.intrsctn)]
## 011
only.HMY127.and.362 <- HMY12.genes[(HMY127.genes %in% HMY362.genes) &
                                         !(HMY127.genes %in% gene.intrsctn)]
## 100
only.HMY12 <- HMY12.genes[!(HMY12.genes %in% HMY127.genes) &
                                !(HMY12.genes %in% HMY362.genes)]
## 010
only.HMY127 <- HMY127.genes[!(HMY127.genes %in% HMY12.genes) &
                                !(HMY127.genes %in% HMY362.genes)]
## 001
only.HMY362 <- HMY362.genes[!(HMY362.genes %in% HMY12.genes) &
                                !(HMY362.genes %in% HMY127.genes)]

DE.sets.df <- rna.annotation %>%
    mutate(in.111 = gene %in% gene.intrsctn) %>%
    mutate(in.110 = gene %in% only.HMY12.and.127) %>%
    mutate(in.101 = gene %in% only.HMY12.and.362) %>%
    mutate(in.011 = gene %in% only.HMY127.and.362) %>%
    mutate(in.100 = gene %in% only.HMY12) %>%
    mutate(in.010 = gene %in% only.HMY127) %>%
    mutate(in.001 = gene %in% only.HMY362)

write.csv(DE.sets.df,file=file.path(proj.dir,"results","DE_sets.csv"))

## examine the intersection in STRING-DB.
DE.111 <- filter(DE.sets.df,in.111==TRUE)
write.csv(DE.111,"../results/DE_111.csv")

## STRING only takes up to 2000 proteins...
## let's filter.
annotated.HMY12.111.results <- annotated.HMY12.test.results %>%
    filter(gene %in% gene.intrsctn) %>%
    arrange(qval) %>% top_n(2000)

write.csv(annotated.HMY12.111.results,"../results/DE_111_HMY12_2000.csv")

####################################################################
## joint analysis of all three strains, combined, over both conditions.

run.full.sleuth <- function(proj.dir, analysis.dir) {

  ## assert that the analysis.dir is valid.
  stopifnot(analysis.dir %in% c("S288c-ref-kallisto-output",
                                "HMY12-ref-kallisto-output",
                                "HMY127-ref-kallisto-output",
                                "HMY362-ref-kallisto-output"))
  if (analysis.dir == "S288c-ref-kallisto-output") {
    print("S288c reference transcriptome-based analysis")
  } else {
    print("GATK reference transcriptome-based analysis")
  }

  ## add path names of kallisto output directories to the metadata table.
  metadata <- metadata %>%
  mutate(path = file.path(proj.dir,"results", analysis.dir, sample, "abundance.h5"))

  ## use partial function application on make.transcript.table.
  make.transcript.tb <- partial(make.transcript.table, proj.dir, analysis.dir)
  
  ## make a df of inferred transcript abundances for each sample.
  transcript.df <- map(metadata$sample, make.transcript.tb) %>%
  reduce(left_join)
  
  my.annotation <- transcript.df %>% select(target_id) %>% distinct()
  
  so <- sleuth_prep(metadata, target_mapping = my.annotation,
                    extra_bootstrap_summary = TRUE) %>%
                    ## Model with GxE interaction.
                    sleuth_fit(~ Clone * Treatment, 'G*E') %>%
                    ## Model without GxE interaction.
                    sleuth_fit(~ Clone + Treatment, 'G+E') %>%
                    sleuth_fit(~ Treatment, 'E') %>%
                    sleuth_fit(~ Clone, 'G') %>%
                    ## test significance of GxE interaction in predicting expression.
                    sleuth_lrt('G+E','G*E') %>%
                    ## test contribution of environment in predicting expression.
                    sleuth_lrt('G','G+E') %>%
                    ## test contribution of genotype in predicting expression.
                    sleuth_lrt('E','G+E')
  return(so)
}

## Harold Pimentel suggests just doing the LRT test, as he thinks the Wald test
## may give too many false positives. If one wants to calculate log-fold-change
## in expression, one can get this from the kallisto normalized values while
## omitting the filtered out values like so:
## SNM <- kallisto_table(so, use_filtered = TRUE, normalized = TRUE, include_covariates = TRUE
## See discussion at: https://www.biostars.org/p/226022/#226190

S288c.ref.so <- run.full.sleuth(proj.dir, "S288c-ref-kallisto-output")
HMY12.ref.so <- run.full.sleuth(proj.dir, "HMY12-ref-kallisto-output")
HMY127.ref.so <- run.full.sleuth(proj.dir, "HMY127-ref-kallisto-output")
HMY362.ref.so <- run.full.sleuth(proj.dir, "HMY362-ref-kallisto-output")

## Figure 3: PCA plot and variance, for all three genotypes analyzed jointly.
Fig3A <- my_plot_pca(S288c.ref.so,units="est_counts") + theme_classic() + guides(color=FALSE, shape=FALSE)
Fig3B <- plot_pc_variance(S288c.ref.so) + theme_classic()
Fig3AB <- plot_grid(Fig3A,Fig3B,labels=c('A','B'),nrow=1)
ggsave("../results/figures/Fig3AB.pdf",Fig3AB)

## double check that labels are correct by using the default sleuth::plot_pca function.
full.pca.plot <- plot_pca(S288c.ref.so,text_labels=TRUE,color_by="Clone") + theme_classic() + xlim(-100000,100000) + ylim(-100000,100000) + guides(color=FALSE, shape=FALSE)
ggsave(filename="../results/figures/full_PCA_plot.pdf", full.pca.plot,height=4,width=6)

## plot PCA using tpm as unit, to show how the PCA depends on scale.
S1Fig <- my_plot_pca(S288c.ref.so,units="tpm") + theme_classic() + guides(color=FALSE,shape=FALSE)
ggsave(filename="../results/figures/S1Fig.pdf", S1Fig, height=3, width=3.5)

## Let's analyze the full sleuth results.
## Three LRT Tests were conducted.
## Here they are using S288c as the reference transcriptome.
S288c.GxE.test.results <- sleuth_results(S288c.ref.so, 'G+E:G*E', 'lrt', show_all = FALSE) %>% dplyr::filter(qval <= THRESHOLD)
S288c.E.test.results <- sleuth_results(S288c.ref.so, 'G:G+E', 'lrt', show_all = FALSE) %>% dplyr::filter(qval <= THRESHOLD)
S288c.G.test.results <- sleuth_results(S288c.ref.so, 'E:G+E', 'lrt', show_all = FALSE) %>% dplyr::filter(qval <= THRESHOLD)

## check robustness in the results doing the same test with the other transcriptomes.
filter.lrt.test.results <- function(S288c.so, HMY12.so, HMY127.so, HMY362.so, test.string='G+E:G*E',threshold=THRESHOLD) {
  S288c.test.results <- sleuth_results(S288c.so, test.string, 'lrt', show_all = FALSE) %>% dplyr::filter(qval <= threshold)
  HMY12.test.results <- sleuth_results(HMY12.so, test.string, 'lrt', show_all = FALSE) %>% dplyr::filter(qval <= threshold)
  HMY127.test.results <- sleuth_results(HMY127.so, test.string, 'lrt', show_all = FALSE) %>% dplyr::filter(qval <= threshold)
  HMY362.test.results <- sleuth_results(HMY362.so, test.string, 'lrt', show_all = FALSE) %>% dplyr::filter(qval <= threshold)

  
  filtered.test.results <- S288c.test.results %>%
  filter(target_id %in% HMY12.test.results$target_id) %>%
  filter(target_id %in% HMY127.test.results$target_id) %>%
  filter(target_id %in% HMY362.test.results$target_id)

  passed <- length(filtered.test.results$target_id)/length(S288c.test.results$target_id)
  print(passed)
  print("percent passed filtering")
  return(filtered.test.results)
}

filtered.GxE.test.results <- filter.lrt.test.results(S288c.ref.so, HMY12.ref.so, HMY127.ref.so, HMY362.ref.so, test.string='G+E:G*E')
filtered.E.test.results <- filter.lrt.test.results(S288c.ref.so, HMY12.ref.so, HMY127.ref.so, HMY362.ref.so, test.string='G:G+E')
filtered.G.test.results <- filter.lrt.test.results(S288c.ref.so, HMY12.ref.so, HMY127.ref.so, HMY362.ref.so, test.string='E:G+E')

annotated.GxE.test.results <- left_join(filtered.GxE.test.results, rna.annotation)
annotated.G.test.results <- left_join(filtered.G.test.results, rna.annotation)
annotated.E.test.results <- left_join(filtered.E.test.results, rna.annotation)

## write results to file.
write.csv(annotated.GxE.test.results,file=file.path(proj.dir,"results","GxE_LRT.csv"))
write.csv(annotated.G.test.results,file=file.path(proj.dir,"results","G_LRT.csv"))
write.csv(annotated.E.test.results,file=file.path(proj.dir,"results","E_LRT.csv"))

################ Repeat, without S288c ref results.
filter.lrt.test.results2 <- function(HMY12.so, HMY127.so, HMY362.so, test.string='G+E:G*E', threshold=THRESHOLD) {

  HMY12.test.results <- sleuth_results(HMY12.so, test.string, 'lrt', show_all = FALSE) %>% dplyr::filter(qval <= threshold)
  HMY127.test.results <- sleuth_results(HMY127.so, test.string, 'lrt', show_all = FALSE) %>% dplyr::filter(qval <= threshold)
  HMY362.test.results <- sleuth_results(HMY362.so, test.string, 'lrt', show_all = FALSE) %>% dplyr::filter(qval <= threshold)
  
  filtered.test.results <- HMY12.test.results %>%
  filter(target_id %in% HMY127.test.results$target_id) %>%
  filter(target_id %in% HMY362.test.results$target_id)

  passed <- length(filtered.test.results$target_id)/length(HMY12.test.results$target_id)
  print(passed)
  print("percent passed filtering")
  return(filtered.test.results)
}

filtered.GxE.test.results2 <- filter.lrt.test.results2(HMY12.ref.so, HMY127.ref.so, HMY362.ref.so, test.string='G+E:G*E')
filtered.E.test.results2 <- filter.lrt.test.results2(HMY12.ref.so, HMY127.ref.so, HMY362.ref.so, test.string='G:G+E')
filtered.G.test.results2 <- filter.lrt.test.results2(HMY12.ref.so, HMY127.ref.so, HMY362.ref.so, test.string='E:G+E')

annotated.GxE.test.results2 <- left_join(filtered.GxE.test.results2, rna.annotation)
annotated.G.test.results2 <- left_join(filtered.G.test.results2, rna.annotation)
annotated.E.test.results2 <- left_join(filtered.E.test.results2, rna.annotation)

## What genes are in the G, E, and GxE sets?
E.G.GxE.intersection.results <- annotated.E.test.results2 %>%
    filter(gene %in% annotated.G.test.results2$gene) %>%
    filter(gene %in% annotated.GxE.test.results2$gene)

## ask if there's a significant overlap between E, G, and GxE gene sets.
## yes! p < 0.00001
calc.intersection.randomization.pval(all.genes,
                                     length(E.G.GxE.intersection.results$gene),
                                     length(annotated.E.test.results2$gene),
                                     length(annotated.G.test.results2$gene),
                                     length(annotated.GxE.test.results2$gene))

write.csv(E.G.GxE.intersection.results,file="../results/E-G-GxE-intersection-results.csv")

## genes
G.genes <- annotated.G.test.results2$gene
E.genes <- annotated.E.test.results2$gene
GxE.genes <- annotated.GxE.test.results2$gene
E.G.GxE.intersection <- Reduce(intersect,list(G.genes,E.genes,GxE.genes))

## make a Venn Diagram for Figure 3C
## TODO: change color scheme!
Fig3C <- venn.diagram(x=list(E.genes, G.genes, GxE.genes),
                      category.names = c(paste0("E",
                                                '(',length(E.genes),')'),
                                         paste0("G "
                                               ,'(',length(G.genes),')'),
                                         paste0("GxE ",
                                                '(',length(GxE.genes),')')),
                     filename = NULL,
                     output = TRUE ,
                     imagetype="svg" ,
                     height = 480 , 
                     width = 480 , 
                     resolution = 300,
                     compression = "lzw",
                     lwd = 1,
                     col=c("blue", 'purple', 'red'),
                     fill = c(alpha("blue",0.3), alpha('purple',0.3), alpha('red',0.3)),
                     cex=2,
                     fontfamily = "sans",
                     cat.cex=1.5,
                     cat.default.pos = "outer",
                     cat.pos = c(-27, 27, 135),
                     cat.dist = c(0.055, 0.055, 0.085),
                     cat.fontfamily = "sans",
                     cat.col = c("blue", 'purple', 'red'),
                     rotation = 1
                     )

ggsave(Fig3C,file=file.path(proj.dir,"results","figures","vennGxE.pdf"))

## These results look correct based on functional genomics.
## Tons of GxE. Verify that this is not an artifact
## based on how the design matrix is coded.
design_matrix(S288c.ref.so, which_model="G*E")
design_matrix(S288c.ref.so, which_model="G+E")

prep.transcript.matrix <- function(so, rna.annotation,logscale=TRUE) {
    ## Here I'm hacking some code from the internals of the plot_transcript_heatmap
    ## function, so that I can use a dataframe rather than a sleuth object to make
    ## the figure.
    ## sleuth::plot_transcript_heatmap
    ## see documentation at:
    ## https://www.rdocumentation.org/packages/sleuth/versions/0.30.0/topics/plot_transcript_heatmap
    ## see general documentation at:
    ## https://www.rdocumentation.org/packages/sleuth/versions/0.30.0/

    all.genes <- sort(rna.annotation$gene)
    all.rnas <- sort(rna.annotation$target_id)

    results.matrix <- sleuth_to_matrix(so, "obs_norm", "est_counts")

    ## change the row names of this matrix for better plotting.
    ## first turn into a dataframe, 
    results.matrix.2 <- data.frame(results.matrix) %>%
        rownames_to_column(var = 'target_id') %>%
        left_join(rna.annotation) %>%
        ## and convert back to plot the heatmap.
        column_to_rownames(var = 'gene') %>%
        select(-target_id, -chromosome, -chr_start, -chr_end, -strand, -gene_biotype,
               -transcript_biotype, -gene_symbol, -description)

    if (logscale) {
        ## transform counts using: log2(1+est_counts)
        trans_mat <- as.matrix(log2(results.matrix.2 + 1))
        ## OR just take the log2? zeros can be problematic...
        ##trans_mat <- as.matrix(log2(results.matrix.2))
    } else {
        trans_mat <- as.matrix(results.matrix.2)
    }

}

### Make a nice heatmap of relevant transcript expression.
plot.transcript.matrix <- function(trans_mat, graph.title, log.colors=TRUE, label.genes=TRUE) {

    paletteLength <- 100
    if (log.colors) {
        ## default from pheatmap:
        ## use this red white blue scheme for log-fold change.
        colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(paletteLength)
    } else { ## default sleuth color scheme, used here for
        ## log(1+x) transformed est_counts.
        colors <- colorRampPalette(brewer.pal(n = 7, name="BuPu"))(paletteLength)
    }

    ## don't cluster rows if there's less than 3 of them!
    clust_rows_bool <- ifelse(nrow(trans_mat) > 2, TRUE, FALSE)

    if (log.colors) {
        ## if log.colors == TRUE, make sure the diverging color scale is
        ## symmetric around 0.
        ## see this solution from stackoverflow:
        ## https://stackoverflow.com/questions/31677923/set-0-point-for-pheatmap-in-r
        myBreaks <- c(seq(min(trans_mat), 0, length.out=ceiling(paletteLength/2) + 1), 
                      seq(max(trans_mat)/paletteLength, max(trans_mat), length.out=floor(paletteLength/2)))
        
        p <- pheatmap::pheatmap(trans_mat,
                                color=colors,
                                breaks=myBreaks,
                                cluster_cols = TRUE,
                                cluster_rows=clust_rows_bool,
                                show_rownames=label.genes,
                                main=graph.title)
    } else {
        p <- pheatmap::pheatmap(trans_mat,
                                color=colors,
                                cluster_cols = TRUE,
                                cluster_rows=clust_rows_bool,
                                show_rownames=label.genes,
                                main=graph.title)
    }
    x_axis_angle <- 50
    p$gtable$grobs[[3]]$rot <- 360 - x_axis_angle

    ## return a Gtable object for arranging with cowplot.
    return(p$gtable)
}

filter.transmat <- function(trans_mat,discussed.genes) {
    ## This one-liner removes rows of the trans_mat which are not in discussed.genes.
    return(subset(trans_mat, rownames(trans_mat) %in% discussed.genes))
}

## IMPORTANT TODO: BUG FIX: make sure the gene sets are CORRECT for publication.
## This is based off genes in S288c.
trans_mat <- prep.transcript.matrix(S288c.ref.so,rna.annotation)
## Supplementary File S3: transcript abundance for each sample.
transcript_counts_df <- data.frame(trans_mat) %>%
    rownames_to_column(var = 'gene') %>% left_join(rna.annotation)
write.csv(transcript_counts_df,file="../results/FileS3-S288c-ref-abundance.csv")

## Let's rank genes by transcript abundance in each condition and each strain.
## take geometric mean over all strains and environments.
gm_mean <- function(x, na.rm=TRUE) {
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

## calculate geometric.mean
gm_abundance <- data.frame("gm_est_counts"=apply(trans_mat, 1, gm_mean)) %>%
    rownames_to_column(var = 'gene') %>% arrange(desc(gm_est_counts)) %>%
    left_join(rna.annotation)
    
head(gm_abundance,15)

## calculate arithmetic mean
am_abundance <- data.frame("am_est_counts"=apply(trans_mat, 1, mean)) %>%
    rownames_to_column(var = 'gene') %>% arrange(desc(am_est_counts)) %>%
    left_join(rna.annotation)

## write to file.
write.csv(gm_abundance,file="../results/gm_abundance.csv")
write.csv(am_abundance,file="../results/am_abundance.csv")

## FLO11 (MUC1) is the 9th most abundant transcript over all conditions!

## let's look at each condition separately (hack for now).

## FLO11 is #9 in LD.
LD_trans_mat <- trans_mat[,c(1,3,5,7,9,11,13,15,17)]
LD_abundance <- data.frame("am_est_counts"=apply(LD_trans_mat, 1, mean)) %>%
    rownames_to_column(var = 'gene') %>% arrange(desc(am_est_counts)) %>%
    left_join(rna.annotation)
head(LD_abundance,15)

## FLO11 is #13 in YPD.
YPD_trans_mat <- trans_mat[,c(2,4,6,8,10,12,14,16,18)]
YPD_abundance <- data.frame("am_est_counts"=apply(YPD_trans_mat, 1, mean)) %>%
    rownames_to_column(var = 'gene') %>% arrange(desc(am_est_counts)) %>%
    left_join(rna.annotation)
head(YPD_abundance,15)

## what is rank correlation between LD and YPD?
## in order to do a rank correlation test:
LD.rank.vec <- match(LD_abundance$gene, rna.annotation$gene)
YPD.rank.vec <- match(YPD_abundance$gene, rna.annotation$gene)
## no rank correlation in abundance between LD and YPD.
## tau = -0.00461, p-value = 0.574.
cor.test(LD.rank.vec, YPD.rank.vec, method="kendall")

################################################################################

## Pair samples by their replicate, and calculate the log-fold-change on est_counts.
## Then, make a heatmap.
prep.log.fold.change.matrix <- function(so, rna.annotation) {
    ## Here I'm hacking some code from the internals of the plot_transcript_heatmap
    ## function, so that I can use a dataframe rather than a sleuth object to make
    ## the figure.
    ## sleuth::plot_transcript_heatmap
    ## see documentation at:
    ## https://www.rdocumentation.org/packages/sleuth/versions/0.30.0/topics/plot_transcript_heatmap
    ## see general documentation at:
    ## https://www.rdocumentation.org/packages/sleuth/versions/0.30.0/

    results.matrix <- sleuth_to_matrix(so, "obs_norm", "est_counts")

    ## change the row names of this matrix for better plotting.
    ## first turn into a dataframe, 
    results.matrix.2 <- data.frame(results.matrix) %>%
        rownames_to_column(var = 'target_id') %>%
        ## filter out all transcripts with an est_count < 1 in any sample.
        filter(HMY127aLD > 1 &
               HMY127aYPD > 1 &
               HMY127bLD > 1 &
               HMY127bYPD > 1 &
               HMY127cLD > 1 &
               HMY127cYPD > 1 &
               HMY12aLD > 1 &
               HMY12aYPD > 1 &
               HMY12bLD > 1 &
               HMY12bYPD > 1 &
               HMY12cLD > 1 &
               HMY12cYPD > 1 &
               HMY362aLD > 1 &
               HMY362aYPD > 1 &
               HMY362bLD > 1 &
               HMY362bYPD > 1 &
               HMY362cLD > 1 &
               HMY362cYPD > 1) %>%
        ## quick and dirty solution to get log fold-change.
        ## for HMY12.
        mutate(HMY12a=log2(HMY12aLD/HMY12aYPD)) %>%
        select(-HMY12aLD,-HMY12aYPD) %>%
        mutate(HMY12b=log2(HMY12bLD/HMY12bYPD)) %>%
        select(-HMY12bLD,-HMY12bYPD) %>%
        mutate(HMY12c=log2(HMY12cLD/HMY12cYPD)) %>%
        select(-HMY12cLD,-HMY12cYPD) %>%
        ## for HMY127.
        mutate(HMY127a=log2(HMY127aLD/HMY127aYPD)) %>%
        select(-HMY127aLD,-HMY127aYPD) %>%
        mutate(HMY127b=log2(HMY127bLD/HMY127bYPD)) %>%
        select(-HMY127bLD,-HMY127bYPD) %>%
        mutate(HMY127c=log2(HMY127cLD/HMY127cYPD)) %>%
        select(-HMY127cLD,-HMY127cYPD) %>%
        ## for HMY362.
        mutate(HMY362a=log2(HMY362aLD/HMY362aYPD)) %>%
        select(-HMY362aLD,-HMY362aYPD) %>%
        mutate(HMY362b=log2(HMY362bLD/HMY362bYPD)) %>%
        select(-HMY362bLD,-HMY362bYPD) %>%
        mutate(HMY362c=log2(HMY362cLD/HMY362cYPD)) %>%
        select(-HMY362cLD,-HMY362cYPD) %>%
        ## join with the annotation.
        left_join(rna.annotation) %>%
        ## and convert back to plot the heatmap.
        column_to_rownames(var = 'gene') %>%
        select(-target_id, -chromosome, -chr_start, -chr_end, -strand, -gene_biotype,
               -transcript_biotype, -gene_symbol, -description)
    
    trans_mat <- as.matrix(results.matrix.2)
}

logfold.mat <- prep.log.fold.change.matrix(S288c.ref.so,rna.annotation)

## what is the expected variance in expression over all genes?
my.logfold.variance <- apply(logfold.mat, 1, var, na.rm=TRUE)
logfold.var.df <- data.frame(logfold.variance=my.logfold.variance) %>%
    rownames_to_column(var = 'gene') %>%
    left_join(rna.annotation) %>%
    arrange(desc(logfold.variance)) %>%
    mutate(index=seq_len(nrow(.)))
## write to file.
write.csv(logfold.var.df,file=file.path(proj.dir,"results","logfold_variance.csv"))
## let's plot the variance ranking.
logfold.variance.plot <- ggplot(logfold.var.df,aes(x=index,y=log(logfold.variance))) +
    theme_classic() + geom_point() + ggtitle("log(logfold-variance)") + ylim(-10,10)
ggsave("../results/figures/log-fold-variance.pdf",logfold.variance.plot)

## HMY12 is sporulation deficient. Remove HMY12 and see what genes are highly variable
## in this response to nutrient stress.
## NOTE: this filtering step to remove HMY12 columns is not robust to a reordering
## of columns. but works for now.
noHMY12.logfold.mat <- logfold.mat[,4:ncol(logfold.mat)]
noHMY12.logfold.variance <- apply(noHMY12.logfold.mat, 1, var, na.rm=TRUE)

noHMY12.logfold.var.df <- data.frame(logfold.variance=noHMY12.logfold.variance) %>%
    rownames_to_column(var = 'gene') %>%
    left_join(rna.annotation) %>%
    arrange(desc(logfold.variance)) %>%
    mutate(index=seq_len(nrow(.)))
## write to file.
write.csv(noHMY12.logfold.var.df,file=file.path(proj.dir,"results","noHMY12_logfold_variance.csv"))
## let's plot the variance ranking.
noHMY12.logfold.variance.plot <- ggplot(noHMY12.logfold.var.df,aes(x=index,y=log(logfold.variance))) +
    theme_classic() + geom_point() + ggtitle("log(logfold-variance)") + ylim(-10,10)
ggsave("../results/figures/noHMY12-log-fold-variance.pdf",noHMY12.logfold.variance.plot)
## Interesting result: STRING-DB still shows sporulation-associated genes,
## as well as cell wall, 
## meiotic spindle shows up: selfish meiotic drivers showing up?
## Genes in this paper are enriched:
## https://www.ncbi.nlm.nih.gov/pubmed/27836908


## what is the identity of genes that are highly variable here?
## we have only 57 genes in this category.
highly.logfold.variable.df <- filter(logfold.var.df,logfold.variance>5)
## write to file.
write.csv(highly.logfold.variable.df,file=file.path(proj.dir,"results","highly_variable_logfold.csv"))
## These are all meiosis associated genes! STRING-db analysis shows this very clearly.
## However, meiosis associated genes are also the most variably expressed genes
## outside of the intersection.


logfoldVariance <- apply(logfold.mat, 1, var, na.rm=TRUE)
logfold.var.df <- data.frame(logfold.variance=logfoldVariance) %>%
    rownames_to_column(var = 'gene') %>%
    left_join(rna.annotation) %>%
    arrange(desc(logfold.variance)) %>%
    mutate(index=seq_len(nrow(.)))

################################################################################
###### querying different genes.
S288c.kallisto.df <- kallisto_table(S288c.ref.so,
                                    use_filtered = TRUE,
                                    normalized = TRUE,
                                    include_covariates = TRUE)

HMY12.kallisto.df <- kallisto_table(HMY12.ref.so,
                                    use_filtered = TRUE,
                                    normalized = TRUE,
                                    include_covariates = TRUE)

## look at genes mentioned by Narode et al. (2014).
## in any case, these papers report that FLO11 is necessary
## but not sufficient for biofilm/mat formation.
## SKN7: YHR206W: down in HMY127, up in HMY12, down in 362.
filter(HMY12.kallisto.df,target_id=='YHR206W_mRNA')
## SLG1: YOR008C: upregulated throughout.
filter(HMY12.kallisto.df,target_id=='YOR008C_mRNA')
## VPS27: YNR006W: upregulated throughout.
filter(HMY12.kallisto.df,target_id=='YNR006W_mRNA')


## IS FLO10 being filtered out of HMY12 due to mapping biases? Let's find out.
filter(HMY12.kallisto.df,target_id=='YKR102W_mRNA')
## look at FLO11 too.
filter(HMY12.kallisto.df,target_id=='YIR019C_mRNA')
## look at FLO8... not present?
filter(HMY12.kallisto.df,target_id=='YER108C_mRNA')
## look at TPK2 (pKA)
filter(HMY12.kallisto.df,target_id=='YPL203W_mRNA')
## look at MIT1.
filter(HMY12.kallisto.df,target_id=='YEL007W_mRNA')
## Look at SFL1.
filter(HMY12.kallisto.df,target_id=='YOR140W_mRNA')

## Let's look at genes in Octavio et al. PLOS Genetics 2009, called
## Epigenetic and Conventional regulation is distributed among activators of FLO11
## allowing tuning of population level heterogeneity in its expression.

## Msn1
filter(HMY12.kallisto.df,target_id=='YOL116W_mRNA')
## Rpd3L
filter(HMY12.kallisto.df,target_id=='YNL330C_mRNA')
## Phd1
filter(HMY12.kallisto.df,target_id=='YKL043W_mRNA')
## Mss11
filter(HMY12.kallisto.df,target_id=='YMR164C_mRNA')
## Ste12
filter(HMY12.kallisto.df,target_id=='YHR084W_mRNA')
## Tec1
filter(HMY12.kallisto.df,target_id=='YBR083W_mRNA')
## IRA1
filter(HMY12.kallisto.df,target_id=='YBR140C_mRNA')
## IRA2
filter(HMY12.kallisto.df,target_id=='YOL081W_mRNA')
## GPA2
filter(HMY12.kallisto.df,target_id=='YER020W_mRNA')
## BCY1
filter(HMY12.kallisto.df,target_id=='YIL033C_mRNA')
## TPK1
filter(HMY12.kallisto.df,target_id=='YJL164C_mRNA')
## TPK3
filter(HMY12.kallisto.df,target_id=='YKL166C_mRNA')
## KRH1 (GPB2)
filter(HMY12.kallisto.df,target_id=='YAL056W_mRNA')
## Gpa2p
filter(HMY12.kallisto.df,target_id=='YER020W_mRNA')
## KRH2
filter(HMY12.kallisto.df,target_id=='YOR371C_mRNA')
## CYR1
filter(HMY12.kallisto.df,target_id=='YJL005W_mRNA')
## Take a look at this subunit of TORC complex Tor1p
filter(HMY12.kallisto.df,target_id=='YPL180W_mRNA')

## let's look at genes that may have been deleted from the end of chromosome XI
## in HMY12.

## FLO10, NFT1, YKR104W, VBA5, GEX2
## FLO10: YKR102W
filter(HMY12.kallisto.df,target_id=='YKR102W_mRNA')
## NFT1: YKR103W
filter(HMY12.kallisto.df,target_id=='YKR103W_mRNA')
## YKR104W
filter(HMY12.kallisto.df,target_id=='YKR104W_mRNA')
## VBA5: YKR105C
filter(HMY12.kallisto.df,target_id=='YKR105C_mRNA')
## GEX2: YKR106W
filter(HMY12.kallisto.df,target_id=='YKR106W_mRNA')

################################################################################
## 1) make heatmaps of expression of genes in pathways and datasets of interest.
## 2) figure out how to do comparisons with, and hypothesis testing of genes
## identified in these other studies.

## split pathway.data on each pathway, then make heatmaps for each.

make_transmat_heatmap <- function(tmat, lg.color, pathway.df) {
    genes.of.interest <- pathway.df$gene
    filtered.mat <- filter.transmat(tmat, genes.of.interest)
    pathway.title <- unique(pathway.df$pathway)
    plot.transcript.matrix(filtered.mat, pathway.title, log.colors=lg.color)
}

## make a heatmap for a pathway using the full expression matrix.
all_sample_pathway_heatmap <- partial(.f=make_transmat_heatmap, trans_mat, FALSE)
## make a heatmap for a pathway using the logfold expression matrix.
logfold_pathway_heatmap <- partial(.f=make_transmat_heatmap, logfold.mat, TRUE)

## pathway data that Helen collated.
pathway.data <- read.csv("../data/Pathways.csv",as.is=TRUE,header=TRUE)

## make heatmaps for every pathway, plotting expression in all samples.
pathway.expr.heatmaps <- pathway.data %>% split(.$pathway) %>%
    map(.f = all_sample_pathway_heatmap)

## make heatmaps for every pathway, plotting logfold-expression.
pathway.logchange.heatmaps <- pathway.data %>% split(.$pathway) %>%
    map(.f = logfold_pathway_heatmap)

## plot the location of each gene in the pathway in the distribution
## of gene expression logfold variance.

plot.logfold.variance <- function(my.logfold.mat, pathway.df) {
    genes.of.interest <- pathway.df$gene
    pathway.name <- unique(pathway.df$pathway)

    my.logfold.variance <- apply(my.logfold.mat, 1, var, na.rm=TRUE)    
    logfold.var.df <- data.frame(logfold.variance=my.logfold.variance) %>%
        rownames_to_column(var = 'gene') %>%
        left_join(rna.annotation) %>%
        arrange(desc(logfold.variance)) %>%
        mutate(index=seq_len(nrow(.))) %>%
        mutate(in.pathway = gene %in% genes.of.interest)
    ## let's plot the variance ranking.
    logfold.variance.plot <- ggplot(logfold.var.df,
                                    aes(x=index,
                                        y=log2(logfold.variance))) +
        scale_color_manual(values = c("grey90","black")) +
        theme_classic() +
        geom_line() + ## to handle overplotting
        ggtitle(paste0(pathway.name," log(logfold-variance)")) +
        ylim(-10,10) +
        xlim(0,6000) +
        geom_point(data=filter(logfold.var.df,in.pathway==TRUE)) +
        geom_text_repel(data=filter(logfold.var.df,in.pathway==TRUE),aes(label=gene)) +
        geom_vline(xintercept=500,linetype='dashed')
    return(logfold.variance.plot)
}

## plot log-fold variance for each pathway.
## Flocculins are what we care about here.
pathway.data %>% split(.$pathway) %>%
    map(.f=partial(plot.logfold.variance,logfold.mat))

Fig10C <- plot.logfold.variance(logfold.mat,filter(pathway.data,pathway=='Flocculins'))


###########################################################################
## let's analyze and compare with the datasets from four papers.

## This includes the hypotheses generated from Aimee Dudley's paper,
## and the hypotheses generated from the Chow filamentous growth paper.

###########################################################################
## Stovicek et al. 2014.
Stovicek.data <- read.csv('../data/Stovicek2014Table1.csv', as.is=TRUE,header=TRUE) %>%
    mutate(direction=ifelse(expected.change.in.biofilm=='upregulated',1,-1))

make_Stovicek_heatmap <- function(tmat, lg.color, express.df) {
    genes.of.interest <- express.df$gene
    filtered.mat <- filter.transmat(tmat, genes.of.interest)
    fname <- unique(express.df$expected.change.in.biofilm)
    title <- paste(fname, 'predictions')
    plot.transcript.matrix(filtered.mat,title,
                           log.colors = lg.color)
}

## make separate heatmaps for genes that are predicted to be upregulated,
## or downregulated.

## first plot expression in all samples.
stovicek.plot.list <- Stovicek.data %>%
    split(.$expected.change.in.biofilm) %>%
    map(.f=partial(.f=make_Stovicek_heatmap, trans_mat, FALSE))

## make heatmaps for every pathway, plotting logfold-expression.
log.stovicek.plot.list <- Stovicek.data %>% split(.$expected.change.in.biofilm) %>%
    map(.f=partial(.f=make_Stovicek_heatmap, logfold.mat,TRUE))

## compare direction of mean logfold change in the Stovicek genes
## to the direction predicted by the Stovicek data.
upregulated.Stovicek.data <- filter(Stovicek.data,direction==1)
downregulated.Stovicek.data <- filter(Stovicek.data,direction==-1)
upregulated.Stovicek.genes.mat <- filter.transmat(logfold.mat,upregulated.Stovicek.data$gene)
downregulated.Stovicek.genes.mat <- filter.transmat(logfold.mat,downregulated.Stovicek.data$gene)

## By inspection, we can reject the hypothesis that the directionality of the
## gene expression change in the Stovicek gene set (reported in bold in their Table 1)
## is conserved in HMY12, 127, 362.

upregulated.directions <- data.frame(upregulated.Stovicek.genes.mat) %>%
    rownames_to_column(var = 'target_id') %>%
    rowwise() %>%
    mutate(HMY12=ifelse(mean(HMY12a,HMY12b,HMY12c) > 0, 1, -1)) %>%
    mutate(HMY127=ifelse(mean(HMY127a,HMY127b,HMY127c) > 0, 1, -1)) %>%
    mutate(HMY362=ifelse(mean(HMY362a,HMY362b,HMY362c) > 0, 1, -1)) %>%
    select(target_id, HMY12, HMY127, HMY362)


downregulated.directions <- data.frame(downregulated.Stovicek.genes.mat) %>%
    rownames_to_column(var = 'target_id') %>%
    rowwise() %>%
    mutate(HMY12=ifelse(mean(HMY12a,HMY12b,HMY12c) > 0, 1, -1)) %>%
    mutate(HMY127=ifelse(mean(HMY127a,HMY127b,HMY127c) > 0, 1, -1)) %>%
    mutate(HMY362=ifelse(mean(HMY362a,HMY362b,HMY362c) > 0, 1, -1)) %>%
    select(target_id, HMY12, HMY127, HMY362)

## Make a Figure 4 and Figure S2, showing results of testing predictions.
Fig4 <- plot_grid(log.stovicek.plot.list$upregulated,
                  log.stovicek.plot.list$downregulated)
Fig4
ggsave("../results/figures/Fig4.pdf",Fig4)


FigS2 <- plot_grid(stovicek.plot.list$upregulated,
                  stovicek.plot.list$downregulated) 
FigS2
ggsave("../results/figures/FigS2.pdf",FigS2,width=9)


################################################
## Voordeckers et al. 2012

## Figure 5A and 5B.

## Voordeckers et al. 2012 found 211 deletions affecting colony morphology.
## 52 cause smooth colony phenotype (code == 1).
## The remaning 159 reduce wrinkliness (code == 2).

## Hence, there's a weak argument that if deletions reduce the phenotype,
## then we expect those genes to be upregulated in YPD.

deletion.screen.data <- read.csv("../data/Voordeckers2012_KnockoutResults.csv",as.is=TRUE)

code1.screen.data <- deletion.screen.data %>% filter(code==1)
code2.screen.data <- deletion.screen.data %>% filter(code==2)

make_heatmap <- function(tmat, relevant.genes, plot.title, lg.colors,label.genes=TRUE) {
    filtered.mat <- filter.transmat(tmat, relevant.genes)
    plot.transcript.matrix(filtered.mat,plot.title, lg.colors,label.genes)
}

code1.log.gtable <- make_heatmap(logfold.mat, code1.screen.data$gene, "smooth when knocked out", TRUE)
code1.gtable <- make_heatmap(trans_mat, code1.screen.data$gene, "smooth when knocked out", FALSE)

code2.log.gtable <- make_heatmap(logfold.mat, code2.screen.data$gene, "less wrinkly when knocked out", TRUE,FALSE)
code2.gtable <- make_heatmap(trans_mat, code2.screen.data$gene, "less wrinkly when knocked out", FALSE,FALSE)

## Figure 5C.
## key genes mentioned by Voordeckers et al.:
##Flo11, Dfg16, Tos1. genes for these are here:
## Flo11 = YIR019C; Dfg16 = YOR030W; Tos1 = YBR162C.
key.genes <- c("YIR019C", "YOR030W", "YBR162C")

## Looks like responses of these "key" genes identified by Voordeckers is not conserved!
key.log.gtable <- make_heatmap(logfold.mat, key.genes, "FLO11, DFG16, TOS1", TRUE)
key.gtable <- make_heatmap(trans_mat, key.genes, "FLO11, DFG16, TOS1", FALSE)


Fig5 <- plot_grid(code1.log.gtable,
                  code2.log.gtable,
                  key.log.gtable,nrow=1)
Fig5
ggsave("../results/figures/Fig5.pdf",Fig5,height=8)

FigS3 <- plot_grid(code1.gtable,
                  code2.gtable,
                  key.gtable,nrow=1)
FigS3
ggsave("../results/figures/FigS3.pdf",FigS3, height=8)

## dig a bit deeper into Dfg16, Tos1, and FLO11.

## Are Dfg16 and Tos1 significant in these data?
"YOR030W" %in% E.G.GxE.intersection ## YES for Dfg16!
"YOR030W" %in% G.genes ## YES for Dfg16
"YOR030W" %in% E.genes ## YES for Dfg16.
"YOR030W" %in% GxE.genes ## YES for Dfg16.
filter(annotated.G.test.results,gene== "YOR030W")
filter(annotated.E.test.results,gene== "YOR030W")
filter(annotated.GxE.test.results,gene== "YOR030W")

## Dfg16 is upregulated overall-- but weakly downregulated in
## HMY362 in LD conditions!
filter(HMY12.kallisto.df,target_id=='YOR030W_mRNA')

"YBR162C" %in% E.G.GxE.intersection ## Yes for Tos1.
"YBR162C" %in% G.genes ## Yes for Tos1
"YBR162C" %in% E.genes ## Yes for Tos1
"YBR162C" %in% GxE.genes ## Yes for Tos1
filter(annotated.G.test.results,gene== "YBR162C")
filter(annotated.E.test.results,gene== "YBR162C")
filter(annotated.GxE.test.results,gene== "YBR162C")

## Tos1 IS DOWNREGULATED IN ALL THREE STRAINS!
filter(HMY12.kallisto.df,target_id=='YBR162C_mRNA')

## Check Flo11.
"YIR019C" %in% E.G.GxE.intersection ## Yes for Flo11
"YIR019C" %in% G.genes ## YES for Flo11
"YIR019C" %in% E.genes ## NO for Flo11.
"YIR019C" %in% GxE.genes ## YES for Flo11.
filter(annotated.G.test.results,gene== "YIR019C")
filter(annotated.E.test.results,gene== "YIR019C")
filter(annotated.GxE.test.results,gene== "YIR019C")

## Are Dfg16 and Tos1 overexpressed in LD conditions?
## NO! Tos1 is downregulated in all three strains!
## And different signs of fold-change in Flo11!

##################################################################################
## Look at genes reported in Chow, Cullen 2019 Genetics paper.
## ste12, ras2, rtg3, spt8.
## these are: YHR084W, YNL098C, YBL103C, YLR055C, respectively.
chow.genes <- c("YHR084W", "YNL098C", "YBL103C", "YLR055C")
log.Chow.gtable <- make_heatmap(logfold.mat,chow.genes,"ste12, ras2, rtg3, and spt8",TRUE)
Chow.gtable <- make_heatmap(trans_mat, chow.genes, "ste12, ras2, rtg3, and spt8", TRUE)

Fig6 <- log.Chow.gtable
ggsave("../results/figures/Fig6.pdf",Fig6)

FigS4 <- Chow.gtable
ggsave("../results/figures/FigS4.pdf",FigS4)

## again-- poor clustering of samples by envirionment or genotype!

## In introduction, Chow mention this gene: YOR008C
## as an independent pathway from FLO11. Is there evidence of this
## in the strain where FLO11 is downregulated?

##################################################################################
## Cromie et al. (2017) G3 paper.

## The authors identified
## DIG1, SFL1, HEK2, TOS8, SAN1, and ROF1/YHR177W
## as regulators of biofilm formation in the F45 strain background.

## DIG1, SFL1, HEK2, TOS8, SAN1, and ROF1/YHR177W:
##YPL049C, YOR140W, YBL032W, YGL096W, YDR143C, YHR177W, respectively.

cromie.regulators <- c("YPL049C", "YOR140W", "YBL032W", "YGL096W", "YDR143C", "YHR177W")

log.cromie.regulator.gtable <- make_heatmap(logfold.mat,cromie.regulators,"DIG1, SFL1, HEK2, TOS8, SAN1, and ROF1",TRUE)

cromie.regulator.gtable <- make_heatmap(trans_mat,cromie.regulators,"DIG1, SFL1, HEK2, TOS8, SAN1, and ROF1",TRUE)

Fig7 <- plot_grid(log.cromie.regulator.gtable,cromie.regulator.gtable,labels=c('A','B'))
ggsave("../results/figures/Fig7.pdf",Fig7)


## They then inferred a latent "Common Factor" that may represent a common
## transcriptional state induced by these regulators (except for TOS8.)

Cromie.S2.up.genes <- read.csv("../data/Cromie-S2-induced.csv")$gene
Cromie.S2.down.genes <- read.csv("../data/Cromie-S2-repressed.csv")$gene
Cromie.S4.up.genes <- read.csv("../data/Cromie-S4-induced.csv")$gene
Cromie.S4.down.genes <- read.csv("../data/Cromie-S4-repressed.csv")$gene

## look at genes upregulated in the Common Factor.
Cromie.S2.up.log.gtable <- make_heatmap(logfold.mat,Cromie.S2.up.genes,'induced by latent factor',TRUE, FALSE)
Cromie.S2.up.gtable <- make_heatmap(trans_mat,Cromie.S2.up.genes,'induced by latent factor',FALSE, FALSE)

## look at genes downregulated in the Common Factor.
Cromie.S2.down.log.gtable <- make_heatmap(logfold.mat,Cromie.S2.down.genes, 'repressed by latent factor', TRUE, FALSE)
Cromie.S2.down.gtable <- make_heatmap(trans_mat,Cromie.S2.down.genes,'repressed by latent factor',FALSE, FALSE)

## look at genes upregulated in both the Common Factor and the TOS8 overexpression strain.
Cromie.S4.up.log.gtable <- make_heatmap(logfold.mat,Cromie.S4.up.genes,'induced by LF and TOS8',TRUE, FALSE)
Cromie.S4.up.gtable <- make_heatmap(trans_mat,Cromie.S4.up.genes,'induced by LF and TOS8',FALSE, FALSE)

## look at genes downregulated in both the Common Factor and the TOS8 overexpression strain.
Cromie.S4.down.log.gtable <- make_heatmap(logfold.mat,Cromie.S4.down.genes,'repressed by LF and TOS8',TRUE, FALSE)
Cromie.S4.down.gtable <- make_heatmap(trans_mat,Cromie.S4.down.genes,'repressed by LF and TOS8',FALSE, FALSE)

## Make Fig 8, 9 and  S5 showing results of testing predictions.
Fig8 <- plot_grid(Cromie.S2.down.log.gtable, Cromie.S2.up.log.gtable)
ggsave("../results/figures/Fig8.pdf",Fig8)

FigS5 <- plot_grid(Cromie.S2.down.gtable, Cromie.S2.up.gtable)
ggsave("../results/figures/FigS5.pdf",FigS5)

Fig9 <- plot_grid(Cromie.S4.down.log.gtable, Cromie.S4.up.log.gtable)
ggsave("../results/figures/Fig9.pdf",Fig9)

FigS6 <- plot_grid(Cromie.S4.down.gtable, Cromie.S4.up.gtable)
ggsave("../results/figures/FigS6.pdf",FigS6)


## Make figures and supplementary figures showing changes of expression of
## selected genes and pathways of interest.

## Easiest to make one figure for each pathway, and throw all pictures into
## Supplement.

## Make Flocculins into main figure.
Fig10 <- plot_grid(
    pathway.logchange.heatmaps$`Flocculins`,
    pathway.expr.heatmaps$`Flocculins`,
    Fig10C,
    labels=c('A','B','C'))
ggsave("../results/figures/Fig10.pdf",Fig10,height=8,width=8)


FigS7 <- plot_grid(
    pathway.logchange.heatmaps$`RIM101 pathway`,
    pathway.expr.heatmaps$`RIM101 pathway`,
    pathway.expr.heatmaps$TORC2,
    pathway.logchange.heatmaps$TORC2,
    nrow=2,
    labels=c('A','B','C','D'))
ggsave("../results/figures/FigS7.pdf",FigS7,height=8,width=8)


FigS8 <- plot_grid(
    pathway.logchange.heatmaps$`Ras pathway`,
    pathway.expr.heatmaps$`Ras pathway`,
    pathway.logchange.heatmaps$`Ras protein signal transduction`,
    pathway.expr.heatmaps$`Ras protein signal transduction`,
    labels=c('A','B','C','D'))
ggsave("../results/figures/FigS8.pdf",FigS8,height=8,width=8)

FigS9 <- plot_grid(
    pathway.logchange.heatmaps$`SAGA complex`,
    pathway.expr.heatmaps$`SAGA complex`,
    pathway.logchange.heatmaps$`SWI-SNF`,
    pathway.expr.heatmaps$`SWI-SNF`,
    labels=c('A','B','C','D'))
ggsave("../results/figures/FigS10.pdf",FigS9,height=8,width=8)

FigS10 <- plot_grid(
    pathway.logchange.heatmaps$`mitochondria-nucleus signaling pathway`,
    pathway.expr.heatmaps$`mitochondria-nucleus signaling pathway`,
    pathway.logchange.heatmaps$Others,
    pathway.expr.heatmaps$Others,
    labels=c('A','B','C','D'))
ggsave("../results/figures/FigS10.pdf",FigS10,height=8,width=11)

FigS11 <- plot_grid(
    pathway.logchange.heatmaps$`TOR signaling`,
    pathway.expr.heatmaps$`TOR signaling`,
    pathway.logchange.heatmaps$`Transcription Factors`,
    pathway.expr.heatmaps$`Transcription Factors`,
    pathway.logchange.heatmaps$`Zinc Reg`,
    pathway.expr.heatmaps$`Zinc Reg`,
    labels=c('A','B','C','D','E','F'),ncol=2)
ggsave("../results/figures/FigS11.pdf",FigS11,height=11,width=8)

FigS12 <- plot_grid(
    pathway.logchange.heatmaps$`fMAPK`,
    pathway.expr.heatmaps$`fMAPK`,
    pathway.logchange.heatmaps$`HOG pathway`,
    pathway.expr.heatmaps$`HOG pathway`,
    pathway.logchange.heatmaps$LiquidResponse,
    pathway.expr.heatmaps$LiquidResponse,
    labels=c('A','B','C','D','E','F'),ncol=2)
ggsave("../results/figures/FigS12.pdf",FigS12,height=11,width=11)
