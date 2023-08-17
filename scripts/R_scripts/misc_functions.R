#R functions for microbial analyses 
#By Cassie Ettinger

#In an effort to cut some repritive code from my Rmd
#I have tried to make some code into functions here

#standard error function for plotting
se <- function(x) sqrt(var(x)/length(x))
getN <- function(x) sum(getUniques(x))

#inspect quality of first six F read files
inspect_qual_F <- function(raw_data) {
  #Sort and get sample names
  fnFs <- sort(list.files(raw_data, pattern="R1.noprimers.fastq.gz"))
  #specify full paths to the data
  fnFs <- file.path(raw_data, fnFs)
  #Inspecting quality of data
  pQF <- plotQualityProfile(fnFs[1:6]) #fwd reads for first 6 samples
  return(pQF)
}

#inspect quality of first six R read files
inspect_qual_R <- function(raw_data) {
  #Sort and get sample names
  fnRs <- sort(list.files(raw_data, pattern="R2.noprimers.fastq.gz"))
  #specify full paths to the data
  fnRs <- file.path(raw_data, fnRs)
  #Inspecting quality of data
  pQR <- plotQualityProfile(fnRs[1:6]) #reverse reads for first 6 samples
  return(pQR)
}

#function to run through whole ITS dada2 pipeline
run_dada2_ITS <- function(raw_data, filt_path, maxEE_F, maxEE_R, Taxa_DB) {

#Sort and get sample names
fnFs <- sort(list.files(raw_data, pattern="R1.noprimers.fastq.gz"))
fnRs <- sort(list.files(raw_data, pattern="R2.noprimers.fastq.gz"))
sample.names <- sapply(strsplit(fnFs, "_ITS"), `[`, 1)

#specify full paths to the data
fnFs <- file.path(raw_data, fnFs)
fnRs <- file.path(raw_data, fnRs)

#specify where to save filtered data that we will generate and what to name the files
#we will mostly filter to remove any 'N's which dada2 cannot handle
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filtered.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filtered.fastq.gz"))

#Using a maxEE of 2 and truncating when read quality falls below 10
#not truncating reads beyond based on quality bc we want to maintain most length due 
#to ITS length variation
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, maxEE=c(maxEE_F,maxEE_R), truncQ=10, matchIDs =TRUE, rm.phix=TRUE, compress=TRUE, multithread=TRUE, verbose = TRUE)
head(out)

#get error rates
errF <- learnErrors(filtFs, multithread=TRUE, randomize=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE, randomize=TRUE)

#checking if there is an output file (e.g. no missing files)
exists <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]


#get file names
sample.names <- sapply(strsplit(basename(filtFs), "_F_filtered.fastq.gz"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_R_filtered.fastq.gz"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Dereplication & Sample inference 
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, trimOverhang=TRUE)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)

# Make a sequence table 
seqtab <- makeSequenceTable(mergers)


#Remove chimeras
seqtab2 <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)


#assign taxonomy 

tax <- assignTaxonomy(seqtab2, Taxa_DB, multithread=TRUE, tryRC = TRUE)

#but how many reads survived to the end?
# set a little function

# making a table
summary_tab <- data.frame(row.names=sample.names, dada2_input=out[,1],
                          filtered=out[,2],merged=sapply(mergers, getN),
                          nonchim=rowSums(seqtab2),
                          final_perc_reads_retained=round(rowSums(seqtab2)/out[,1]*100, 1))

#returns a list of the fwd error, rev error, sequence table, chimera-filtered seq table, taxonomy, summary table from processing
return(list(errF, errR, seqtab, seqtab2, tax, summary_tab))

}


#function to run through whole 16S dada2 pipeline
run_dada2_16S <- function(raw_data, filt_path, maxEE_F, maxEE_R, trun_F, trun_R, Taxa_DB) {
  
  #Sort and get sample names
  fnFs <- sort(list.files(raw_data, pattern="R1.noprimers.fastq.gz"))
  fnRs <- sort(list.files(raw_data, pattern="R2.noprimers.fastq.gz"))
  sample.names <- sapply(strsplit(fnFs, "_16S"), `[`, 1)
  
  #specify full paths to the data
  fnFs <- file.path(raw_data, fnFs)
  fnRs <- file.path(raw_data, fnRs)
  
  #specify where to save filtered data that we will generate and what to name the files
  #we will mostly filter to remove any 'N's which dada2 cannot handle
  filtFs <- file.path(filt_path, paste0(sample.names, "_F_filtered.fastq.gz"))
  filtRs <- file.path(filt_path, paste0(sample.names, "_R_filtered.fastq.gz"))
  
  #Using a maxEE of 2 and truncating when read quality falls below 10
  #not truncating reads beyond based on quality bc we want to maintain most length due 
  #to ITS length variation
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, maxEE=c(maxEE_F,maxEE_R), truncLen = c(trun_F, trun_R), truncQ=10, matchIDs =TRUE, rm.phix=TRUE, compress=TRUE, multithread=TRUE, verbose = TRUE)
  head(out)
  
  #get error rates
  errF <- learnErrors(filtFs, multithread=TRUE, randomize=TRUE)
  errR <- learnErrors(filtRs, multithread=TRUE, randomize=TRUE)
  
  
  #checking if there is an output file (e.g. no missing files)
  exists <- file.exists(filtFs) & file.exists(filtRs)
  filtFs <- filtFs[exists]
  filtRs <- filtRs[exists]
  
  
  #get file names
  sample.names <- sapply(strsplit(basename(filtFs), "_F_filtered.fastq.gz"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
  sample.namesR <- sapply(strsplit(basename(filtRs), "_R_filtered.fastq.gz"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
  if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names
  
  #Dereplication & Sample inference 
  mergers <- vector("list", length(sample.names))
  names(mergers) <- sample.names
  
  for(sam in sample.names) {
    cat("Processing:", sam, "\n")
    derepF <- derepFastq(filtFs[[sam]])
    ddF <- dada(derepF, err=errF, multithread=TRUE)
    derepR <- derepFastq(filtRs[[sam]])
    ddR <- dada(derepR, err=errR, multithread=TRUE)
    merger <- mergePairs(ddF, derepF, ddR, derepR, trimOverhang=TRUE)
    mergers[[sam]] <- merger
  }
  rm(derepF); rm(derepR)
  
  # Make a sequence table 
  seqtab <- makeSequenceTable(mergers)
  
  
  #Remove chimeras
  seqtab2 <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
  
  
  #assign taxonomy 
  
  tax <- assignTaxonomy(seqtab2, Taxa_DB, multithread=TRUE, tryRC = TRUE)
  
  #but how many reads survived to the end?
  # set a little function
  
  # making a table
  summary_tab <- data.frame(row.names=sample.names, dada2_input=out[,1],
                            filtered=out[,2],merged=sapply(mergers, getN),
                            nonchim=rowSums(seqtab2),
                            final_perc_reads_retained=round(rowSums(seqtab2)/out[,1]*100, 1))
  
  #returns a list of the fwd error, rev error, sequence table, chimera-filtered seq table, taxonomy, summary table from processing
  return(list(errF, errR, seqtab, seqtab2, tax, summary_tab))
  
}


#function to run decontam and return identified contaminant ASVs
run_decontam_threshold <- function(dataframe, threshold, type) {

  sample_data(dataframe)$is.neg <- sample_data(dataframe)$FilePrefix == type
  contamdf.prev <- isContaminant(dataframe, method="prevalence", neg="is.neg",  threshold=threshold)
  
  contams <- row.names(contamdf.prev)[which(contamdf.prev$contaminant)]
  
  ps.pa <- transform_sample_counts(dataframe, function(abund) 1*(abund>0))
  ps.pa.neg <- prune_samples(sample_data(ps.pa)$FilePrefix == type, ps.pa)
  ps.pa.pos <- prune_samples(sample_data(ps.pa)$FilePrefix != type, ps.pa)
  
  df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                          contaminant=contamdf.prev$contaminant)
  contam_plot = ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
    xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
  
#returns a list of contaminant ASVs, plot of contaminant prevalence  
  return(list(contams, contam_plot))
  
}

#Custom script from phyloseq_extended
#not written by me!
#https://github.com/mahendra-mariadassou/phyloseq-extended/blob/master/R/richness.R
ggrare <- function(physeq, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  ## Args:
  ## - physeq: phyloseq class object, from which abundance data are extracted
  ## - step: Step size for sample size in rarefaction curves
  ## - label: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - color: (Optional). Default ‘NULL’. Character string. The name of the
  ##          variable to map to colors in the plot. This can be a sample
  ##          variable (among the set returned by
  ##          ‘sample_variables(physeq)’ ) or taxonomic rank (among the set
  ##          returned by ‘rank_names(physeq)’).
  ##
  ##          Finally, The color scheme is chosen automatically by
  ##          ‘link{ggplot}’, but it can be modified afterward with an
  ##          additional layer using ‘scale_color_manual’.
  ## - color: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - plot:  Logical, should the graphic be plotted.
  ## - parallel: should rarefaction be parallelized (using parallel framework)
  ## - se:    Default TRUE. Logical. Should standard errors be computed. 
  ## require vegan
  x <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) { x <- t(x) }
  
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    #cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  
  ## Get sample data 
  if (!is.null(sample_data(physeq, FALSE))) {
    sdf <- as(sample_data(physeq), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  ## Add, any custom-supplied plot-mapped variables
  if( length(color) > 1 ){
    data$color <- color
    names(data)[names(data)=="color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  if( length(label) > 1 ){
    labels$label <- label
    names(labels)[names(labels)=="label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot(data = data, aes_string(x = "Size", y = ".S", group = "Sample", color = color))
  p <- p + labs(x = "Sample Size", y = "Species Richness")
  if (!is.null(label)) {
    p <- p + geom_text(data = labels, aes_string(x = "x", y = "y", label = label, color = color),
                       size = 4, hjust = 0)
  }
  p <- p + geom_line()
  if (se) { ## add standard error if available
    p <- p + geom_ribbon(aes_string(ymin = ".S - .se", ymax = ".S + .se", color = NULL, fill = color), alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}

#function that takes in df, category of interest, alpha div metric
#and calculates alpha diversity, then tests whether signficantly different
test_alpha_div <- function(alpha_div_df, category, alphatype) {

    #make a factor
  sample_data(alpha_div_df)[[category]] <- factor(sample_data(alpha_div_df)[[category]])
  
  ## Stats ##
  
  #merged
  alpha_est <- estimate_richness(alpha_div_df, measures = alphatype)
  alpha_est <- cbind(alpha_est, sample_data(alpha_div_df))
  
  # kruskal tests overall to look at sample type
  form <- as.formula(paste(alphatype, category, sep = "~"))
  result <- kruskal_test(form, distribution = approximate(nresample = 9999), 
               data = alpha_est)
  
  
  
  posthoc <- tryCatch({dunnTest(form, data = alpha_est, method = "bh")}, error = {function(x) return(FALSE)})
  
#returns the KW test, posthoc test, and alpha diversity estimates  
  return(list(result, posthoc, alpha_est))
}


#removes prefixes from ITS taxonomy file
fix_ITS_taxonomy <- function(phyloseq_obj) {
  df.ITS.tax <- data.frame(tax_table(phyloseq_obj))
  
  df.ITS.tax %<>% 
    mutate(Phylum = fct_explicit_na(Phylum, na_level = "p__Unclassified"), 
           Class = fct_explicit_na(Class, na_level = "c__Unclassified"), 
           Order = fct_explicit_na(Order, na_level = "o__Unclassified"), 
           Family = fct_explicit_na(Family, na_level = "f__Unclassified"), 
           Genus = fct_explicit_na(Genus, na_level = "g__Unclassified"), 
           Species = fct_explicit_na(Species, na_level = "s__Unclassified"))
  
  tax.list <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax.header <- c(Phylum = "p__", Class = "c__", Order = "o__", 
                  Family = "f__", Genus = "g__", Species = "s__")
  
  for (i in tax.list) {
    names <- sapply(strsplit(as.character(df.ITS.tax[[i]]), as.character(tax.header[[i]])), `[`, 2)
    df.ITS.tax[[i]] <- names 
  }
  
  #row.names(df.ITS.tax) <- row.names(tax_table(ps_OF_nz_rT))
  ITS.tax <- as.matrix(df.ITS.tax)
  
  return(ITS.tax)
}

#function calculates beta diversity, an calculates ordinantion
calculate_beta_div <- function(beta_div_df, beta_metric, rare, method) {
  #takes in phyloseq obj unrareified, beta diversity metric to calculate, rarefaction level, ordination method
  
  dist <- avgdist(as.data.frame(beta_div_df@otu_table), sample=rare, dmethod=beta_metric)
  
  ordin_dist <- ordinate(
    physeq = beta_div_df, 
    method = method, 
    dist = dist)

#returns ordination, beta-div distance matrix
  return(list(ordin_dist, dist))
}


#takes in adonis results and extracts pvalues, Fmodel and R2
get_beta_pvals <- function(results, n) {
  
  if (n > 1) {
    pvals = list()
    Fmods = list()
    Rs = list()
  
    for (i in 1:n) {
    pvals = append(pvals, results$"Pr(>F)"[i])
    Fmods = append(Fmods, results$F[i])
    Rs = append(Rs, results$R2[i])}
  }
  else {
    pvals = results$"Pr(>F)"[n]
    Fmods = results$F[n]
    Rs = results$R2[n] 
    }
    
  return(list(pvals, Fmods, Rs))
}

#takes in betadisper results and extracts pvalues, Fmodel
get_disp_pvals <- function(results) {
  
  pvals = results$tab$"Pr(>F)"[1]
  Fmods = results$tab$F[1]
  
  return(list(pvals, Fmods))
  
}

get_kw_results_rel_abun <- function(DataSet, TaxaLevel, CategoryToTest) {
  #takes in a tibble from phyloseq of raw abundance, needs column with taxonomy level and the category to test
  #returns kw results
  
  taxa_to_test = unique(DataSet[,TaxaLevel])
  
  all.results.df <- data.frame()
  
  for (cat in taxa_to_test) {
    new_df <- DataSet %>%
      mutate(CategoryToTest = as.factor(DataSet[,CategoryToTest])) %>%
      filter(DataSet[,TaxaLevel] == cat) 
    kw = kruskal.test(Abundance ~ CategoryToTest, data = new_df)
    res.df <- data.frame(TaxonomicLevel = TaxaLevel, 
                         Category = CategoryToTest,
                         Taxon = cat,
                         chisq = kw$statistic,
                         pval = kw$p.value)
    
    all.results.df <- bind_rows(res.df, all.results.df)
  }
    
  all.results.df <- all.results.df %>%
    mutate(padj = p.adjust(pval, method = "bonferroni")) 
    
    
  return(all.results.df)
    
}



get_dunn_results_rel_abun <- function(DataSet, TaxaLevel, CategoryToTest) {
  #get post-hoc test results for signifcant KW tests
  
  all.dunn.res <- data.frame()

  taxa_to_test = unique(DataSet[,TaxaLevel])
  
  for (cat in taxa_to_test) {
    new_df <- DataSet %>%
      mutate(CategoryToTest = as.factor(DataSet[,CategoryToTest])) %>%
      filter(DataSet[,TaxaLevel] == cat) 
    
    dT = dunnTest(Abundance ~ CategoryToTest, data = new_df, method = "bonferroni")
    
    for (i in 1:length(dT$res$Comparison)) {
      
      res.df <- data.frame(TaxonomicLevel = TaxaLevel, 
                           Category = CategoryToTest,
                           Taxon = cat,
                           comparison = as.character(dT$res$Comparison)[i], 
                           Zsc = dT$res$Z[i],
                           pval = dT$res$P.unadj[i],
                           padj = dT$res$P.adj[i])
      
      all.dunn.res <- bind_rows(res.df, all.dunn.res)
      
    }
    
  }
  

  return(all.dunn.res)
  }



