
################################################################################
###
###  A script for generating graphs of protein domain structure 
###  and adding various annotations including location of fixed sites,
###  sliding window amino acid diversity plots etc...
###  Requires: 1. An amino acid alignment file in fasta format
###            2. a .csv text file specifying taxonomic groups for species of interest
###            3. Uniprot accession number for protein of interest
###
################################################################################

### v. 1.0.3 bug fixes
### fixed issue with conversion to character matrix in polymorphic sites function
### fixed issue with loop indexing in seqdifs.multi function related to use of character matrix input rather than AAStringSets
### corrected function inputs to as required for AAStrings, AAStringSets, or character matrix
### general code tidying and improved comments

### v 1.0.4 changes
### added script for fixing indels

### Link to Uniprot: https://www.uniprot.org/

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")

####Run this if you need to install necessary packages
library(BiocManager)

aBiocManager::install("drawProteins")
BiocManager::install("Biostrings")
BiocManager::install("ggmsa")
#BiocManager::install.packages("msa")  ## note this is different from ggmsa - but might cause conflicts?
install.packages("Rtools")
install.packages("ggplot2")
install.packages("statebins")
install.packages("RColorBrewer")
install.packages("patchwork")
install.packages("paletteer")

####Open libraries

library(ggplot2)          ### plotting
library(Biostrings)       ### DNA/AA sequence functions
library(drawProteins)     ### Protein drawing functions 
library(statebins)        ### Needed for plotting rectangles with nice rounded corners
library(RColorBrewer)     ### Colour palette functions
library(patchwork)        ### Easy layouts for multiple plot objects
library(ggmsa)            ### Nice sequence alignment plots
library(paletteer)        ### Pretty Colours
#library(msa)             ### Functions for multiple sequence alignments - might cause conflicts with Biostrings?

############ Links to vignettes for key packages ########

# ggplot       - https://ggplot2-book.org/ (free online book)  
# RcolorBrewer - https://jmsallan.netlify.app/blog/the-brewer-palettes/
# DrawProteins - https://bioconductor.org/packages/devel/bioc/vignettes/drawProteins/inst/doc/drawProteins_BiocStyle.html
#              - https://rforbiochemists.blogspot.com/2018/02/my-drawproteins-package-has-been.html
# Biostrings   - http://127.0.0.1:27877/library/Biostrings/doc/BiostringsQuickOverview.pdf
#              - https://bioconductor.org/packages/release/bioc/html/Biostrings.html
#              - https://www.rdocumentation.org/packages/Biostrings/versions/2.40.2
# ggmsa        - http://yulab-smu.top/ggmsa/index.html
# msa          - https://bioconductor.riken.jp/packages/3.8/bioc/vignettes/msa/inst/doc/msa.pdf

# See also:
#browseVignettes(package = "Biostrings")

########## Some other packages that might be helpful for future development further analyses

## protr - R package for generating various numerical representation schemes of protein sequences 
##         e.g. biochemical properties.
##         Also has functions for retrieving AA sequences and 3D model data from Uniprot
## https://rdrr.io/cran/protr/f/vignettes/protr.Rmd

## UniprotR - looks similar to above
##            https://www.sciencedirect.com/science/article/pii/S1874391919303859?via%3Dihub

## Bio3d    - Compartive analysis of 3D structures
##            http://thegrantlab.org/bio3d/articles/online/pca_vignette/Bio3D_pca.html


########## Important! 
########## Currently the script does not account for length and indel differences between
########## the positions of features in the Uniprot ref and the alignment being analysed
########## (unless the alignment includes the species the feature info is drawn from as 
##########  as the longest sequence with no indels). 
########## Thus there may be discrepancies between data plotted from the alignment file
########## and the domain structure graph. Depending on the size/indel differences, these
########## discrepancies may be too small to notice, but BEWARE!!!
########## A future task is to add a function that can adjust feature start/end points
########## in the feature data frame generated from Uniprot data by the 
########## drawProteins::get_features() function.


########## Helper functions

########## Window sites - general function for summing a vector by sliding window ##########

window.sites<-function(x, y){                              ### x is the results vector to analyse, y is window size
  
  window_size      <- y
  num_positions    <- length(x)                            ### total residue positions in alignment
  num_windows      <- floor(num_positions / window_size)   ### number of windows of size 'y'
  num_polymorphic  <- numeric(num_windows)                 ### empty results vector number of poylmorphic sites within window 'i'
  
  # Sum the results vector in non-overlapping blocks of size window_size
  for (i in 1:num_windows) {
    start <- (i - 1) * window_size + 1
    end <- i * window_size
    num_polymorphic[i] <- sum(x[start:end])
  }
  return(num_polymorphic)  }

### End function


####### Function to generate the consensus sequence with random choices in the event of ties
generate_consensus <- function(alignment) {
  
  # Generate the consensus matrix
  cons_matrix <- consensusMatrix(alignment)
  
  # Initialize a vector for the consensus sequence
  consensus_seq <- character(ncol(cons_matrix))
  
  # Loop through each column in the consensus matrix
  for (i in 1:ncol(cons_matrix)) {
    # Get the amino acids with the highest counts
    max_count <- max(cons_matrix[, i])
    possible_aa <- rownames(cons_matrix)[cons_matrix[, i] == max_count]
    
    # If there is a tie (multiple amino acids with max count), choose one randomly
    if (length(possible_aa) > 1) {
      consensus_seq[i] <- sample(possible_aa, 1)  # Randomly select one
    } else {
      consensus_seq[i] <- possible_aa  # Single most frequent amino acid
    }
  }
  
  # Return the generated consensus sequence
  consensus <- paste(consensus_seq, collapse = "")
  return(consensus)
}


########## Amino acid string to character matrix ###########

###Needed as part of calculating numbers variable and fixed sites, 
###converts AAstringSet object to matrix of individual characters 

aas2cm<-function(y){t(sapply(y, function(x) unlist(strsplit(as.character(x), ""))))}

###Function end

########## Polymorphic sites ##################

## Calculates number of polymorphic sites in non-overlapping sliding window for an alignment object 
## x - AAstringSet object for alignment
## Set window size to 1 return vector with 1 recorded for each polymorphic site

polymorphicSites<-function(x, window_size){ 
  
  num_positions <- nchar(x[[1]])                        ### total residue positions in alignment 
  num_windows <- floor(num_positions / window_size)     ### number of windows of size 'x'
  num_polymorphic <- numeric(num_windows)               ### empty results vector number of poylmorphic sites within window 'x'
  results_vector <- rep(0, num_positions)               ### empty results vector for whether residue position 'x' is polymorphic (0,1)
  
  # Use the first sequence as the reference
  alignment_matrix <- aas2cm(x)
  reference_seq    <- alignment_matrix[1, ]
  
  # Iterate over each residue position
  for (i in 1:num_positions) {
    # Compare residue identities at position i with the reference sequence
    differences <- alignment_matrix[-1, i] != reference_seq[i]                     ### Here the '-1' drops the first sequence which is being used as the ref
    if (any(differences)) {
      results_vector[i] <- 1 }
  }
  
  if (window_size==1){return(results_vector)}
  else{  
    # Sum the results vector in non-overlapping blocks of size window_size
    for (i in 1:num_windows) {
      start <- (i - 1) * window_size + 1
      end <- i * window_size
      num_polymorphic[i] <- sum(results_vector[start:end])
    }
    
    return(num_polymorphic)}
  
}

###End function


########## Calculate positions of variable sites between single sequence pairs##############

seqdifs.single <- function(query.seq, ref.seq) {  #input sequences for query and reference (AAstrings)
  
  query.mat<-aas2cm(query.seq)  ##First convert to character matrix
  ref.mat  <-aas2cm(ref.seq)
  
  num_positions  <- nchar(query.seq[[1]])
  results_vector <- rep(0, num_positions)  ### empty results vector for whether residue position 'x' is polymorphic (0,1)
  
  # Iterate over each residue position
  for (i in 1:num_positions) {
    differences <- query.mat[i] != ref.mat[i]
    if (any(differences)) {
      results_vector[i] <- 1 }}
  return(results_vector)}

###End function


########## Calculate positions of variable sites for multiple alignment ####################
########## (Internal function for 'fixed.sites')  

seqdifs.multi <- function(query.seq, ref.seq) {        ### expects character matrix input not AAStringSets
  
  num_positions  <- sum(nchar(query.seq[1,]))          ### uses this syntax as working from character matrix not AAStringSet object
  results_vector <- rep(0, num_positions)              ### empty results vector for whether residue position 'x' is polymorphic (0,1)
  
  # Iterate over each residue position
  for (i in 1:num_positions) {
    differences <- query.seq[-1, i] != ref.seq[i]     ### Here the '-1' drops the first sequence which is being used as the ref
    if (any(differences)) {
      results_vector[i] <- 1 }}
  return(results_vector)}

###End function


########## Generate list of fixed sites between 2 specified groups ##############

#### Input is 2 different AAstringsets to calculate fixed differences between

## The logic behind this is:
## i)   First calculate the consensus sequences for the 2 sequence sets being compared
## ii)  Then create a vector for sites where these consensus sequences differ - sites with
##      fixed differences between groups will have different consensus states and be coded 1 
## iii) Next within each set identify polymorphic sites, and invert to give monomorphic sites
##      i.e. those that are fixed within each lineage
## iv)  Then for the focal lineage multiply the monomorphic sites vector, by the vector
##      of differences between the consensus sequences. Since everything is coded as 0 or 1
##      the monomorphic sites in the focal lineage that have the same state as the outgroup
##      will multiplied by 0, and revert to 0, leaving only monomorphic (fixed) sites in the 
##      in the focal lineage that differ in state from the outgroup coded as 1 (since they
##      would be multiplied by 1). Similarly sites that are different in state between the 
##      consensus sequences, but are polymorphic within groups will also be multiplied by 0,
##      reverting to 0 in the final results vector. 

fixed.sites <- function(x, y) {  
  # x = sequence set 1, y = sequence set 2
  
  ## First calculate consensus sequence comparison for the 2 sets
  set1.cons  <- generate_consensus(x)  # Use the new consensus function
  set2.cons  <- generate_consensus(y)  # Use the new consensus function
  
  # Convert to AAStringSet
  set1.cons  <- AAStringSet(set1.cons)
  set2.cons  <- AAStringSet(set2.cons)
  
  # Compute consensus sequence differences
  set1.set2.cons.difs  <- seqdifs.single(set2.cons, set1.cons)
  
  ## Monomorphic (fixed) sites in set 1
  set.1.mat  <- aas2cm(x)
  reference_seq.1  <- set.1.mat[1, ]
  set.1.difs  <- seqdifs.multi(set.1.mat, reference_seq.1)
  
  ## Need to invert differences so 1 represents monomorphic sites in the query group (1)
  set.1.monomorphic  <- (set.1.difs - 1) * -1
  
  ## Monomorphic (fixed) sites in set 2  
  set.2.mat  <- aas2cm(y)
  reference_seq.2  <- set.2.mat[1, ]
  set.2.difs  <- seqdifs.multi(set.2.mat, reference_seq.2)
  set.2.monomorphic  <- (set.2.difs - 1) * -1
  
  ########## Lastly multiply monomorphic and consensus diffs to get sites fixed in set 1 relative to set 2
  ###        This gives sites with strict fixed differences for alternate variants
  ###        i.e. sites which are exclusive for some polymorphic variants are not counted 
  
  fixed <- set1.set2.cons.difs * set.1.monomorphic
  fixed <- fixed * set.2.monomorphic  # Removes any sites polymorphic in set 2
  
  return(fixed) 
}

###End function
###End function


########## Calculate amino acid diversity ##############

########## Calculates AA div along length of AA sequence alignment 
#          for non-overlapping windows of 'y' residues

calc.AAdiv <- function(x, y) {   # x = sequence alignment object, y = window size
  
  ### Set up variables  
  
  num_positions         <- nchar(x[[1]])                ### total residue positions in alignment 'x'
  num_windows           <- floor(num_positions / y)     ### number of windows of size 'y'
  num_seqs              <- length(alignment_seq)        ### number of sequences in the alignment 
  amino_acid_diversity  <- numeric(num_windows)         ### variable for amino acid diversity
  ad1                   <- numeric(num_windows)          ###          for sum pairwise AA differences    
  
  ### AAdiv calculation
  
  for (i in 1:num_windows) {
    
    start <- (i - 1) * y + 1  ### specifies start/end point for non-overlapping sliding window... 
    end <- i * y                        ### ... at iteration 'i' 
    
    window_seq <- substr(x, start = start, stop = end)                                    ### subsets window 'i' from alignment
    # print(paste(c(i," ",i+y-1)))                                                        ### diagnostic (prints window start/end)  
    amino_acid_diversity[i]<-mean(stringDist(window_seq, method="hamming")/y)   ### calculates distance matrix (num of pw AA difs) and then calcs AA div for window'i'  
    ad1[i]                 <- sum(stringDist(window_seq, method="hamming"))               ### returns sum of pw difs for window 'i'
  }
  
  ### Prepares results data frame
  
  aadiv_df_out <- data.frame(
    window.num  = seq(1:num_windows),             #Window number
    residue.num  = seq(1:num_windows)*y, #gives the end residue position for window 'i' 
    aadiv        = amino_acid_diversity,           #amino acid diversity
    aadifs       = ad1                             #amino acid differences
  )
  return(aadiv_df_out)
  
}

### End function

########## Rescale vector #############

##Rescales a vector to a new range

rescale_vector <- function(x, new_min, new_max) {
  # Find the current minimum and maximum of the vector
  old_min <- min(x)
  old_max <- max(x)
  
  # Rescale the vector to the new range
  scaled_vector <- ((x - old_min) / (old_max - old_min)) * (new_max - new_min) + new_min
  
  return(scaled_vector)
}

### End function

########## Plot protein domains ################

## Basic plot of protein domain structures with nicer aesthetics compared to drawProtein
## Resulting ggplot object can then be augmented by theme adjustments and annotations
## Input data:
##           p.dat - A data frame of protein feature info generated by functions from the drawProtein package
##             x   - A negative integer value to set the minimum xlim value.
##                   This adds space to the left of the graph to add text annotations for 
##                   feature tracks such as active sites, fixed site locations, indels etc.
##                   Recommend specify 10-15% of total chain length. e.g. if protein 2500 aa,
##                   x would be -250 or -300. 
##             y's - The plot y lims are 0.5 to 1.5, with the plotting of features centred on y=1
##                   The 3 y values allow the relative heights on the chain, domains and regions
##                   to be adjusted

domainPlot<-function(p.dat, x, y.chain = 0.05 , y.domain = 0.125 , y.region = 0.1 ) {
  
  ## Plot the main chain
  
  dp<-ggplot() +
    
    geom_rect(    data = p.dat[which(p.dat$type=="CHAIN"),],
                  aes(xmin = begin, xmax = end, ymin = 1-y.chain, ymax = 1+y.chain)) + 
    xlab("Amino acid position")+
    xlim(x, max(p.dat$end)) + ylim(0.5, 1.5)
  
  ## Plot the domains
  
  if("DOMAIN" %in% as.factor(prot_data$type)==TRUE){
    dp<-dp+
      statebins:::geom_rrect(        ###uses statebins funtion to add rectangle with rounded corners
        data = p.dat[which(p.dat$type=="DOMAIN"),],
        aes(xmin = begin, xmax = end, ymin = 1-y.domain, ymax = 1+y.domain, fill = description),
        color = "darkgrey" ) }
  
  ## Plot the regions
  
  if("REGION" %in% as.factor(prot_data$type)==TRUE){
    dp<-dp+
      statebins:::geom_rrect(
        data = p.dat[which(p.dat$type=="REGION"),],
        aes(xmin = begin, xmax = end, ymin = 1-y.region, ymax = 1+y.region, fill = description),
        color = "darkgrey" )              }
  return(dp)
}

# End function

##########


########## Analysis

########## Read in data ##########

###Remember to check or update the taxonomy file/object so it matches the species composition and order of the alignment

taxonomy       <-read.csv("/path/to/DGAT1_taxonomy.csv")      ### taxonomy info
alignment      <-readAAStringSet("/path/to/dgat1_trimmed.fa")   ### read in FASN data (AA fasta file)
alignment_seq  <-AAStringSet(alignment)        ### make AAStringSet object from the alignment
setwd("/path/to/wd")

###If going to analyse a different gene later save copies

#lalba.taxonomy        <-taxonomy                      
#lalba.alignment_seq   <-alignment_seq

##### XDH  - skip if analysing FASN

#taxonomy       <-read.csv("taxonomy.csv")
#alignment      <-readAAStringSet("XDH.fas")   
#alignment_seq  <-AAStringSet(alignment)       

### Monachus monachus is missing from XDH so need to drop from taxonomy

#XDH.taxonomy<-taxonomy[-9,]                               ### drops row 9 with Monachus monachus
#XDH.taxonomy$seq.no<-seq(1:length(XDH.taxonomy$species))  ### Need to update seq nos due to row deletion
#taxonomy<-XDH.taxonomy

####

########## Define sequence groups ########

outgroup <- taxonomy$seq.no[taxonomy$tax4=="Outgroup"]  ### vector of sequence nos labelled as outgroup
pinnipedia <-taxonomy$seq.no[taxonomy$tax1=="Pinnipedia"]  ### same for other taxa
otariidae  <-taxonomy$seq.no[taxonomy$tax2=="Otariidae"]   ### Note '$tax2' column changes to switch to correct taxonmic level specification
phocidae   <-taxonomy$seq.no[taxonomy$tax2=="Phocidae"]
phocinae   <-taxonomy$seq.no[taxonomy$tax3=="Phocinae"]
monachinae <-taxonomy$seq.no[taxonomy$tax3=="Monachinae"]

outgroup.seqs    <-alignment_seq[c(outgroup)]              ### subsets to outgroup seqs
pinnipedia.seqs  <-alignment_seq[c(pinnipedia)]            ### subsets alignment to pinnipedia seqs
otariid.seqs     <-alignment_seq[c(otariidae)] 
phocid.seqs      <-alignment_seq[c(phocidae)]
phocinae.seqs    <-alignment_seq[c(phocinae)]
monachinae.seqs  <-alignment_seq[c(monachinae)]



#outgroup <- taxonomy_mcl1$seq.no[taxonomy$tax4=="Outgroup"]  ### vector of sequence nos labelled as outgroup
#pinnipedia <-taxonomy_mcl1$seq.no[taxonomy$tax1=="Pinnipedia"]  ### same for other taxa
#otariidae  <-taxonomy_mcl1$seq.no[taxonomy$tax2=="Otariidae"]   ### Note '$tax2' column changes to switch to correct taxonmic level specification
#phocidae   <-taxonomy_mcl1$seq.no[taxonomy$tax2=="Phocidae"]
#phocinae   <-taxonomy_mcl1$seq.no[taxonomy$tax3=="Phocinae"]
#monachinae <-taxonomy_mcl1$seq.no[taxonomy$tax3=="Monachinae"]

#outgroup.seqs    <-alignment_seq_mcl1[c(outgroup)]              ### subsets to outgroup seqs
#pinnipedia.seqs  <-alignment_seq_mcl1[c(pinnipedia)]            ### subsets alignment to pinnipedia seqs
#otariid.seqs     <-alignment_seq_mcl1[c(otariidae)] 
#phocid.seqs      <-alignment_seq_mcl1[c(phocidae)]
#phocinae.seqs    <-alignment_seq_mcl1[c(phocinae)]
#monachinae.seqs  <-alignment_seq_mcl1[c(monachinae)]
############


########## Calculate consensus sequences ##########

outgroup.cons <- generate_consensus(outgroup.seqs)
pinnipedia.cons <- generate_consensus(pinnipedia.seqs)
otariid.cons <- generate_consensus(otariid.seqs)
phocid.cons <- generate_consensus(phocid.seqs)
phocinae.cons <- generate_consensus(phocinae.seqs)
monachinae.cons <- generate_consensus(monachinae.seqs)


outgroup.cons    <-AAStringSet(outgroup.cons)        ### converts to AAstringSets 
pinnipedia.cons  <-AAStringSet(pinnipedia.cons)
otariid.cons     <-AAStringSet(otariid.cons)
phocid.cons      <-AAStringSet(phocid.cons)
phocinae.cons    <-AAStringSet(phocinae.cons)
monachinae.cons  <-AAStringSet(monachinae.cons)


consensus.seq    <-AAStringSet(c(outgroup.cons,
                                 pinnipedia.cons,
                                 otariid.cons,
                                 phocid.cons,
                                 phocinae.cons,
                                 monachinae.cons))   ### convert to an AAStringSet

names(consensus.seq) <-c( "Outgroup consensus",
                          "Pinniped consensus",
                          "Otariid consensus",
                          "Phocid consensus",
                          "Phocinae consensus",
                          "Monachinae consensus")   ### Adds names

consensus.seq


#lalba.consensus.seq <- consensus.seq
# FASN.consensus.seq<-consensus.seq    ### save for poserity
#  XDH.consensus.seq<-consensus.seq

sum(seqdifs.single(phocid.cons,otariid.cons))        ###Number of sites different between phocic & otariid consensus

which(seqdifs.single(phocid.cons,otariid.cons)==1)   ###List of sites with differences 
###Can cross reference this list with alignment view in MEGA

########## Calculating polymorphic and fixed sites ##########

## Pattern of polymorphic sites in the whole alignment

num_polymorphic<-polymorphicSites(alignment_seq, 1)  #No sliding window

sum(num_polymorphic)           ##gives total polymorphic sites
which(num_polymorphic==1)        ##lists positions of polymorphics sites

num_polymorphic<-polymorphicSites(alignment_seq, 20) # This time with non-overlapping sliding window of 20

# Create a data frame for the plot data

window_size   <- 20
num_positions <- nchar(alignment_seq[[1]])                ### total residue positions in alignment 
num_windows   <- floor(num_positions / window_size)

hist_df <- data.frame(position = seq(1:num_windows),             #Window number
                      res_num  = seq(1:num_windows)*window_size, #gives the end residue position for window 'i' 
                      count    = num_polymorphic                 #number of polymorphic sites in window 'i'
)
#XDH.polymorphic.df<-hist_df
#FASN.polymorphic.df<-hist_df                                     ### Keep for posterity 

# Plots of polymorphic residues (sliding window)

hist_plot.polysites <- ggplot(hist_df, aes(x = res_num, y = count)) +
  geom_bar(stat = "identity", fill = "lightblue" , color = "darkgrey", 
  ) +
  labs(x = "Amino acid position", y = "Number of polymorphic sites") +
  theme_light()  

hist_plot.polysites

#FASN.hist_plot.polysites<-hist_plot.polysites
# XDH.hist_plot.polysites<-hist_plot.polysites

line_plot.polysites <- ggplot(hist_df, aes(x = res_num, y = count)) +
  geom_line(color = "orange") +
  labs(x = "Amino acid position", y = "Number of polymorphic sites") +
  theme_light()
line_plot.polysites

#FASN.line_plot.polysites<-line_plot.polysites
# XDH.line_plot.polysites<-line_plot.polysites

hist_plot.polysites/line_plot.polysites  ###  for illustration: plot layout using syntax from the patchwork package 
hist_plot.polysites+line_plot.polysites

############ Pattern of fixed sites in different groups

### Note...
### This gives sites with strict fixed diffs for alternate variants between both groups
### i.e. Sites which are exclusive for some polymorphic variants are not counted here 

f.pinnipeds<-fixed.sites(pinnipedia.seqs, outgroup.cons)  ###fixed sites in all pinnipeds relative to outgroup
sum(f.pinnipeds)                                          ###numnber of sites and list
which(f.pinnipeds==1)

f.phocids<-fixed.sites(phocid.seqs, outgroup.cons)        ###Note this will also includes a subset of sites that are fixed in all pinnipeds
sum(f.phocids)
which(f.phocids==1)

f.otariids<-fixed.sites(otariid.seqs, outgroup.cons)
sum(f.otariids)
which(f.otariids==1)

### Can cross reference list of sites from 'which output' with MEGA alignment  

### This works but it includes sites that fixed between phocids and otarids, where one of the other also shares a variant present in the outgroup
### Therefore need an alternative to identify sites that are unique to phocids and or otariids, but distinct from outgroup
### Do this by specifying sequence sets that combine OG & either phocids or otariids, then compare with the other group 

og.seqs        <-taxonomy$seq.no[taxonomy$tax1!="Pinnipedia"]     ### vector of sequence nos labelled as outgroup
og.ots         <-append(og.seqs, taxonomy$seq.no[taxonomy$tax2=="Otariidae"])
og.phoc        <-append(og.seqs, taxonomy$seq.no[taxonomy$tax2=="Phocidae"])

og.ots.seqs      <-alignment_seq[c(og.ots)]                ### subsets
og.phoc.seqs     <-alignment_seq[c(og.phoc)]               ### subsets

unique.phocids   <-fixed.sites(phocid.seqs, og.ots.seqs)
unique.otariids  <-fixed.sites(otariid.seqs, og.phoc.seqs)

sum(unique.phocids)
which(unique.phocids==1)
sum(unique.otariids)
which(unique.otariids==1)

###works, but excludes sites where ots and phoc are fixed for different aas relative to og at same site

### Make a data frame of fixed sites results

fixed.df<-data.frame(f.pinnipeds,
                     f.phocids  = f.phocids-f.pinnipeds, ###Subtracting to remove sites fixed in all pinnipeds also fixed in phocids
                     f.otariids = f.otariids-f.pinnipeds,
                     unique.phocids,
                     unique.otariids)



########## Plotting of domain graphs #######################  
library(ggplot2)  
### Get the data for protein of interest from Uniprot

### P49327 is the Uniprot accession number for FASN, text must be in quotes

# prot_data <- drawProteins::get_features("P47989")         ### XDH
# prot_data <- drawProteins::get_features("P00709")         ### lalba
#prot_data <- drawProteins::get_features("Q13410")
#prot_data     <- drawProteins::feature_to_dataframe(prot_data)

# Fetch and format Uniprot feature data
prot_data <- drawProteins::get_features("O75907") # change to whatever Uniprot accession you need
prot_data <- drawProteins::feature_to_dataframe(prot_data)

library(dplyr)

# Copy data
adjusted_prot_data <- prot_data

# Convert AAStringSet to a character matrix (assuming 'alignment_seq' is already loaded)
alignment_matrix <- as.matrix(alignment_seq)


# Select human reference (row index as appropriate)
human_sequence <- alignment_matrix[1, ]  # 1 row vector of chars (length = ncol(alignment))

# Optional: if you want to ignore columns where human = '-' AND all pinnipeds = '-'
# (but the mapping below ignores human '-' anyway)
pinniped_sequences <- alignment_matrix[24:37, ]
# outgroup_insertions <- colSums(pinniped_sequences == "-") == nrow(pinniped_sequences) & human_sequence == "-"

# ---- Build mapping: human residue index -> alignment column ----
L <- length(human_sequence)
res_idx <- rep(NA_integer_, L)  # residue index at each column (NA if human gap)
r <- 0L
for (c in seq_len(L)) {
  aa <- human_sequence[c]
  if (!is.na(aa) && aa != "-") {
    r <- r + 1L
    res_idx[c] <- r
  }
}

# inverse map: for each human residue position (1..r), what's the alignment column?
max_res <- max(res_idx, na.rm = TRUE)
align_col_for_res <- rep(NA_integer_, max_res)
align_col_for_res[res_idx[!is.na(res_idx)]] <- which(!is.na(res_idx))

# Helper to convert a human (ungapped) feature coordinate to an alignment column
map_res_to_col <- function(pos_res) {
  if (is.na(pos_res)) return(NA_integer_)
  if (pos_res < 1L) return(NA_integer_)
  if (pos_res > length(align_col_for_res)) return(NA_integer_)
  align_col_for_res[pos_res]
}

# ---- Adjust features by lookup (no cumulative shift) ----
adjusted_prot_data <- prot_data

adj_begin <- integer(nrow(adjusted_prot_data))
adj_end   <- integer(nrow(adjusted_prot_data))

for (i in seq_len(nrow(adjusted_prot_data))) {
  fs <- adjusted_prot_data$begin[i]
  fe <- adjusted_prot_data$end[i]
  
  # Map human residue positions to alignment columns
  sc <- map_res_to_col(fs)
  ec <- map_res_to_col(fe)
  
  # If mapping failed (e.g., feature extends past known residues), try nearest non-NA
  if (is.na(sc)) {
    # nearest alignment column where res_idx == fs (should exist if fs <= max_res)
    if (!is.na(fs) && fs >= 1L && fs <= max_res) sc <- align_col_for_res[fs]
  }
  if (is.na(ec)) {
    if (!is.na(fe) && fe >= 1L && fe <= max_res) ec <- align_col_for_res[fe]
  }
  
  # Final safety: clamp and order
  if (is.na(sc)) sc <- 1L
  if (is.na(ec)) ec <- L
  if (sc > ec) { tmp <- sc; sc <- ec; ec <- tmp }
  
  adj_begin[i] <- max(1L, min(L, sc))
  adj_end[i]   <- max(1L, min(L, ec))
}

adjusted_prot_data$begin <- adj_begin
adjusted_prot_data$end   <- adj_end

# Optional: log changes
adjustment_log <- prot_data %>%
  mutate(Adjusted_Start = adjusted_prot_data$begin,
         Adjusted_End   = adjusted_prot_data$end) %>%
  select(description, begin, end, Adjusted_Start, Adjusted_End)




# Adjust the overall chain length
#adjusted_prot_data$begin[1] <- 1
#adjusted_prot_data$end[1] <- 600


##troubleshooting
#print(gap_indices)
#print(gap_counts)
#print(adjusted_prot_data$begin[i])
#print(adjusted_prot_data$end[i])

#thi_start <- adjusted_prot_data$begin[10]
#thi_end <- adjusted_prot_data$end[10]
#gap_counts_within_feature <- gap_counts[thi_start:thi_end]
#print(gap_counts_within_feature)


# prot_data.FASN <- prot_data    # save for posterity
#  prot_data.XDH  <-prot_data

#use 'fix(prot_data)' to view the data frame... remember to close the window!

domain_plot <-
  domainPlot(adjusted_prot_data, -70) +
  paletteer::scale_fill_paletteer_d("RColorBrewer::Set2") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = c(0.5, 0.9),
    legend.direction = "horizontal"
  ) +
  labs(fill = "")

domain_plot

########## Add short labels to domains, make a new data frame for convenience      

###For FASN, may not be necessary for proteins with simpler structures and can be skipped

#adjusted_prot_data1 <-adjusted_prot_data[6:9,]        #### selects rows 1:10 this depends on the number of features etc
#adjusted_prot_data1$short.lab<-c("extracellular",
                                 "extracellular",
                                 "XDH binding region",
                                 "intracellular")
#domain_plot<-
  
#  domain_plot+geom_label(data = adjusted_prot_data1[which(adjusted_prot_data1$type=="REGION"),], 
#                         aes(x=begin+(length/2), y=order, label=short.lab), size=3)+
#  geom_label(data = adjusted_prot_data1[which(adjusted_prot_data1$type=="DOMAIN"),], 
#             aes(x=begin+(length/2), y=order+0.025, label=short.lab), size=3)


#domain_plot


########## Annotating plot #################

#fixed.df<-FASN.fixed.df      ### Remember to select correct fixed sites df if working with multiple genes

###this plots fixed sites between pinnipeds and outgroup on gene structure 

###set up temporary data frames for each site class to be plotted

df1<-data.frame(begin = which(fixed.df$f.pinnipeds==1),
                order = 1)
df2<-data.frame(begin = which(fixed.df$f.phocids==1),
                order = 1)
df3<-data.frame(begin = which(fixed.df$unique.phocids==1),
                order = 1)
df4<-data.frame(begin = which(fixed.df$f.otariids==1),
                order = 1)
df5<-data.frame(begin = which(fixed.df$unique.otariids==1),
                order = 1)

### add some annotation lines for each site track

domain_plot <-  domain_plot +
  
  annotate("segment", x = 0, xend = max(adjusted_prot_data$end), y = adjusted_prot_data$order+0.15,  yend = adjusted_prot_data$order+0.15, color="grey")+
  annotate("segment", x = 0, xend = max(adjusted_prot_data$end), y = adjusted_prot_data$order+0.175, yend = adjusted_prot_data$order+0.175, color="grey")+
  annotate("segment", x = 0, xend = max(adjusted_prot_data$end), y = adjusted_prot_data$order+0.2,   yend = adjusted_prot_data$order+0.2, color="grey")+
  annotate("segment", x = 0, xend = max(adjusted_prot_data$end), y = adjusted_prot_data$order+0.225, yend = adjusted_prot_data$order+0.225, color="grey")
  #annotate("segment", x = 0, xend = max(adjusted_prot_data$end), y = adjusted_prot_data$order+0.25,  yend = adjusted_prot_data$order+0.25, color="grey")


### Now add points

domain_plot <-  domain_plot +           ###order+0.15 etc specifies vertical position relative to horizontal centre line of plot

  geom_point(data = df1, aes(x = begin, y = order+0.15),
  shape = 21, colour = "black", fill = "orange", size=3, alpha=0.5) +
  geom_point(data = df2, aes(x = begin, y = order+0.175),
  shape = 21, colour = "black", fill = "cyan", size=3, alpha=0.5) +
  geom_point(data = df4, aes(x = begin, y = order+0.2),
  shape = 21, colour = "black", fill = "darkcyan", size=3, alpha=0.5) +
  geom_point(data = df5, aes(x = begin, y = order+0.225),
  shape = 21, colour = "black", fill = "mediumorchid1", size=3, alpha=0.5)
  #geom_point(data = df5, aes(x = begin, y = order+0.25),
 # shape = 21, colour = "black", fill = "mediumorchid4", size=3, alpha=0.5) 

    ##add functional sites
  
  #annotate("segment", x = 0, xend = max(adjusted_prot_data$end), y = adjusted_prot_data$order+0.325, yend = adjusted_prot_data$order+0.325, color="grey")+
  #annotate("segment", x = 0, xend = max(adjusted_prot_data$end), y = adjusted_prot_data$order+0.3,  yend = adjusted_prot_data$order+0.3, color="grey")+
  
  
  #geom_point(data = adjusted_prot_data[adjusted_prot_data$type == "ACT_SITE",], aes(x = begin, y = order+0.325), 
  #shape = 24, colour = "black", fill = "lightblue", size=2.5)+
  #geom_point(data = adjusted_prot_data[adjusted_prot_data$type == "BINDING",],  aes(x = begin, y = order+0.3), 
  #shape = 23, colour = "black", fill = "lightgreen", size=2.5)
  
  
#rm(df4,df5) ### clean up  

domain_plot

domain_plot<-
  domain_plot + 
  
  
  annotate("text", x = -25, y = 1.15,  label= "Fixed Pinnipeds", size=3, hjust=1, vjust=0.5)+
  annotate("text", x = -25, y = 1.175,   label= "Fixed Phocids", size=3, hjust=1, vjust=0.5)+
  #annotate("text", x = -25, y = 1.20,  label= "Unique Phocids", size=3, hjust=1, vjust=0.5)+
  annotate("text", x = -25, y = 1.20,   label= "Fixed Otariids", size=3, hjust=1, vjust=0.5)+
  annotate("text", x = -25, y = 1.225,  label= "Unique Otariids", size=3, hjust=1, vjust=0.5)
  #annotate("text", x = -25, y = 1.30,  label= "Binding sites", size=3, hjust=1, vjust=0.5)+
  #annotate("text", x = -25, y = 1.325,    label=  "Active sites", size=3, hjust=1, vjust=0.5)

domain_plot


### If desired add a label for the gene symbol/name to the Chain (central bar) 
### Update 'label' to set the desired text, adjust the x value to set horizontal position 
### Can also add domain or region text/labels in the same way if desired

#domain_plot<-

# domain_plot+annotate("text", x = 250, y = 1, label = "BTN1A1",      
# size =5, color = "grey95" , vjust = 0.5)


#  domain_plot+annotate("text", x = 350, y = 1, label = "JAK1",      
                       size =5, color = "grey95" , vjust = 0.5)

### Plotting positions of deletions

### remember to select correct alignment if working with multiple genes

# alignment_seq<-FASN.alignment_seq

### use code like this to see which sequences have deletions
### need to look at all relevant sequences by changing the sequence no. in the second line
### Decide which sequences, if any, to plot deletions for
### use code like this to see which sequences have deletions
### need to look at all relevant sequences by changing the sequence no. in the bottom line
### Decide which sequences, if any, to plot deletions for

### Convert AAStringSet to character matrix
alignment_matrix <- aas2cm(alignment_seq)

# Identify columns where **all sequences except human** have gaps
all_except_human_gaps <- colSums(alignment_matrix[-1, ] == "-") == (nrow(alignment_matrix) - 1)

# Filter only **pinniped sequences** (rows 5-20)
pinniped_seqs <- alignment_matrix[26:39, ]

# Identify **true deletions** in pinnipeds
true_deletions <- matrix(FALSE, nrow = nrow(pinniped_seqs), ncol = ncol(pinniped_seqs))

for (i in seq_len(nrow(pinniped_seqs))) {
  for (j in seq_len(ncol(pinniped_seqs))) {
    # A "true deletion" occurs when:
    # 1. The pinniped sequence has a gap ("-")
    # 2. The equivalent position in the human reference **isn't** a gap
    # 3. The position is NOT an alignment artifact (not in `all_except_human_gaps`)
    if (pinniped_seqs[i, j] == "-" && alignment_matrix[1, j] != "-" && !all_except_human_gaps[j]) {
      true_deletions[i, j] <- TRUE
    }
  }
}

# Store the filtered deletion matrix
alignment_matrix_filtered <- alignment_matrix  # Copy original matrix

# Apply filtering only for **pinniped sequences**
alignment_matrix_filtered[24:37, ] <- ifelse(true_deletions, "-", "X")

# Find positions of **true** deletions in any pinniped species
which(alignment_matrix_filtered[37,] == "-")  # Example for sequence 9

# Vector of 0/1 for gaps in human
human_gap_counts <- as.numeric(alignment_matrix[1, ] == "-")

# Build map: original alignment positions to adjusted coords (i.e., only non-gap sites in human)
original_to_adjusted <- rep(NA, length(human_gap_counts))
adjusted_pos <- 1

for (i in seq_along(human_gap_counts)) {
  if (human_gap_counts[i] == 0) {
    original_to_adjusted[i] <- adjusted_pos
    adjusted_pos <- adjusted_pos + 1
  }
}

# Function to map deletion positions into adjusted space
map_deletions <- function(deletion_positions, mapping_vector, order_val) {
  # Map deletion positions (some may be NA if they hit gaps)
  adjusted <- mapping_vector[deletion_positions]
  adjusted <- adjusted[!is.na(adjusted)]  # remove unmappable positions
  data.frame(begin = adjusted, order = order_val)
}

# Example for 5 pinniped species (adjust the row numbers accordingly)
df1 <- map_deletions(which(alignment_matrix_filtered[26, ] == "-"), original_to_adjusted, 1)     # walrus
df2 <- map_deletions(which(alignment_matrix_filtered[39, ] == "-"), original_to_adjusted, 1)     # hooded
#df3 <- map_deletions(which(alignment_matrix_filtered[34, ] == "-"), original_to_adjusted, 1)     # hooded
#df4 <- map_deletions(which(alignment_matrix_filtered[28, ] == "-"), original_to_adjusted, 1)     # hooded seal
#df5 <- map_deletions(which(alignment_matrix_filtered[30, ] == "-"), original_to_adjusted, 1)     # grey seal
#df6 <- map_deletions(which(alignment_matrix_filtered[28, ] == "-"), original_to_adjusted, 1)     # nes
#df7 <- map_deletions(which(alignment_matrix_filtered[30, ] == "-"), original_to_adjusted, 1)     # ses
# Create an empty data frame to hold results

domain_plot <- domain_plot +
  geom_point(data = df1, aes(x = begin, y = order - 0.15), colour = "skyblue", size = 1.5) +
  geom_point(data = df2, aes(x = begin, y = order - 0.175), colour = "hotpink", size = 1.5)+
  #geom_point(data = df3, aes(x = begin, y = order - 0.2), colour = "hotpink", size = 1.5)+
  
  

  

    
  annotate("text", x = -20, y = 1 - 0.15,   label = "Deletion Pinnipeds",            size = 2.7, hjust = 1) +
  annotate("text", x = -20, y = 1 - 0.175,  label = "Deletion hooded seal",   size = 2.7, hjust = 1)
  #annotate("text", x = -20, y = 1 - 0.2,  label = "Deletion hooded seal",   size = 2.7, hjust = 1)



domain_plot


rm(df1,df2,df3,df4,df5,df6,df7) ### clean up 


#### Add amino acid diversity track
#BiocManager::install("pwalign")

aadiv_df_pins<-calc.AAdiv(pinnipedia.seqs, 10)        ### AA diversity pinnipeds only, 10 AA window
head(aadiv_df_pins)                                   ### take a quick look

### Make a plot 

p_aadiv_pins   <-ggplot(aadiv_df_pins, aes(x = residue.num, y = aadiv)) +
  geom_line(color = "orange") +
  xlim(0,max(adjusted_prot_data$end))+
  labs(x = "Amino acid position", y = "Amino acid diversity") +
  theme_light()

p_aadiv_pins

#### Now add the AA diversity track to the plot

# Set diversity cap
aadiv_cap <- 0.25
y_offset <- 0.525

domain_plot <- domain_plot +
  
  ### Shaded rectangle for 95% quantile (capped at 0.17)
  annotate("rect", 
           xmin = 0, 
           xmax = max(aadiv_df_pins$residue.num), 
           ymin = min(quantile(aadiv_df_pins$aadiv, probs = 0.025), aadiv_cap) + y_offset, 
           ymax = min(quantile(aadiv_df_pins$aadiv, probs = 0.975), aadiv_cap) + y_offset, 
           fill = "lightgrey", alpha = 0.2) +
  
  ### The AA diversity line (capped)
  geom_line(data = aadiv_df_pins, 
            aes(x = residue.num, y = pmin(aadiv, aadiv_cap) + y_offset), 
            color = "orange") +
  
  ### Axes and threshold lines
  annotate("segment", x = 0, xend = max(aadiv_df_pins$residue.num), y = 0 + y_offset, yend = 0 + y_offset, color = "darkgrey") +
  annotate("segment", x = 0, xend = 0, 
           y = 0 + y_offset, 
           yend = aadiv_cap + y_offset, 
           color = "darkgrey") +
  annotate("segment", x = 0, xend = max(aadiv_df_pins$residue.num), 
           y = min(mean(aadiv_df_pins$aadiv), aadiv_cap) + y_offset, 
           yend = min(mean(aadiv_df_pins$aadiv), aadiv_cap) + y_offset,
           color = "darkgrey", linetype = "dashed") +
  
  ### Axis labels (capped at 0.17)
  annotate("text", x = -20, y = 0 + y_offset, label = "0", size = 2.5, hjust = 1) +
  annotate("text", x = -20, y = aadiv_cap + y_offset, label = aadiv_cap, size = 2.5, hjust = 1) +
  annotate("text", x = -20, y = (aadiv_cap / 2) + y_offset, label = round(aadiv_cap / 2, 2), size = 2.5, hjust = 1) +
  
  ### Mean label
  annotate("text", x = 50, 
           y = (aadiv_cap * 0.85) + y_offset, 
           label = paste("Mean =", round(mean(aadiv_df_pins$aadiv), 4)), 
           size = 2.5, hjust = 0) +
  
  ### Vertical label
  annotate("text", x = -110, 
           y = (aadiv_cap / 2) + y_offset, 
           label = "Amino acid diversity", size = 2.5, angle = 90)

domain_plot

# XDH.plot<-domain_plot

# FASN.plot/XDH.plot     # example of stacked plot layout using patchwork package syntax

########## Plotting alignments with ggmas() #############

### Nice msa of focal region with all kinds of cool options
### http://yulab-smu.top/ggmsa/index.html

available_msa() # tells which type of data objects can be plotted

### This makes some alignment plots for the FASN Acyl Carrier Protein region (ACP)
### In the current test alignment this lies at residues ~2076-2153
### But be sure to identify the right region, e.g. search for motifs in MEGA to help
### find start-end positions. You can get the sequence of specific regions/domains in Uniprot

### Note the plotting of this can be very, very, very slow!!! So wait... can be >2mins
### Also can take a while switching between plots in RStudio



BTN1A1.acp.alignment.vis<-ggmsa(alignment_seq, 450, 510, seq_name = TRUE, char_width = 0.5) + 
  geom_seqlogo(color = "Chemistry_AA") + geom_msaBar()

BTN1A1.acp.alignment.vis

### Other optons with consensus view and shading for hypdrophobicity of residues

### Consensus view

BTN1A1.acp.alignment.vis.cons<-
  
  ggmsa(alignment_seq, 450, 510, seq_name = TRUE, char_width = 0.5, 
        consensus_views=TRUE , use_dot=TRUE) + 
  geom_seqlogo(color = "Chemistry_AA") + geom_msaBar()

### hydrophobicity colouring

ggmsa(FASN_sequence_with_indels, 2076, 2153, color = "Hydrophobicity", 
      seq_name = TRUE, char_width = 0.5) + 
  geom_seqlogo(color = "Hydrophobicity") + geom_msaBar()









#
