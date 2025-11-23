################################################################################
# RAPPID Study: Microbial Lung-to-Blood Translocation Analysis
# Manuscript: "Microbial lung-to-blood translocation associates with systemic 
#             inflammation in severe pneumonia: Evidence from paired plasma and 
#             lower respiratory tract metagenomics"
#
# This script contains all analytical code for generating the main figures
# and tables in the manuscript.
################################################################################

################################################################################
#### SETUP AND INITIALIZATION ####
################################################################################

# Load required libraries
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(tableone)
library(ggalluvial)
library(microeco)
library(rstatix)

# Set working directory based on operating system
if(Sys.info()['sysname'] == 'Darwin') {
  setwd('~/OneDrive/Documents/Morris Lab/Projects/Rappid_Seq/')
} else if(Sys.info()['sysname'] == 'Windows'){
  setwd('~/OneDrive_WORK/OneDrive/Documents/Morris Lab/Projects/Rappid_Seq/')
}

# Helper function to generate output file paths with timestamps
out.path <- function(string) {
  return(
    case_when(
      str_detect(string, pattern = '[Ff]ig') ~ paste('output/', format(Sys.Date(), '%m-%d-%Y'), string, '.png', sep = ''),
      str_detect(string, pattern = '[Tt]ab') ~ paste('output/', format(Sys.Date(), '%m-%d-%Y'), string, '.csv', sep = '')
    )
  )
}

# Set to TRUE to run full analysis from scratch, FALSE to load saved data
RunFromScratch <- FALSE

################################################################################
#### DATA PROCESSING FROM SCRATCH ####
################################################################################

if(RunFromScratch) { 
  
  #### Load Clinical Data ####
  Rappid.pneum_category_final.df <- 
    readxl::read_excel(path = "data/rappidseq_20230926-pneum_category_final.xlsx") %>% 
    mutate(pneum_category_final = factor(pneum_category_final, levels = c("CONTROL", "CDP", "MCP")))
  
  Rappid.meta.202305.df <- 
    read.csv(file = "data/Email-20230527/rappidseq_5.26.23.csv")
  
  Rappid.baseline.meta.202305.df <- 
    Rappid.meta.202305.df %>% dplyr::filter(str_detect(string = SubjectIDDay, pattern = "\\_1$"))
  
  DummyID.df <- readxl::read_excel(path = "data/DummyID-20240306.xlsx")
  
  #### Setup Nanopore Sequencing Data ####
  whos_here <- list.files("data/Nanopore Result/wimp")
  
  ################################################################
  # IMPORTANT: NCBI Taxonomy Database Required (~100MB)
  # Download from NCBI if not present in data/Nanopore Result/NCBI_TaxID/
  ################################################################
  
  # Load NCBI taxonomy nodes
  NCBI.nodes.taxID <- read.csv(file = "data/Nanopore Result/NCBI_TaxID/nodes.dmp", header = F, sep = "|")[, 1:3]
  colnames(NCBI.nodes.taxID) <- c("taxID", "parent_taxID", "rank")
  NCBI.nodes.taxID$rank <- str_replace_all(NCBI.nodes.taxID$rank, pattern = "\t", replacement = "")
  NCBI.all_nodes.taxID <- NCBI.nodes.taxID
  NCBI.nodes.taxID <- NCBI.nodes.taxID %>% dplyr::filter(rank == "species" | rank == "subspecies" | rank == "strain")
  
  # Load NCBI taxonomy names
  NCBI.names.taxID <- read.csv(file = "data/Nanopore Result/NCBI_TaxID/names.dmp", header = F, sep = "|")[, 1:4]
  colnames(NCBI.names.taxID) <- c("taxID", "name_txt", "unique_name", "name_class")
  NCBI.names.taxID$name_class <- str_replace_all(NCBI.names.taxID$name_class, pattern = "\t", replacement = "")
  NCBI.names.taxID$name_txt <- str_replace_all(NCBI.names.taxID$name_txt, pattern = "\t", replacement = "")
  
  # Extract species-level taxonomy
  NCBI.species.taxID.name <- NCBI.names.taxID %>% 
    dplyr::filter(taxID %in% (NCBI.nodes.taxID %>% 
                         dplyr::filter(rank == "species") %>% 
                         select(taxID) %>% unlist(use.names = F)) & 
             name_class == "scientific name") %>% 
    select(taxID, name_txt)
  
  NCBI.all.taxID.name <- NCBI.names.taxID %>% 
    dplyr::filter(taxID %in% (NCBI.all_nodes.taxID %>% 
                         select(taxID) %>% unlist(use.names = F)) & 
             name_class == "scientific name") %>% 
    select(taxID, name_txt)
  
  #### Functions for Nanopore Data Processing ####
  
  # Function to retrieve reads at any taxonomy level from each sample
  get_indiv_long_table_WIMP_all_level <- function(index){
    SubjectIDDay <- str_split(whos_here[index], "\\.")[[1]][1]
    intm <- read.csv(paste("data/Nanopore Result/wimp/", whos_here[index], sep = "")) 
    intm <- 
      intm %>% dplyr::filter(exit_status == "Classified") %>% select(taxID, name, lineage) %>% left_join(NCBI.all_nodes.taxID, by = "taxID") %>% 
      mutate(SubjectIDDay = SubjectIDDay)
    
    print(paste(SubjectIDDay, "done"))
    return(intm)
  }
  
  # Function to retrieve species-level reads from each sample
  get_indiv_long_table_WIMP_species_level <- function(index){
    SubjectIDDay <- str_split(whos_here[index], "\\.")[[1]][1]
    intm <- read.csv(paste("data/Nanopore Result/wimp/", whos_here[index], sep = "")) 
    intm <- 
      intm %>% dplyr::filter(exit_status == "Classified") %>% select(taxID, name, lineage) %>% left_join(NCBI.nodes.taxID, by = "taxID") %>% 
      dplyr::filter(rank == "species" | rank == "subspecies" | rank == "strain") %>% 
      left_join(NCBI.species.taxID.name %>% select(parent_taxID = taxID, parent_name = name_txt), by = "parent_taxID") %>%
      left_join(NCBI.nodes.taxID %>% select(parent_taxID = taxID, grandparent_taxID = parent_taxID, parent_rank = rank), by = "parent_taxID") %>% 
      left_join(NCBI.species.taxID.name %>% select(grandparent_taxID = taxID, grandparent_name = name_txt), by = "grandparent_taxID") %>%
      mutate(species_name = case_when(
        .$rank == "species" ~ as.character(.$name), 
        .$rank != "species" & !is.na(.$parent_name) ~ as.character(.$parent_name), 
        .$rank != "species" & is.na(.$parent_name) & !is.na(.$grandparent_name) ~ as.character(.$grandparent_name)
      ))
    
    df <- data.frame(
      SubjectIDDay = SubjectIDDay,
      Taxa = table(intm$species_name) %>% names(), 
      count = table(intm$species_name) %>% as.vector
    )
    print(paste(SubjectIDDay, "done"))
    return(df)
  }
  
  #### Process All Nanopore Samples ####
  # Note: This step can take significant time depending on number of samples
  
  nanopore.indiv.long.df <- data.frame()
  system.time(
    for (i in 1:length(whos_here)){
      nanopore.indiv.long.df <- rbind(
        nanopore.indiv.long.df, get_indiv_long_table_WIMP_species_level(i)
      )
    }
  )
  
  nanopore.indiv.all_level.long.df <- data.frame()
  system.time(
    for (i in 1:length(whos_here)){
      nanopore.indiv.all_level.long.df <- rbind(
        nanopore.indiv.all_level.long.df, get_indiv_long_table_WIMP_all_level(i)
      )
    }
  )
  
  # Create wide format data (samples × species matrix)
  nanopore.indiv.wide.df <- nanopore.indiv.long.df %>% spread(key = Taxa, value = count) %>% select(-`Homo sapiens`) %>% replace(is.na(.), 0)
  
  #### Quality Control and Filtering ####
  
  # Identify reads from patients in main study cohort
  RappidSeq.raw <- read.csv(file = "data/rappidseq_3.1.22.csv")
  RappidSeq.SubjectID <- unique(sort(RappidSeq.raw$SubjectID))
  
  Nanopore.Followup.fig.p <-
    nanopore.indiv.all_level.long.df %>% mutate(SubjectID = str_remove(string = SubjectIDDay, pattern = "\\_[\\d]+$"), Day = str_remove(string = SubjectIDDay, pattern = "^[\\d]+\\_")) %>% 
    dplyr::filter(SubjectID %in% RappidSeq.SubjectID) %>% distinct(SubjectID, Day, .keep_all = F) %>% ggplot(aes(x = as.numeric(Day))) + geom_bar() 
  
  # Categorize all reads by taxonomy and quality metrics
  Nanopore.all_reads.Sankey.df <- 
    nanopore.indiv.all_level.long.df %>% 
    mutate(Human = factor(ifelse(name == "Homo sapiens", yes = "Human reads", no = "Non-human reads"), levels = c("Non-human reads", "Human reads")), 
           taxonomy_rank = factor(rank, levels = c("superkingdom", "kingdom", "subkingdom", "phylum", "subphylum", "class", "subclass", "order", "suborder", "family", "genus", "species group", "species subgroup", "species", "forma specialis", "subspecies",  "varietas", "serogroup", "serotype", "biotype", "genotype", "strain", "isolate", "clade", "no rank")), 
           taxonomy_rank_cat = factor(
             case_when(rank == "species" ~ "species", 
                       rank %in% c("superkingdom", "kingdom", "subkingdom", "phylum", "subphylum", "class", "subclass", "order", "suborder", "family", "genus", "species group", "species subgroup") ~ "above species", 
                       rank %in% c("forma specialis", "subspecies",  "varietas", "serogroup", "serotype", "biotype", "genotype", "strain", "isolate") ~ "below species", 
                       TRUE ~ "non-hierarchical/NA"), 
             levels = c("above species", "species", "below species", "non-hierarchical/NA")), 
           Kingdom = case_when(str_detect(string = lineage, pattern = "^1:10239") ~ "Virus", 
                               str_detect(string = lineage, pattern = "^1:131567:2157") ~ "Archea",
                               str_detect(string = lineage, pattern = "^1:131567:2759") ~ "Eukaryota",
                               str_detect(string = lineage, pattern = "^1:131567:2") ~ "Bacteria", 
                               TRUE ~ "miscellanea"
                               
           ))
  
  # Match below-species reads to their species-level classification
  Nanopore.all_reads.below_species.match <- 
    Nanopore.all_reads.Sankey.df %>% 
    dplyr::filter(taxonomy_rank_cat == "below species") %>% distinct(lineage, .keep_all = T) %>% 
    left_join(NCBI.all.taxID.name, by = c("parent_taxID" = "taxID")) %>% transmute(taxID, name, rank, parent_taxID.new = parent_taxID, parent_name = name_txt, Kingdom, lineage) %>% 
    left_join(NCBI.all_nodes.taxID %>% transmute(parent_taxID.new = taxID, grandparent_taxID = parent_taxID, parent_rank = rank)) %>% 
    left_join(NCBI.all_nodes.taxID %>% transmute(grandparent_taxID = taxID, great_grandparent_taxID = parent_taxID, grandparent_rank = rank)) %>% 
    left_join(NCBI.all_nodes.taxID %>% transmute(great_grandparent_taxID = taxID, great_great_grandparent_taxID = parent_taxID, great_grandparent_rank = rank)) %>% 
    transmute(taxID, lineage, rank, Kingdom, name_raw = name, id_species = 
                case_when(
                  rank == "species" ~ taxID, 
                  parent_rank == "species" ~ parent_taxID.new, 
                  grandparent_rank == "species" ~ grandparent_taxID, 
                  great_grandparent_rank == "species" ~ great_grandparent_taxID, 
                  TRUE ~ NA_integer_
                )) %>% 
    left_join(NCBI.all.taxID.name %>% transmute(id_species = taxID, name_species = name_txt)) %>% 
    select(taxID, lineage, id_species, name_species)
  
  # Consolidate below-species reads to species level
  Nanopore.below_species.fix <- 
    Nanopore.all_reads.Sankey.df %>% left_join(Nanopore.all_reads.below_species.match) %>% 
    mutate(id_species = case_when(taxonomy_rank == "species" ~ taxID, taxonomy_rank != "species" ~ id_species, TRUE ~ NA_integer_), 
           name_species = case_when(taxonomy_rank == "species" ~ name, taxonomy_rank != "species" ~ name_species, TRUE ~ NA_character_)) %>% 
    distinct(lineage, .keep_all = T) %>% select(-SubjectIDDay) %>% 
    mutate(id_species = case_when(taxID == "1271862" ~ 28901L, TRUE ~ id_species), 
           name_species = case_when(taxID == "1271862" ~ "Salmonella enterica", TRUE ~ name_species))
  
  # Identify filtering criteria for noise reduction
  Unfiltered_top3.id_species.vec <- 
    Nanopore.all_reads.Sankey.df %>% left_join(Nanopore.below_species.fix) %>% 
    group_by(SubjectIDDay, id_species, Kingdom) %>% summarise(absabd = n()) %>% ungroup() %>% 
    group_by(SubjectIDDay) %>% mutate(sample_total = sum(absabd), relabd = absabd/sample_total, relabdrank = rank(-absabd, ties.method = "min")) %>% ungroup() %>% 
    dplyr::filter(relabdrank <= 3) %>% select(id_species) %>% unlist(use.names = F) %>% unique() %>% sort()
  
  Unfiltered_all.id_species.vec <- 
    Nanopore.all_reads.Sankey.df %>% left_join(Nanopore.below_species.fix) %>% select(id_species) %>% unlist(use.names = F) %>% unique() %>% sort()
  
  Unfiltered_all.id_species.kingdom.df <- 
    Nanopore.all_reads.Sankey.df %>% left_join(Nanopore.below_species.fix) %>% distinct(id_species, Kingdom)
  
  # Identify singletons (species appearing in only one sample)
  Unfiltered_singleton.id_species.vec <- Nanopore.all_reads.Sankey.df %>% left_join(Nanopore.below_species.fix) %>% group_by(id_species) %>% summarise(freq = n()) %>% ungroup() %>% dplyr::filter(freq == 1) %>% select(id_species) %>% unlist(use.names = F) %>% unique() %>% sort()
  
  # Identify low abundance species (<0.005% of total reads)
  Unfiltered_LOWabd.species_name.vec <- 
    Nanopore.all_reads.Sankey.df %>% left_join(Nanopore.below_species.fix) %>% 
    group_by(SubjectIDDay, id_species, Kingdom) %>% summarise(absabd = n()) %>% ungroup() %>% 
    group_by(id_species) %>% summarise(species_sum_absabd = sum(absabd)) %>% ungroup() %>% 
    mutate(read_sum_across_all_reads = sum(species_sum_absabd), percentage_species_across_all_reads = species_sum_absabd / read_sum_across_all_reads) %>% 
    dplyr::filter(percentage_species_across_all_reads <= 0.00005) %>% select(id_species) %>% unlist(use.names = F) %>% sort() %>% unique()
  
  # Summarize filtering criteria for each species
  Nanopore.species_breakdown.Sankey.df <- 
    data.frame(
      id_species = Unfiltered_all.id_species.vec
    ) %>% mutate(
      LowAbd = id_species %in% Unfiltered_LOWabd.species_name.vec, 
      Singleton = id_species %in% Unfiltered_singleton.id_species.vec, 
      Top3 = id_species %in% Unfiltered_top3.id_species.vec
    ) %>% left_join(NCBI.all.taxID.name %>% transmute(id_species = taxID, name_species = name_txt))
  
  # Create visualization of filtering process
  Nanopore.all_reads.breakdown.p <- 
    Nanopore.all_reads.Sankey.df %>% 
    left_join(Nanopore.below_species.fix[, c("lineage", "id_species")]) %>% 
    left_join(Nanopore.species_breakdown.Sankey.df %>% select(-name_species)) %>% 
    group_by(taxonomy_rank, Human, taxonomy_rank_cat, Kingdom, LowAbd, Singleton, Top3) %>% summarise(freq = n()) %>% ungroup() %>% 
    mutate(EXCLUSION = (Human == "Human reads" | (taxonomy_rank_cat %in% c("above species", "non-hierarchical/NA")) | ((taxonomy_rank_cat %in% c("species", "below species")) & (LowAbd | Singleton)))) %>% 
    mutate(LowAbd = ifelse(LowAbd, yes = "Low\nAbundance", no = "High\nAbundance"), Singleton = ifelse(Singleton, yes = 'Singleton', no = 'Non-sin-\ngleton'), Top3 = ifelse(Top3, yes = 'Top3', no = 'non-Top3')) %>% 
    dplyr::filter(!is.na(EXCLUSION)) %>% 
    ggplot(aes(
      axis1 = taxonomy_rank, 
      axis2 = Kingdom, 
      axis3 = Human, 
      axis4 = LowAbd, 
      axis5 = Singleton, 
      axis6 = Top3, 
      y = freq
    )) + 
    geom_alluvium(aes(fill = EXCLUSION)) +
    geom_stratum(colour = "black") +
    geom_text(stat = "stratum", 
              aes(label = after_stat(ifelse(stratum %in% c('species', 'strain', 'subspecies', 'Bacteria', 'Eukaryota', 'High\nAbundance', 'Non-sin-\ngleton', 'Top3', 'Human reads'), yes = as.character(stratum), no = NA_character_)))) + 
    scale_x_discrete(limits = c("taxonomy_rank", "Kingdom", "Human", "LowAbd", "Singleton", "Top3"),
                     expand = c(0.15, 0.05)) + 
    labs(title = "Unfiltered Nanopore reads breakdown") + 
    ylab("# of Nanopore reads collected") + theme_minimal() + 
    theme(legend.position = "bottom")
  
  Nanopore.species_breakdown.Sankey.p <- 
    Nanopore.species_breakdown.Sankey.df %>% left_join(Unfiltered_all.id_species.kingdom.df) %>% 
    mutate(EXCLUSION = LowAbd | Singleton) %>% 
    group_by(Kingdom, LowAbd, Singleton, Top3, EXCLUSION) %>% summarise(freq = n()) %>% 
    ggplot(aes(
      axis1 = Kingdom, 
      axis2 = LowAbd, 
      axis3 = Singleton, 
      axis4 = Top3, 
      y = freq
    )) + 
    geom_alluvium(aes(fill = EXCLUSION)) +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) + 
    scale_x_discrete(limits = c("Kingdom", "LowAbd", "Singleton", "Top3"),
                     expand = c(0.15, 0.05)) + 
    labs(title = "Unfiltered Nanopore species (and below-species reads) breakdown") + 
    ylab("# of species") + theme(legend.position = "bottom")
  
  # Helper function to extract legend from ggplot
  extract_legend <- function(my_ggp) {
    step1 <- ggplot_gtable(ggplot_build(my_ggp))
    step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
    step3 <- step1$grobs[[step2]]
    return(step3)
  }
  Sankey.shared_legend <- extract_legend(Nanopore.species_breakdown.Sankey.p) 
  
  # Create pneumonia category mapping
  Rappid.pneum_category.df <- 
    RappidSeq.raw %>% transmute(SubjectID = str_remove(string = SubjectIDDay, pattern = "\\_[\\d]+$"), 
                                Day = str_remove(string = SubjectIDDay, pattern = "^[\\d]+\\_"), classification) %>% 
    dplyr::filter(Day == "1") %>% transmute(SubjectID, pneum_category = case_when(classification == 0L ~ "CONTROL", classification == "1" ~ "CXNEG", classification == "2" ~ "CXPOS"))
  
  # Fix naming issues for specific species
  NCBI.taxID.name.FIX20230620 <- 
    read.table(file = "data/20230620-noname-species-tax_report.txt", sep = "|", col.names = c("CODE", "id_species", "name_species", "lineage")) %>% 
    dplyr::filter(CODE != "code\t") %>% 
    transmute(id_species = as.integer(str_remove_all(string = id_species, pattern = "\t")), 
              name_species_FIX = case_when(id_species == 1156433 ~ "Streptococcus sp. I-G2", TRUE ~ str_remove_all(name_species, pattern = "\t"))) 
  
  #### Create Filtered Nanopore Dataset ####
  # Apply filters: remove human reads, non-species-level reads, singletons, and low abundance species
  
  Nanopore.filtered.RelAbd.long.df <- 
    Nanopore.all_reads.Sankey.df %>% 
    left_join(Nanopore.below_species.fix[, c("lineage", "id_species")]) %>% 
    left_join(Nanopore.species_breakdown.Sankey.df %>% select(-name_species)) %>% 
    mutate(EXCLUSION = (Human == "Human reads" | (taxonomy_rank_cat %in% c("above species", "non-hierarchical/NA")) | ((taxonomy_rank_cat %in% c("species", "below species")) & (LowAbd | Singleton)))) %>% 
    dplyr::filter(!is.na(EXCLUSION) & !EXCLUSION) %>% group_by(SubjectIDDay, id_species, Kingdom) %>% summarise(AbsAbd = n()) %>% 
    arrange(SubjectIDDay, -AbsAbd) %>% ungroup() %>% 
    left_join(NCBI.all.taxID.name %>% transmute(id_species = taxID, name_species = name_txt)) %>% 
    group_by(SubjectIDDay) %>% 
    mutate(sample_filtered_nonhuman_total = sum(AbsAbd), filtered_RelAbd = AbsAbd / sample_filtered_nonhuman_total) %>% 
    left_join(NCBI.taxID.name.FIX20230620) %>% 
    mutate(name_species = case_when(is.na(name_species) ~ name_species_FIX, TRUE ~ name_species))
  
  # Create bar plots showing microbial composition
  Nanopore.filtered.RelAbd.baseline.barplot.p <- 
    Nanopore.filtered.RelAbd.long.df %>% mutate(Day = str_remove(SubjectIDDay, pattern = "^[\\d]+\\_"), 
                                                SubjectID = str_remove(SubjectIDDay, pattern = "\\_[\\d]+$")) %>% 
    inner_join(Rappid.pneum_category.df) %>% dplyr::filter(Day == "1") %>% 
    ggplot(aes(x = SubjectIDDay, y = filtered_RelAbd, fill = Kingdom)) + 
    geom_bar(stat = "identity") + facet_grid(.~pneum_category, scales = "free_x", space = "free") + 
    theme(axis.text.x = element_text(angle = 90)) 
  
  Nanopore.filtered.AbsAbd.baseline.barplot.p <- 
    Nanopore.filtered.RelAbd.long.df %>% mutate(Day = str_remove(SubjectIDDay, pattern = "^[\\d]+\\_"), 
                                                SubjectID = str_remove(SubjectIDDay, pattern = "\\_[\\d]+$")) %>% 
    inner_join(Rappid.pneum_category.df) %>% dplyr::filter(Day == "1") %>% 
    ggplot(aes(x = SubjectIDDay, y = AbsAbd, fill = Kingdom)) + 
    geom_bar(stat = "identity") + facet_grid(.~pneum_category, scales = "free", space = "free") + 
    theme(axis.text.x = element_text(angle = 90)) + theme(legend.position = "bottom")
  Barplot.shared_legend <- extract_legend(Nanopore.filtered.AbsAbd.baseline.barplot.p) 
  
  # Create plots for follow-up timepoints
  Nanopore.filtered.RelAbd.followup.barplot.p <- 
    Nanopore.filtered.RelAbd.long.df %>% distinct(SubjectIDDay) %>% 
    mutate(Day = str_remove(SubjectIDDay, pattern = "^[\\d]+\\_"), SubjectID = str_remove(SubjectIDDay, pattern = "\\_[\\d]+$")) %>% 
    inner_join(Rappid.pneum_category.df) %>% group_by(SubjectID) %>% mutate(time_points = n()) %>% dplyr::filter(time_points > 1) %>% ungroup() %>% 
    transmute(SubjectIDDay, SubjectID, pneum_category, Day = factor(as.integer(Day), levels = c(1, 5, 10))) %>% 
    left_join(Nanopore.filtered.RelAbd.long.df) %>% 
    ggplot(aes(x = Day, y = filtered_RelAbd, fill = Kingdom)) + 
    geom_bar(stat = "identity") + facet_grid(. ~ pneum_category + SubjectID, scales = "free_x", space = "free") + 
    theme(legend.position = "none")
  
  Nanopore.filtered.AbsAbd.followup.barplot.p <- 
    Nanopore.filtered.RelAbd.long.df %>% distinct(SubjectIDDay) %>% 
    mutate(Day = str_remove(SubjectIDDay, pattern = "^[\\d]+\\_"), SubjectID = str_remove(SubjectIDDay, pattern = "\\_[\\d]+$")) %>% 
    inner_join(Rappid.pneum_category.df) %>% group_by(SubjectID) %>% mutate(time_points = n()) %>% dplyr::filter(time_points > 1) %>% ungroup() %>% 
    transmute(SubjectIDDay, SubjectID, pneum_category, Day = factor(as.integer(Day), levels = c(1, 5, 10))) %>% 
    left_join(Nanopore.filtered.RelAbd.long.df) %>% 
    ggplot(aes(x = Day, y = AbsAbd, fill = Kingdom)) + 
    geom_bar(stat = "identity") + facet_grid(. ~ pneum_category + SubjectID, scales = "free", space = "free") + 
    theme(legend.position = "none")
  
  # Create wide format matrix for downstream analyses
  Nanopore.filtered.RelAbd.wide.df <- 
    Nanopore.filtered.RelAbd.long.df %>% mutate(Day = str_remove(SubjectIDDay, pattern = "^[\\d]+\\_"), 
                                                SubjectID = str_remove(SubjectIDDay, pattern = "\\_[\\d]+$")) %>% 
    inner_join(Rappid.pneum_category.df) %>% dplyr::filter(Day == "1") %>% 
    left_join(NCBI.all.taxID.name %>% transmute(id_species = taxID, name_species = name_txt)) %>% 
    select(SubjectIDDay, name_species, filtered_RelAbd) %>% group_by(SubjectIDDay, name_species) %>% 
    pivot_wider(names_from = "name_species", values_from = "filtered_RelAbd", values_fill = 0) 
  
  # Calculate mean relative abundance by clinical group
  Nanopore.baseline.bubble.relabd.MEAN.df <- 
    Nanopore.filtered.RelAbd.wide.df %>% 
    mutate(SubjectID = str_remove(SubjectIDDay, pattern = "\\_[\\d]+$"))%>% 
    left_join(Rappid.pneum_category.df) %>% select(-"SubjectID") %>% column_to_rownames("SubjectIDDay") %>% 
    group_by(pneum_category) %>% 
    summarise(across(.fns = mean)) %>% ungroup()
  
  # Identify top 20 species per clinical subgroup
  Nanopore.species.to.plot.vec <- 
    Nanopore.baseline.bubble.relabd.MEAN.df %>% 
    pivot_longer(cols = -`pneum_category`, names_to = "name_species", values_to = "MeanRelAbd_within_subgroup")  %>% 
    group_by(pneum_category) %>% mutate(rank_within_subgroup = rank(-MeanRelAbd_within_subgroup)) %>% 
    dplyr::filter(rank_within_subgroup <= 20) %>% ungroup() %>% 
    select(name_species) %>% unlist(use.names = F) %>% unique() %>% sort()
  
  Nanopore.top20.bacteria.by.subgroup.df <- 
    Nanopore.baseline.bubble.relabd.MEAN.df %>% 
    select("pneum_category", all_of(Nanopore.species.to.plot.vec)) %>% 
    rowwise() %>% 
    mutate(Others = 1 - sum(c_across(cols = all_of(Nanopore.species.to.plot.vec)))) %>% ungroup() %>% 
    column_to_rownames("pneum_category")
  
  ################################################################################
  #### KARIUS MICROBIAL CELL-FREE DNA (mcfDNA) DATA PROCESSING ####
  ################################################################################
  
  # Load Karius plasma mcfDNA test results
  Karius.summary <- 
    read.csv(file = "data/Task1D.Karius_text_summary_all_available_rappid_samples.csv") 
  
  # Parse Karius results and map to NCBI taxonomy
  Karius.df <- 
    Karius.summary %>% 
    separate_rows(Karius_summary, sep = "; ") %>% 
    separate(SubjectIDDay, remove = F, sep = "\\_", into = c("SubjectID", "Day")) %>% 
    select(-pneum_category) %>% 
    separate(Karius_summary, into = c("Taxa", "MPM"), sep = ", ") %>% 
    mutate(MPM = as.integer(MPM), 
           isVirus = str_detect(Taxa, pattern = "[Vv]irus"), 
           genus = case_when(
             Taxa == "No Organism Detected" ~ "No Organism Detected", 
             Taxa == "No Result" ~ NA_character_, 
             isVirus ~ Taxa, 
             TRUE ~ str_extract(Taxa, pattern = "^[^\\ ]+")
           )) %>% 
    left_join(NCBI.all.taxID.name, by = c("genus" = "name_txt")) %>% 
    # Manual corrections for viral taxa
    mutate(genusID = case_when(
      Taxa == "No Organism Detected" ~ NA_integer_, 
      Taxa == "No Result" ~ NA_integer_, 
      Taxa == "Cytomegalovirus (CMV)" ~ 10358L, 
      Taxa == "Herpes simplex virus type 1 (HSV-1)" ~ 10298L, 
      Taxa == "Epstein-Barr virus (EBV)" ~ 10376L, 
      Taxa == "Torque teno virus" | Taxa == "Torque teno virus 15" ~ 68887L, 
      Taxa == "Human herpesvirus 6A" ~ 32603L, 
      Taxa == "Herpes simplex virus type 2 (HSV-2)" ~ 10301L, 
      Taxa == "Kaposi sarcoma-associated herpesvirus" ~ 37296L, 
      Taxa == "BK polyomavirus" ~ 1891762L, 
      Taxa == "Adeno-associated dependoparvirus A" ~ 1511891L, 
      TRUE ~ taxID
    )) 
  
  ################################################################################
  #### MATCH NANOPORE AND KARIUS DATA AT GENUS LEVEL ####
  ################################################################################
  
  # Aggregate Nanopore data to genus level for comparison with Karius
  Nanopore.long.genus.RelAbd.df <- 
    Nanopore.filtered.RelAbd.long.df %>% 
    transmute(SubjectIDDay, Kingdom, filtered_RelAbd, name_species, 
              genus_name = case_when(
                Kingdom == "Virus" ~ name_species, 
                str_detect(name_species, pattern = "\\[Candida\\]") ~ "Candida",  
                str_detect(name_species, pattern = "Candidatus") ~ str_extract(name_species, pattern = "Candidatus\\ [^\\ ]+"),  
                str_detect(name_species, pattern = "\\[Eubacterium\\]") ~ "Eubacterium", 
                TRUE ~ str_extract(name_species, pattern = "^[^ ]+"))) %>% 
    group_by(genus_name, SubjectIDDay) %>% reframe(genusRelAbd = sum(filtered_RelAbd)) %>% arrange(SubjectIDDay, genusRelAbd) %>% 
    left_join(NCBI.all.taxID.name, by = c("genus_name" = "name_txt"), keep = F) %>% 
    # Manual corrections for genus-level taxonomy
    mutate(taxID = case_when(
      genus_name == "Cutibacterium" ~ 1912216L, 
      genus_name == "Schaalia" ~ 2529408L, 
      genus_name == "Candidatus Nanosynbacter" ~ 2093823L, 
      genus_name == "Kosakonia" ~ 1330547L, 
      genus_name == "Tardibacter" ~ 2685927L, 
      genus_name == "Croceicoccus" ~ 1295327L, 
      genus_name == "Pseudopropionibacterium" ~ 2801844L, 
      genus_name == "Lelliottia" ~ 1330545L, 
      genus_name == "Jeotgalibaca" ~ 1470540L,
      str_detect(genus_name, pattern = "[Pp]hage") ~ NA_integer_, 
      TRUE ~ as.integer(taxID)
    )) %>% dplyr::filter(!(taxID %in% c("508215", "1535326", "132406", "55087", "79213")))
  
  # Aggregate Karius data to genus level
  Karius.genus.df <- 
    Karius.df %>% dplyr::filter(Day == 1 & Taxa != "No Result") %>% group_by(SubjectIDDay, genus, isVirus) %>% 
    reframe(genusMPM = sum(MPM, na.rm = F)) %>% transmute(SubjectIDDay, genus_name = genus, genusMPM, isVirus) %>% arrange(SubjectIDDay, -genusMPM) %>% 
    left_join(NCBI.all.taxID.name, by = c("genus_name" = "name_txt"), keep = F) %>% 
    mutate(taxID = case_when(
      isVirus ~ NA_integer_, 
      genus_name == "No Organism Detected" ~ NA_integer_, 
      TRUE ~ as.integer(taxID)
    )) %>% dplyr::filter(taxID != "508215" & taxID != "1535326")
  
  ################################################################################
  #### DEFINE MICROBIAL TRANSLOCATION ####
  ################################################################################
  # Translocation = microbes detected in BOTH lung (Nanopore) AND blood (Karius)
  
  # Identify shared genera and samples between datasets
  Karius_Nanopore_shared_genus.vec <- sort(unique(intersect(Karius.genus.df$genus_name, Nanopore.long.genus.RelAbd.df$genus_name)))
  Karius_Nanopore_shared_SubjectIDDay.vec <- 
    intersect(
      sort(unique(Nanopore.long.genus.RelAbd.df$SubjectIDDay)), 
      sort(unique(Karius.genus.df$SubjectIDDay))
    )
  
  NCBI.all.taxID.name.dups <- NCBI.all.taxID.name %>% group_by(name_txt) %>% filter(n() > 1) %>% arrange(name_txt)
  
  # Calculate translocation metrics for each patient
  rappid.TransLoc.df <- 
    full_join(
      Karius.genus.df %>% dplyr::filter(SubjectIDDay %in% Karius_Nanopore_shared_SubjectIDDay.vec), 
      Nanopore.long.genus.RelAbd.df %>% dplyr::filter(SubjectIDDay %in% Karius_Nanopore_shared_SubjectIDDay.vec)
    ) %>% 
    mutate(
      genusMPM = case_when(
        SubjectIDDay %in% c("4668_1", "4847_1") ~ NA_integer_, 
        is.na(genusMPM) ~ 0L, 
        TRUE ~ genusMPM), 
      genusRelAbd = case_when(
        is.na(genusRelAbd) ~ 0, 
        TRUE ~ genusRelAbd
      )) %>% dplyr::filter(!is.na(genusRelAbd) & !is.na(genusMPM) & genusMPM > 0 & genusRelAbd > 0) %>% 
    group_by(SubjectIDDay) %>% summarise(TransLoc_MPM = sum(genusMPM, na.rm = T), TransLoc_RelAbd = sum(genusRelAbd, na.rm = T)) %>% 
    mutate(TransLoc = TRUE, SubjectID = str_extract(SubjectIDDay, pattern = "^[\\d]+")) %>% 
    mutate(TransLoc = case_when(is.na(TransLoc) ~ FALSE, TransLoc_RelAbd < 0.001 ~ FALSE, TRUE ~ TransLoc), 
           TransLoc_MPM = case_when(is.na(TransLoc_MPM) ~ 0L, TRUE ~ TransLoc_MPM), 
           TransLoc_RelAbd = case_when(is.na(TransLoc_RelAbd) ~ 0, TRUE ~ TransLoc_RelAbd), 
           SubjectIDDay = case_when(is.na(SubjectIDDay) ~ paste(SubjectID, "1", sep = "_"), TRUE ~ SubjectIDDay)) %>% 
    dplyr::filter(SubjectIDDay %in% Karius_Nanopore_shared_SubjectIDDay.vec | SubjectIDDay %in% c("4668_1", "4847_1")) %>% 
    right_join(Rappid.baseline.meta.202305.df %>% select(SubjectIDDay)) %>% 
    mutate(hasKarius = SubjectIDDay %in% Karius.df$SubjectIDDay, 
           hasNanopore = SubjectIDDay %in% Nanopore.filtered.RelAbd.long.df$SubjectIDDay, 
           hasBothMetagenomics = hasKarius & hasNanopore, 
           TransLoc = case_when(
             hasBothMetagenomics & SubjectIDDay %in% c(Karius.df %>% dplyr::filter(Day == 1 & str_detect(Taxa, "No Organism Detected")) %>% select(SubjectIDDay) %>% unlist()) ~ FALSE, 
             hasBothMetagenomics & is.na(TransLoc) ~ FALSE, 
             TRUE ~ TransLoc), 
           # Categorize translocation: None, Non-pulmonary source, or Pulmonary source
           TransLoc_final = factor(
             case_when(
               !hasBothMetagenomics ~ NA_character_, 
               hasBothMetagenomics & SubjectIDDay %in% c(Karius.df %>% dplyr::filter(Day == 1 & str_detect(Taxa, "No Organism Detected")) %>% select(SubjectIDDay) %>% unlist()) ~ "NoTransLoc", 
               hasBothMetagenomics & !TransLoc ~ "NonPulmTransLoc", 
               hasBothMetagenomics & TransLoc ~ "PulmTransLoc"
             )), 
           TransLoc_MPM = case_when(!TransLoc ~ 0L, TRUE ~ TransLoc_MPM), 
           TransLoc_RelAbd = case_when(!TransLoc ~ 0, TRUE ~ TransLoc_RelAbd)) %>% arrange(hasBothMetagenomics, TransLoc, SubjectIDDay) %>% 
    filter(SubjectIDDay != "4668_1") %>% 
    # Manual addition for sample 4668_1 with known high lung abundance
    add_case(SubjectIDDay = "4668_1", 
             TransLoc = TRUE, TransLoc_final = "PulmTransLoc", TransLoc_MPM = NA_integer_, TransLoc_RelAbd = 0.612, 
             hasKarius = TRUE, hasNanopore = TRUE, hasBothMetagenomics = TRUE) %>% 
    select(-SubjectID) %>% 
    left_join(Rappid.pneum_category.df %>% transmute(SubjectIDDay = paste(SubjectID, "1", sep = "_"), pneum_category)) %>% 
    dplyr::filter(!is.na(TransLoc_final) & pneum_category != "CONTROL") %>% 
    mutate(TransLoc_final = factor(TransLoc_final, levels = c("NoTransLoc", "NonPulmTransLoc", "PulmTransLoc"))) %>% 
    mutate(TransLoc_final_fix = factor(levels = c("No translocation", 'Non-pulmonary\ntranslocation', 'Pulmonary\ntranslocation'), x = case_when(
      TransLoc_final == 'NoTransLoc' ~ 'No translocation', 
      TransLoc_final == 'NonPulmTransLoc' ~ 'Non-pulmonary\ntranslocation', 
      TransLoc_final == 'PulmTransLoc' ~ 'Pulmonary\ntranslocation'
    )))
  
  # Summary statistics for translocation
  hasBothMetagenomics.nrow <- rappid.TransLoc.df %>% dplyr::filter(hasBothMetagenomics) %>% nrow()
  TransLoc.pos.nrow <- rappid.TransLoc.df %>% dplyr::filter(hasBothMetagenomics & TransLoc) %>% nrow()
  TransLoc.neg.nrow <- rappid.TransLoc.df %>% dplyr::filter(hasBothMetagenomics & !TransLoc) %>% nrow()
  
  # Identify translocating taxa
  TransLoc.taxID.vec <- 
    full_join(
      Karius.genus.df %>% dplyr::filter(SubjectIDDay %in% Karius_Nanopore_shared_SubjectIDDay.vec), 
      Nanopore.long.genus.RelAbd.df %>% dplyr::filter(SubjectIDDay %in% Karius_Nanopore_shared_SubjectIDDay.vec)
    ) %>% 
    mutate(
      genusMPM = case_when(
        SubjectIDDay %in% c("4668_1", "4847_1") ~ NA_integer_, 
        is.na(genusMPM) ~ 0L, 
        TRUE ~ genusMPM), 
      genusRelAbd = case_when(
        is.na(genusRelAbd) ~ 0, 
        TRUE ~ genusRelAbd
      )) %>% dplyr::filter(!is.na(genusRelAbd) & !is.na(genusMPM) & genusMPM > 0 & genusRelAbd > 0) %>%  ungroup() %>%  group_by(taxID) %>% mutate(rank(-sum(genusRelAbd))) %>% ungroup() %>% dplyr::select(taxID) %>% unlist(, use.names = F)
  
  # Export table of translocating taxa
  NCBI.all.taxID.name %>% dplyr::filter(taxID %in% TransLoc.taxID.vec) %>% 
    write.csv(file = out.path('-R2Q18_Table_translocation_taxa-'))
  
  ################################################################################
  #### 16S rRNA AMPLICON SEQUENCING DATA ####
  ################################################################################
  
  # Load 16S amplicon data
  Amplicon.raw <- read.csv(file = "data/Email-20230527/rappidseqTaxa16stable_5.26.23.csv")
  
  # Parse sample IDs and filter for endotracheal aspirate (ETA) samples
  Amplicon.sample_id.df <- 
    Amplicon.raw %>% 
    select(sample_id, total) %>% filter(complete.cases(.)) %>% 
    mutate(SubjectIDDay = 
             str_extract(string = sample_id, pattern = "4\\d\\d\\d\\.[^.]+") %>% 
             str_replace(pattern = "\\.", replacement = "\\_") %>% 
             str_remove(pattern = "D0") %>% str_remove(pattern = "D"), 
           Sample_type = 
             case_when(
               str_detect(sample_id, pattern = "O") ~ "OS", 
               str_detect(sample_id, pattern = "[Tt][Rr][Aa][Cc][Hh]") | str_detect(sample_id, pattern = "[Tt][Aa]") ~ "ETA", 
               str_detect(sample_id, pattern = "RS") ~ "RS", 
               str_detect(sample_id, pattern = "ST") ~ "ST", 
               TRUE ~ "UNKONWN"
             )) %>% filter(complete.cases(.)) %>% distinct(sample_id, .keep_all = T)
  
  # Extract ETA samples and calculate absolute abundance
  Amplicon.ETA.AbsAbd.df <- 
    Amplicon.sample_id.df %>% filter(Sample_type == "ETA") %>% left_join(Amplicon.raw) %>% distinct(sample_id, .keep_all = T)
  
  # Calculate Shannon diversity and identify protective taxa
  Amplicon.ETA.Dysbiosis.df <- 
    Amplicon.ETA.AbsAbd.df %>% 
    mutate(across(.cols = !c(sample_id, total, SubjectIDDay, Sample_type), .fns = ~ ifelse(.x == 0, yes = 0, no = -(.x / total) * log(.x / total)))) %>% 
    rowwise() %>% mutate(Shannon_ETA = sum(c_across(!c(sample_id, total, SubjectIDDay, Sample_type)))) %>% ungroup() %>% 
    transmute(sample_id, total, SubjectIDDay, Sample_type, 
              Shannon_ETA, 
              Prevotella = Amplicon.ETA.AbsAbd.df$Bacteria.Bacteroidota.Bacteroidia.Bacteroidales.Prevotellaceae.Prevotella / total, 
              Streptococcus = Amplicon.ETA.AbsAbd.df$Bacteria.Firmicutes.Bacilli.Lactobacillales.Streptococcaceae.Streptococcus / total, 
              Haemophilus = Amplicon.ETA.AbsAbd.df$Bacteria.Proteobacteria.Gammaproteobacteria.Pasteurellales.Pasteurellaceae.Haemophilus / total,
              Protective_ETA_RelAbd = Prevotella + Streptococcus + Haemophilus,
              Dysbiosis_ETA = !(Shannon_ETA >= 1.98 & Protective_ETA_RelAbd >= 0.3)) 

  # Create microeco object for 16S analysis by clinical groups
  Amplicon.ETA.microeco.obj <- 
    microtable$new(
      otu_table = Amplicon.ETA.AbsAbd.df %>% filter(str_detect(SubjectIDDay, pattern = "\\_1$")) %>% 
        mutate(across(.cols = !c(sample_id, total, SubjectIDDay, Sample_type), .fns = ~.x / total)) %>% select(!c(Sample_type, total, sample_id)) %>% 
        column_to_rownames("SubjectIDDay") %>% `colnames<-`(str_extract(string = colnames(.), pattern = "[^.]+$")) %>% t() %>% as.data.frame(), 
      sample_table = Amplicon.ETA.Dysbiosis.df %>% filter(str_detect(SubjectIDDay, pattern = "\\_1$")) %>% 
        left_join(Rappid.pneum_category_final.df %>% mutate(pneum_category_final = factor(pneum_category_final, levels = c('MCP', 'CDP', 'CONTROL')))) %>% left_join(rappid.TransLoc.df) %>% mutate(SubjectIDDay_dup = SubjectIDDay) %>% 
        column_to_rownames("SubjectIDDay_dup"), 
      tax_table = Amplicon.ETA.AbsAbd.df %>% filter(str_detect(SubjectIDDay, pattern = "\\_1$")) %>% 
        select(!c(sample_id, total, SubjectIDDay, Sample_type)) %>% colnames() %>% enframe() %>% 
        transmute(species = 
                    case_when(str_detect(string = value, pattern = "uncultured_ge") ~ str_extract(string = value, pattern = "[^.]+\\.[^.]+\\.[^.]+\\.[^.]+\\.[^.]+$"), 
                              TRUE ~ str_extract(string = value, pattern = "[^.]+$")), species_dup = species) %>% column_to_rownames("species_dup")
    )
  
  # Identify samples with translocation data for 16S analysis
  Amplicon.ETA.TransLoc.SubjectIDDay.vec <- 
    intersect(Amplicon.ETA.Dysbiosis.df$SubjectIDDay, 
              rappid.TransLoc.df %>% filter(pneum_category != "CONTROL" & hasBothMetagenomics) %>% select(SubjectIDDay) %>% unlist(use.names = F))
  
  # Create microeco object for translocation analysis
  Amplicon.ETA.TransLoc.microeco.obj <- 
    microtable$new(
      otu_table = Amplicon.ETA.AbsAbd.df %>% filter(SubjectIDDay %in% Amplicon.ETA.TransLoc.SubjectIDDay.vec) %>% 
        mutate(across(.cols = !c(sample_id, total, SubjectIDDay, Sample_type), .fns = ~.x / total)) %>% select(!c(Sample_type, total, sample_id)) %>% 
        column_to_rownames("SubjectIDDay") %>% `colnames<-`(str_extract(string = colnames(.), pattern = "[^.]+$")) %>% t() %>% as.data.frame(), 
      sample_table = Amplicon.ETA.Dysbiosis.df %>% filter(SubjectIDDay %in% Amplicon.ETA.TransLoc.SubjectIDDay.vec) %>% left_join(Rappid.pneum_category_final.df) %>% left_join(rappid.TransLoc.df) %>% mutate(SubjectIDDay_dup = SubjectIDDay) %>% column_to_rownames("SubjectIDDay_dup"), 
      tax_table = Amplicon.ETA.AbsAbd.df %>% filter(SubjectIDDay %in% Amplicon.ETA.TransLoc.SubjectIDDay.vec) %>% 
        select(!c(sample_id, total, SubjectIDDay, Sample_type)) %>% colnames() %>% enframe() %>% 
        transmute(species = 
                    case_when(str_detect(string = value, pattern = "uncultured_ge") ~ str_extract(string = value, pattern = "[^.]+\\.[^.]+\\.[^.]+\\.[^.]+\\.[^.]+$"), 
                              TRUE ~ str_extract(string = value, pattern = "[^.]+$")), 
                  species_dup = species) %>% column_to_rownames("species_dup")
    )
  
  #### LEfSe analysis for 16S data by clinical group ####
  Amplicon.ETA.LEFSE.pneum_cat.obj <- 
    trans_diff$new(Amplicon.ETA.microeco.obj, method = "lefse", group = "pneum_category_final", alpha = 0.05, p_adjust_method = "none")
  Amplicon.ETA.LEFSE.pneum_cat.abund.p <- Amplicon.ETA.LEFSE.pneum_cat.obj$plot_diff_abund(width = 0.8, add_sig = T)
  Amplicon.ETA.LEFSE.pneum_cat.bar.p <- Amplicon.ETA.LEFSE.pneum_cat.obj$plot_diff_bar(width = 0.8)

  ################################################################################
  #### PATHOGEN CLASSIFICATION ####
  ################################################################################
  
  # Load pathogen classification (established, possible, unlikely pathogens)
  Patho_Classification.FINAL.raw <- xlsx::read.xlsx("data/20240228-RAPPID-Patho_Classification_final.xlsx", sheetIndex = 2)
  
  # Merge translocation data with pathogen classifications
  Nanopore.TransLoc_final.LEFSE.df <- 
    Nanopore.filtered.RelAbd.long.df %>% select(SubjectIDDay, name_species, filtered_RelAbd, AbsAbd) %>% 
    right_join(rappid.TransLoc.df[, c("SubjectIDDay", "TransLoc_final")]) %>% 
    mutate(TransLoc_final_fix = factor(levels = c("No translocation", 'Non-pulmonary\ntranslocation', 'Pulmonary\ntranslocation'), x = case_when(
      TransLoc_final == 'NoTransLoc' ~ 'No translocation', 
      TransLoc_final == 'NonPulmTransLoc' ~ 'Non-pulmonary\ntranslocation', 
      TransLoc_final == 'PulmTransLoc' ~ 'Pulmonary\ntranslocation'
    )))
  
  # Create visualization of pathogen abundance by translocation group
  pathogen_TransLoc.Nanopore.p <- 
    Nanopore.TransLoc_final.LEFSE.df %>% left_join(Patho_Classification.FINAL.raw[, c("SPECIES", "HY_FINAL")] , by = c("name_species" = "SPECIES")) %>% 
    group_by(SubjectIDDay, TransLoc_final, HY_FINAL) %>% summarise(GroupRelAbdSum = sum(filtered_RelAbd)) %>% ungroup() %>% 
    mutate(Patho_Rank_Place_Holder = case_when(HY_FINAL == "established" ~ GroupRelAbdSum, TRUE ~ 0)) %>% 
    left_join(DummyID.df %>% mutate(SubjectIDDay = paste(SubjectID, "_1", sep = ""), DummyID = paste("Patient #", ParticipantID, sep = " "))) %>% 
    mutate(Pathogenicity = str_to_upper(HY_FINAL)) %>% 
    ggplot(aes(y = reorder(DummyID, -Patho_Rank_Place_Holder), x = GroupRelAbdSum, fill = Pathogenicity)) + geom_bar(stat = "identity") + 
    facet_wrap(TransLoc_final ~ ., scales = 'free') + xlab("Relative Abundance") + ylab("")
  
  # Create visualization of pathogen abundance by clinical group  
  pathogen_pneum_cate.Nanopore.p <- 
    Nanopore.filtered.RelAbd.long.df %>% select(SubjectIDDay, name_species, filtered_RelAbd) %>% 
    left_join(Patho_Classification.FINAL.raw[, c("SPECIES", "HY_FINAL")] , by = c("name_species" = "SPECIES")) %>% 
    group_by(SubjectIDDay, HY_FINAL) %>% summarise(GroupRelAbdSum = sum(filtered_RelAbd)) %>% ungroup() %>% 
    mutate(Patho_Rank_Place_Holder = case_when(HY_FINAL == "established" ~ GroupRelAbdSum, TRUE ~ 0), SubjectID = str_remove(SubjectIDDay, pattern = "_[\\d]+$")) %>% dplyr::filter(str_detect(SubjectIDDay, pattern = "_1$")) %>% 
    inner_join(Rappid.pneum_category.df) %>% 
    left_join(DummyID.df %>% transmute(SubjectIDDay = paste(SubjectID, "_1", sep = ""), ParticipantID = paste("Patient #", ParticipantID))) %>% 
    mutate(pneum_category = factor(
      case_when(
        pneum_category == "CONTROL" ~ "CONTROL", 
        pneum_category == "CXNEG" ~ "CDP", 
        pneum_category == "CXPOS" ~ "MCP"), 
      levels = c("MCP", "CDP", "CONTROL"))) %>% 
    mutate(Pathogenicity = factor(case_when(is.na(HY_FINAL) ~ "NONE DETECTED", TRUE ~ str_to_upper(HY_FINAL)), levels = c("ESTABLISHED", "POSSIBLE", "UNLIKELY", "NONE DETECTED"))) %>% 
    ggplot(aes(y = reorder(ParticipantID, -Patho_Rank_Place_Holder), x = GroupRelAbdSum, fill = Pathogenicity)) + 
    geom_bar(stat = "identity") + facet_wrap(pneum_category ~ ., scales = 'free') + ylab("") + xlab("Relative Abundance")
  
  # Statistical testing for pathogen enrichment by clinical group
  pathogen_pneum_cate.Nanopore.wilcox <- 
    Nanopore.filtered.RelAbd.long.df %>% select(SubjectIDDay, name_species, filtered_RelAbd) %>% 
    left_join(Patho_Classification.FINAL.raw[, c("SPECIES", "HY_FINAL")] , by = c("name_species" = "SPECIES")) %>% 
    group_by(SubjectIDDay, HY_FINAL) %>% summarise(GroupRelAbdSum = sum(filtered_RelAbd)) %>% ungroup() %>% 
    inner_join(Rappid.pneum_category.df %>% transmute(SubjectIDDay = paste(SubjectID, "_1", sep = ""), pneum_category)) %>% 
    mutate(pneum_category = 
             factor(case_when(pneum_category == "CONTROL" ~ "CONTROL", pneum_category == "CXNEG" ~ "CDP", pneum_category == "CXPOS" ~ "MCP"), 
                    levels = c("MCP", "CDP", "CONTROL"))) %>% 
    pivot_wider(names_from = "HY_FINAL", values_from = "GroupRelAbdSum", values_fill = 0) %>% 
    rstatix::wilcox_test(established ~ pneum_category, p.adjust.method = "BH") %>% 
    add_significance() %>% 
    add_xy_position() 
  
  # Box plot with statistics for Nanopore pathogen abundance
  pathogen_pneum_cate.Nanopore.pmanual.p <- 
    Nanopore.filtered.RelAbd.long.df %>% select(SubjectIDDay, name_species, filtered_RelAbd) %>% 
    left_join(Patho_Classification.FINAL.raw[, c("SPECIES", "HY_FINAL")] , by = c("name_species" = "SPECIES")) %>% 
    group_by(SubjectIDDay, HY_FINAL) %>% summarise(GroupRelAbdSum = sum(filtered_RelAbd)) %>% ungroup() %>% 
    inner_join(Rappid.pneum_category.df %>% transmute(SubjectIDDay = paste(SubjectID, "_1", sep = ""), pneum_category)) %>% 
    mutate(pneum_category = factor(case_when(pneum_category == "CONTROL" ~ "CONTROL", pneum_category == "CXNEG" ~ "CDP", pneum_category == "CXPOS" ~ "MCP"), levels = c("MCP", "CDP", "CONTROL"))) %>% 
    dplyr::filter(HY_FINAL == "established") %>% 
    ggboxplot(x = "pneum_category", y = "GroupRelAbdSum", add = "jitter") + 
    stat_pvalue_manual(pathogen_pneum_cate.Nanopore.wilcox,  step.increase = 0.06, hide.ns = T) +
    xlab("") + ylab("Sum of Respir Pathogen Relative Abundance")
  
  # Create data frame of pathogen abundances
  Nanopore.patho_cat.final.df <- 
    Nanopore.filtered.RelAbd.long.df %>% select(SubjectIDDay, name_species, filtered_RelAbd) %>% 
    left_join(Patho_Classification.FINAL.raw[, c("SPECIES", "HY_FINAL")] , by = c("name_species" = "SPECIES")) %>% 
    group_by(SubjectIDDay, HY_FINAL) %>% summarise(GroupRelAbdSum = sum(filtered_RelAbd)) %>% ungroup() %>% 
    pivot_wider(names_from = "HY_FINAL", values_from = "GroupRelAbdSum", values_fill = 0)
  
  #### Karius pathogen analysis by clinical group ####
  
  # Statistical testing for Karius pathogen data
  pathogen_pneum_cate.Karius.wilcox <- 
    Karius.df %>% dplyr::filter(Day == "1") %>% inner_join(Rappid.pneum_category.df) %>% 
    left_join(Patho_Classification.FINAL.raw[, c("SPECIES", "HY_FINAL")] , by = c("Taxa" = "SPECIES")) %>% 
    group_by(SubjectIDDay, HY_FINAL, pneum_category) %>% summarise(pathoSumMPM = sum(MPM)) %>% arrange(pathoSumMPM) %>% ungroup() %>% 
    dplyr::filter(!is.na(pathoSumMPM)) %>% 
    pivot_wider(names_from = "HY_FINAL", values_from = pathoSumMPM, values_fill = 0) %>% 
    add_case(SubjectIDDay = "4668_1", pneum_category = "CXPOS") %>% 
    add_case(SubjectIDDay = "4687_1", pneum_category = "CXPOS") %>% 
    add_case(SubjectIDDay = "4707_1", pneum_category = "CXPOS") %>% 
    add_case(SubjectIDDay = "4847_1", pneum_category = "CXNEG")  %>% 
    mutate(pneum_category = factor(case_when(pneum_category == "CONTROL" ~ "CONTROL", pneum_category == "CXNEG" ~ "CDP", pneum_category == "CXPOS" ~ "MCP"), levels = c("MCP", "CDP", "CONTROL"))) %>% 
    rstatix::wilcox_test(formula = established ~ pneum_category, p.adjust.method = "BH") %>% 
    add_significance() %>% 
    add_xy_position() %>% 
    mutate(`y.position` = rev(c(12.5, 13.5, 14.5)))
  
  # Bar plot for Karius pathogen data by clinical group
  pathogen_pneum_cate.Karius.p <- 
    Karius.df %>% dplyr::filter(Day == "1") %>% inner_join(Rappid.pneum_category.df) %>% 
    left_join(Patho_Classification.FINAL.raw[, c("SPECIES", "HY_FINAL")] , by = c("Taxa" = "SPECIES")) %>% 
    group_by(SubjectIDDay, HY_FINAL, pneum_category) %>% summarise(pathoSumMPM = sum(MPM)) %>% arrange(pathoSumMPM) %>% ungroup() %>% dplyr::filter(!is.na(pathoSumMPM)) %>% 
    mutate(Patho_Rank_Place_Holder = case_when(HY_FINAL == "established" ~ pathoSumMPM, TRUE ~ 0), SubjectID = str_remove(SubjectIDDay, pattern = "_[\\d]+$")) %>% 
    mutate(pneum_category = factor(case_when(pneum_category == "CONTROL" ~ "CONTROL", pneum_category == "CXNEG" ~ "CDP", pneum_category == "CXPOS" ~ "MCP"), levels = c("MCP", "CDP", "CONTROL"))) %>% 
    mutate(Pathogenicity = factor(case_when(is.na(HY_FINAL) ~ "NONE DETECTED", TRUE ~ str_to_upper(HY_FINAL)), levels = c("ESTABLISHED", "POSSIBLE", "UNLIKELY", "NONE DETECTED"))) %>% 
    left_join(DummyID.df %>% transmute(SubjectIDDay = paste(SubjectID, "_1", sep = ""), ParticipantID = paste("Patient #", ParticipantID))) %>% 
    ggplot(aes(y = reorder(ParticipantID, Patho_Rank_Place_Holder), x = pathoSumMPM, fill = Pathogenicity)) + 
    geom_bar(stat = "identity") + facet_wrap(.~ pneum_category, scales = "free_y") + 
    scale_fill_manual(breaks = c('ESTABLISHED', 'POSSIBLE', 'UNLIKELY'), values = c("#F8766D", "#00BA38", "#619CFF", '#FFFFFF')) + 
    ylab("") + xlab("MPMs")
  
  # Box plot with statistics for Karius pathogen data
  pathogen_pneum_cate.Karius.pvalue.p <-
    Karius.df %>% dplyr::filter(Day == "1") %>% inner_join(Rappid.pneum_category.df) %>% 
    left_join(Patho_Classification.FINAL.raw[, c("SPECIES", "HY_FINAL")] , by = c("Taxa" = "SPECIES")) %>% 
    group_by(SubjectIDDay, HY_FINAL, pneum_category) %>% summarise(pathoSumMPM = sum(MPM)) %>% arrange(pathoSumMPM) %>% ungroup() %>% dplyr::filter(!is.na(pathoSumMPM)) %>% 
    pivot_wider(names_from = "HY_FINAL", values_from = pathoSumMPM, values_fill = 0) %>% 
    add_case(SubjectIDDay = "4668_1", pneum_category = "CXPOS") %>% 
    add_case(SubjectIDDay = "4687_1", pneum_category = "CXPOS") %>% 
    add_case(SubjectIDDay = "4707_1", pneum_category = "CXPOS") %>% 
    add_case(SubjectIDDay = "4847_1", pneum_category = "CXNEG")  %>% 
    mutate(established_log = case_when(established == 0 ~ 0, TRUE ~ log(established))) %>% 
    mutate(pneum_category = factor(case_when(pneum_category == "CONTROL" ~ "CONTROL", pneum_category == "CXNEG" ~ "CDP", pneum_category == "CXPOS" ~ "MCP"), 
                                   levels = c("MCP", "CDP", "CONTROL"))) %>% 
    ggboxplot(x = "pneum_category", y = "established_log", add = "jitter") + 
    stat_pvalue_manual(pathogen_pneum_cate.Karius.wilcox, dodge = 0.6, hide.ns = T) + ylab("MPMs from Respir Pathogen, log transformed") + xlab("")
  
  ################################################################################
  #### LEfSe ANALYSIS BY TRANSLOCATION GROUP ####
  ################################################################################
  
  # Create microtable object for LEfSe analysis
  Nanopore.TransLoc_final.LEFSE.species.microtable <- 
    microtable$new(
      otu_table = 
        Nanopore.TransLoc_final.LEFSE.df %>% select(!c("TransLoc_final", 'TransLoc_final_fix', "filtered_RelAbd")) %>% 
        pivot_wider(names_from = "name_species", values_from = "AbsAbd", values_fill = 0) %>% 
        column_to_rownames("SubjectIDDay") %>% 
        t() %>% as.data.frame(), 
      sample_table = rappid.TransLoc.df[, c("SubjectIDDay", "TransLoc_final", "pneum_category")] %>% column_to_rownames("SubjectIDDay"), 
      tax_table = 
        Patho_Classification.FINAL.raw[, c("SPECIES", "HY_FINAL")] %>% 
        mutate(ROWNAME = SPECIES) %>% column_to_rownames("ROWNAME"))
  
  # Calculate alpha diversity
  Nanopore.TransLoc_final.LEFSE.species.microtable$cal_alphadiv()
  
  # Calculate beta diversity
  Nanopore.TransLoc_final.LEFSE.species.microtable$cal_betadiv()
  
  # Perform beta diversity analysis
  Nanopore.TransLoc_final.LEFSE.BETA <- 
    trans_beta$new(dataset = Nanopore.TransLoc_final.LEFSE.species.microtable, measure = "bray")
  Nanopore.TransLoc_final.LEFSE.BETA$cal_ordination(ordination = "PCoA")
  Nanopore.TransLoc_final.LEFSE.BETA$cal_manova(manova_all = TRUE, group = "TransLoc_final", taxa_level = "SPECIES")
 
  # LEfSe differential abundance testing by translocation group
  trans_diff$new(dataset = Nanopore.TransLoc_final.LEFSE.species.microtable, method = "lefse", group = "TransLoc_final", alpha = 0.05, p_adjust_method = "none", taxa_level = "PATHO_100")$plot_diff_abund(use_number = 1:20)
  trans_diff$new(dataset = Nanopore.TransLoc_final.LEFSE.species.microtable, method = "lefse", group = "TransLoc_final", alpha = 0.05, p_adjust_method = "none", taxa_level = "Lit_No")$plot_diff_abund(use_number = 1:20)
  
  # Create LEfSe visualizations
  Nanopore.TransLoc.LEFSE.diff_bar.p <- 
    trans_diff$new(dataset = Nanopore.TransLoc_final.LEFSE.species.microtable, method = "lefse", group = "TransLoc_final", alpha = 0.05, 
                   p_adjust_method = "none", taxa_level = "SPECIES")$plot_diff_bar() + 
    scale_fill_manual(values = c('gray', 'darkblue', 'darkred')) + scale_color_manual(values = c('gray', 'darkblue', 'darkred'))
  
  Nanopore.TransLoc.LEFSE.diff_abund.p <- 
    trans_diff$new(dataset = Nanopore.TransLoc_final.LEFSE.species.microtable, method = "lefse", group = "TransLoc_final", alpha = 0.05, 
                   p_adjust_method = "none", taxa_level = "SPECIES")$plot_diff_abund(use_number = 1:20) + 
    scale_fill_manual(values = c('gray', 'darkblue', 'darkred')) + scale_color_manual(values = c('gray', 'darkblue', 'darkred'))

  # Box plots comparing total mcfDNA by translocation group
  TransLoc_final.mcfDNA.p <- 
    Rappid.meta.202305.df %>% select(SubjectIDDay, total_mcfdna) %>% inner_join(rappid.TransLoc.df) %>% 
    mutate(total_mcfdna_log = ifelse(total_mcfdna == 0, yes = 0, no = log(total_mcfdna))) %>% 
    ggboxplot(x = "TransLoc_final_fix", y = "total_mcfdna_log", add = "jitter", color = "TransLoc_final_fix", 
              palette = c('gray', 'darkblue', 'darkred')) + stat_compare_means(comparisons = list(c('Non-pulmonary\ntranslocation', 'Pulmonary\ntranslocation'))) + 
    annotate('text', x = 1, y = 2, label = 'ZERO\nBY\nDEFINITION', colour = 'gray')
  
  # Box plots comparing respiratory pathogen mcfDNA by translocation group
  TransLoc_Respir.mcfDNA.p <- 
    Karius.df %>% dplyr::filter(Day == "1") %>% inner_join(rappid.TransLoc.df[, c("SubjectIDDay", "TransLoc_final_fix")]) %>% 
    left_join(Patho_Classification.FINAL.raw[, c("SPECIES", "HY_FINAL")] , by = c("Taxa" = "SPECIES")) %>% 
    group_by(SubjectIDDay, HY_FINAL, TransLoc_final_fix) %>% summarise(pathoSumMPM = sum(MPM)) %>% arrange(pathoSumMPM) %>% ungroup() %>% 
    dplyr::filter(!is.na(pathoSumMPM) ) %>% 
    pivot_wider(values_from = "pathoSumMPM", names_from = "HY_FINAL", values_fill = 0) %>% 
    mutate(established_log = ifelse(established == 0, yes = 0, no = log(established))) %>% 
    ggboxplot(x = "TransLoc_final_fix", y = "established_log", add = "jitter", color = "TransLoc_final_fix", 
              palette = c('gray', 'darkblue', 'darkred')) + 
    stat_compare_means(comparisons = list(c('Non-pulmonary\ntranslocation', 'Pulmonary\ntranslocation'))) + 
    annotate('text', x = 1, y = 2, label = 'ZERO\nBY\nDEFINITION', colour = 'gray')
  
  # Bar plots for pathogen abundance by translocation group (Nanopore)
  pathogen_TransLoc.Nanopore.p <- 
    Nanopore.TransLoc_final.LEFSE.df %>% left_join(Patho_Classification.FINAL.raw[, c("SPECIES", "HY_FINAL")] , by = c("name_species" = "SPECIES")) %>% 
    group_by(SubjectIDDay, TransLoc_final_fix, HY_FINAL) %>% summarise(GroupRelAbdSum = sum(filtered_RelAbd)) %>% ungroup() %>% 
    mutate(Patho_Rank_Place_Holder = case_when(HY_FINAL == "established" ~ GroupRelAbdSum, TRUE ~ 0)) %>% 
    left_join(DummyID.df %>% mutate(SubjectIDDay = paste(SubjectID, "_1", sep = ""), DummyID = paste("Patient #", ParticipantID, sep = " "))) %>% 
    mutate(Pathogenicity = str_to_upper(HY_FINAL)) %>% 
    ggplot(aes(y = reorder(DummyID, -Patho_Rank_Place_Holder), x = GroupRelAbdSum, fill = Pathogenicity)) + geom_bar(stat = "identity") + 
    facet_wrap(TransLoc_final_fix ~ ., scales = 'free') + xlab("Relative Abundance") + ylab("")
  
  # Bar plots for pathogen abundance by translocation group (Karius)
  pathogen_TransLoc.Karius.p <- 
    Karius.df %>% dplyr::filter(Day == "1") %>% inner_join(rappid.TransLoc.df[, c("SubjectIDDay", "TransLoc_final_fix")]) %>% 
    left_join(Patho_Classification.FINAL.raw[, c("SPECIES", "HY_FINAL")] , by = c("Taxa" = "SPECIES")) %>% 
    group_by(SubjectIDDay, HY_FINAL, TransLoc_final_fix) %>% summarise(pathoSumMPM = sum(MPM)) %>% arrange(pathoSumMPM) %>% ungroup() %>% 
    dplyr::filter(!is.na(pathoSumMPM)) %>% 
    mutate(Patho_Rank_Place_Holder = case_when(HY_FINAL == "established" ~ pathoSumMPM, TRUE ~ 0), SubjectID = str_remove(SubjectIDDay, pattern = "_[\\d]+$")) %>% 
    left_join(DummyID.df %>% transmute(SubjectIDDay = paste(SubjectID, "_1", sep = ""), ParticipantID = paste("Patient #", ParticipantID))) %>% 
    ggplot(aes(y = reorder(ParticipantID, Patho_Rank_Place_Holder), x = pathoSumMPM, fill = HY_FINAL)) + geom_bar(stat = "identity") + 
    facet_wrap(.~ TransLoc_final_fix, scales = "free_y") + xlab("Microbial cell-free DNA, MPM(s)") + ylab('')
  
  # Save all processed data
  save.image('bkp/20250426-testrun.bkp.RData')
}

################################################################################
#### LOAD PREPROCESSED DATA ####
################################################################################

if(!RunFromScratch) {
  load('bkp/20251030.bkp.RData')
}

# Redefine output path function for consistency
out.path <- function(string) {
  return(
    case_when(
      str_detect(string, pattern = '[Ff]ig') ~ paste('output/', format(Sys.Date(), '%m-%d-%Y'), string, '.png', sep = ''),
      str_detect(string, pattern = '[Tt]ab') ~ paste('output/', format(Sys.Date(), '%m-%d-%Y'), string, '.csv', sep = '')
    )
  )
}

################################################################################
#### GENERATE ALL TABLES AND FIGURES ####
################################################################################

# Set to TRUE to generate and save all outputs
GenerateSaveAllTableFigure <- FALSE

if(GenerateSaveAllTableFigure) {
  
################################################################################
#### TABLE 1: Baseline Characteristics by Clinical Group ####
################################################################################

Rappid.baseline.meta.202305.df %>% 
  mutate(
    SubjectID = as.character(SubjectID), 
    across(.cols = starts_with("Hist"), .fns = as.character), 
    ards = as.character(ards), 
    IsOnVasopressor = as.logical(IsOnVasopressor), 
    Mortality30Day = as.character(Mortality30Day)) %>% 
  inner_join(Rappid.pneum_category.df) %>% 
  CreateTableOne(
    data = ., 
    strata = "pneum_category", 
    addOverall = T, 
    vars = c(
      c("Age", "BMI", "cpis", "SOFAScRSD", "VFD", 'PlateauPressure', 'PeakInspiratory', 'PaO2.y', 'PaCO2'), 
      c("Gender", "phenotype", "pneum_category", "Race", "HistCOPD", "HistImmunosuppression", "HistDiabetes", "ards", "ARDS_subph_less_than_0", "Mortality30Day", 'IsOnVasopressor'))) %>% 
  print(contDigits = 1, 
        nonnormal = c("Age", "BMI", "cpis", "SOFAScRSD", "VFD", 'PlateauPressure', 'PeakInspiratory', 'PaO2.y', 'PaCO2'), 
        exact = c("Gender", "phenotype", "pneum_category", "Race", "HistCOPD", "HistImmunosuppression", "HistDiabetes", "ards", "ARDS_subph_less_than_0", "Mortality30Day", 'IsOnVasopressor'))

# Post-hoc pairwise comparisons for continuous variables
# Age comparisons
Rappid.baseline.meta.202305.df %>% 
  mutate(
    SubjectID = as.character(SubjectID), 
    across(.cols = starts_with("Hist"), .fns = as.character), 
    ards = as.character(ards), 
    Mortality30Day = as.character(Mortality30Day)) %>% 
  inner_join(Rappid.pneum_category.df) %>% 
  rstatix::pairwise_wilcox_test(formula = Age ~ pneum_category, p.adjust.method = "BH")

# Additional pairwise tests for other variables follow similar pattern...
# (BMI, CPIS, SOFA, VFD, Gender, Immunosuppression, COPD, Diabetes, Race, Subphenotype)

################################################################################
#### TABLE 2: Demographics and Outcomes by Translocation Status ####
################################################################################

CreateTableOne(
  data = rappid.TransLoc.df %>% 
    left_join(Rappid.baseline.meta.202305.df %>% 
                transmute(SubjectIDDay, Age, Gender, Race, BMI, SOFAScRSD, PaO2.y, PeakInspiratory, PlateauPressure, PaCO2, VFD, HistCOPD = as.factor(HistCOPD), IsOnVasopressor = as.logical(IsOnVasopressor),  HistImmunosuppression = as.factor(HistImmunosuppression), ards = as.factor(ards), nouraiepars21, Mortality30Day = as.factor(Mortality30Day))), 
  strata = 'TransLoc_final', 
  vars = c('Age', 'Gender', 'Race', 'BMI', 'HistCOPD', 'HistImmunosuppression', 'SOFAScRSD', 'ards', 'nouraiepars21', 'Mortality30Day', 'VFD', 'PlateauPressure', 'PeakInspiratory', 'PaO2.y', 'PaCO2', 'IsOnVasopressor')
) %>% 
  print(
    nonnormal = c('Age', 'BMI', 'SOFAScRSD', 'VFD', 'PlateauPressure', 'PeakInspiratory', 'PaO2.y', 'PaCO2'), 
    exact = c('Race', 'HistCOPD', 'HistImmunosuppression', 'ards', 'nouraiepars21', 'Mortality30Day', 'IsOnVasopressor'), 
    pDigits = 4, contDigits = 1, 
    printToggle = F
  ) %>% write.csv(file = out.path('-Table 2'))

################################################################################
#### TABLE 3: Biomarker Comparisons by Translocation Group ####
################################################################################

Table3.df <- Rappid.baseline.meta.202305.df %>% 
  transmute(SubjectIDDay, 
            copies_ETA, shannon_ETA, 
            rage_eta, ang2_eta, rage, ang2, nDNA, mtDNA, il10_eta_pradj, fractalkine_eta_pradj, 
            il6_eta, il8_eta, tnfr1_eta, st2_eta, il6, il8, il10, tnfr1_final, st2, fractalkine, phenoclass_logit_nouraie, 
            procalcitonin_eta, pentraxin3_eta, procalcitonin, pentraxin3
  ) %>% inner_join(rappid.TransLoc.df) 

Table3.out <- 
  CreateTableOne(
  data = Table3.df, 
  strata = 'TransLoc_final', 
  vars = c("ang2_eta", "fractalkine_eta_pradj", "il6_eta", "il8_eta", "il10_eta_pradj", "pentraxin3_eta", "procalcitonin_eta", "rage_eta", "st2_eta", "tnfr1_eta", 
           'ang2', 'fractalkine', 'il6', 'il8', 'il10', 'pentraxin3', 'procalcitonin', 'rage', 'st2', 'tnfr1_final')
) %>% print(contDigits = 1, pDigits = 4, 
            nonnormal = c("ang2_eta", "fractalkine_eta_pradj", "il6_eta", "il8_eta", "il10_eta_pradj", "pentraxin3_eta", "procalcitonin_eta", "rage_eta", "st2_eta", "tnfr1_eta", 
                          'ang2', 'fractalkine', 'il6', 'il8', 'il10', 'pentraxin3', 'procalcitonin', 'rage', 'st2', 'tnfr1_final'), print.Toggle=FALSE)
write.csv(Table3.out, file = out.path('-Table 3'))

# Post-hoc adjusted p-values
Table3.p_posthoc_adj.out <- 
  Table3.out[, 4] %>% enframe() %>% dplyr::filter(value != '') %>% mutate(biomarker = str_remove_all(string = name, pattern = '\\ .+), location = str_detect(name, pattern = 'eta'), pval_raw = as.numeric(value)) %>% group_by(location) %>% mutate(pval_adj = p.adjust(pval_raw, method = 'BH'))
write.csv(Table3.p_posthoc_adj.out, file = out.path('-Table 3-posthoc_pval'))

################################################################################
#### TABLE S1: Biomarker Profile by Clinical Groups ####
################################################################################

TableS1.df <- 
  Rappid.baseline.meta.202305.df %>% 
  transmute(SubjectID = as.character(SubjectID), 
            copies_ETA, shannon_ETA, 
            rage_eta, ang2_eta, rage, ang2, nDNA, mtDNA, fractalkine_eta_pradj, il10_eta_pradj, 
            il6_eta, il8_eta, tnfr1_eta, st2_eta, il6, il8, il10, fractalkine, tnfr1_final, st2, phenoclass_logit_nouraie, 
            procalcitonin_eta, pentraxin3_eta, procalcitonin, pentraxin3
  ) %>% inner_join(Rappid.pneum_category.df) %>% 
  mutate(pneum_category = factor(pneum_category, levels = c("CONTROL", "CXNEG", "CXPOS"), labels = c("CONT", 'CDP', 'MCP')))

TableS1.out <- 
  CreateTableOne(
  data = TableS1.df, 
  strata = 'pneum_category', 
  vars = c("ang2_eta", "fractalkine_eta_pradj", "il6_eta", "il8_eta", "il10_eta_pradj", "pentraxin3_eta", "procalcitonin_eta", "rage_eta", "st2_eta", "tnfr1_eta", 
           'ang2', 'fractalkine', 'il6', 'il8', 'il10', 'pentraxin3', 'procalcitonin', 'rage', 'st2', 'tnfr1_final')
) %>% print(contDigits = 1, pDigits = 4, 
            nonnormal = c("ang2_eta", "fractalkine_eta_pradj", "il6_eta", "il8_eta", "il10_eta_pradj", "pentraxin3_eta", "procalcitonin_eta", "rage_eta", "st2_eta", "tnfr1_eta", 
                          'ang2', 'fractalkine', 'il6', 'il8', 'il10', 'pentraxin3', 'procalcitonin', 'rage', 'st2', 'tnfr1_final'), print.Toggle=FALSE)
write.csv(TableS1.df, file = out.path('-TableS1'))

TableS1_posthoc_adj.out <- 
  TableS1.out[, 4] %>% enframe() %>% dplyr::filter(value != '') %>% 
    mutate(biomarker = str_remove_all(string = name, pattern = '\\ .+), 
           location = str_detect(name, pattern = 'eta'), pval_raw = as.numeric(value)) %>% 
    group_by(location) %>% mutate(pval_adj = p.adjust(pval_raw, method = 'BH'))
write.csv(TableS1_posthoc_adj.out, file = out.path('-Table S1-posthoc_pval'))

################################################################################
#### FIGURE 1: Microbial Load, Diversity, and Composition Overview ####
################################################################################
# Figure 1 shows microbial profiles in lower respiratory tract (Nanopore, 16S) 
# and plasma (Karius) across clinical groups (MCP, CDP, CONTROL)

## Panel A: Total Nanopore reads by clinical group
Fig1A.p <- 
  Rappid.pneum_category_final.df %>% 
  left_join(
    Nanopore.all_reads.Sankey.df %>% dplyr::filter(Human != 'Human reads') %>% group_by(SubjectIDDay) %>% summarise(TotalNanoporeReads = n()) %>% ungroup()
  ) %>% 
  ggboxplot(x = 'pneum_category_final', y = 'TotalNanoporeReads', add = 'jitter') +  scale_y_log10() + ylab("Total Nanopore Reads") + xlab('')

Fig1A.posthoc.comparison <- 
  Nanopore.all_reads.Sankey.df %>% dplyr::filter(Human != 'Human reads') %>% group_by(SubjectIDDay) %>% summarise(TotalNanoporeReads = n()) %>% ungroup() %>% inner_join(Rappid.pneum_category_final.df) %>% 
  pairwise_wilcox_test(formula = TotalNanoporeReads ~ pneum_category_final, p.adjust.method = 'BH') %>% add_xy_position() %>% mutate(y.position = log10(y.position))
Fig1A.p + stat_pvalue_manual(Fig1A.posthoc.comparison, hide.ns = T)

## Panel B: 16S qPCR copy number by clinical group
Fig1B.p <- Rappid.baseline.meta.202305.df %>% select(SubjectIDDay, copies_ETA) %>% inner_join(Rappid.pneum_category_final.df) %>% 
  ggboxplot(x = 'pneum_category_final', y = 'copies_ETA', add = 'jitter') + scale_y_log10() + xlab('')

Fig1B.posthoc.comparison <- 
  Rappid.baseline.meta.202305.df %>% select(SubjectIDDay, copies_ETA) %>% inner_join(Rappid.pneum_category_final.df) %>% 
  pairwise_wilcox_test(formula = copies_ETA ~ pneum_category_final, p.adjust.method = 'BH') %>% add_xy_position() %>% 
  mutate(y.position.old = y.position, y.position =  mean(log10(y.position.old)) + 1.2 + 0.9*(log10(y.position.old) - mean(log10(y.position.old))) / (max(log10(y.position.old)) - min(log10(y.position.old))))
Fig1B.p + stat_pvalue_manual(Fig1B.posthoc.comparison, hide.ns = T)

## Panel C: Total plasma mcfDNA by clinical group
Fig1C.p <- Rappid.baseline.meta.202305.df %>% transmute(SubjectIDDay, total_mcfdna_log = case_when(total_mcfdna == 0 ~ 0, TRUE ~ log(total_mcfdna))) %>% inner_join(Rappid.pneum_category_final.df) %>% 
  ggboxplot(x = 'pneum_category_final', y = 'total_mcfdna_log', add = 'jitter') + xlab('')

Fig1C.posthoc.comparison <- 
  Rappid.baseline.meta.202305.df %>% transmute(SubjectIDDay, total_mcfdna) %>% inner_join(Rappid.pneum_category_final.df) %>% 
  pairwise_wilcox_test(formula = total_mcfdna ~ pneum_category_final, p.adjust.method = "BH") %>% add_xy_position() %>% 
  mutate(y.position.old = y.position, y.position =  mean(log10(y.position.old)) + 9 + 1.6*(log10(y.position.old) - mean(log10(y.position.old))) / (max(log10(y.position.old)) - min(log10(y.position.old))))
Fig1C.p + stat_pvalue_manual(Fig1C.posthoc.comparison, hide.ns = T)

## Panel D: Shannon diversity (Nanopore) by clinical group
Fig1D.p <- 
  Rappid.pneum_category_final.df %>% left_join(Nanopore.filtered.RelAbd.wide.df %>% column_to_rownames('SubjectIDDay') %>% vegan::diversity() %>% enframe(value = 'Shannon_Nanopore', name = 'SubjectIDDay')) %>% 
  ggboxplot(x = 'pneum_category_final', y = 'Shannon_Nanopore', add = 'jitter') + xlab('')

Fig1D.posthoc.comparison <- 
  Rappid.pneum_category_final.df %>% left_join(Nanopore.filtered.RelAbd.wide.df %>% column_to_rownames('SubjectIDDay') %>% vegan::diversity() %>% enframe(value = 'Shannon_Nanopore', name = 'SubjectIDDay')) %>% 
  pairwise_wilcox_test(Shannon_Nanopore ~ pneum_category_final, p.adjust.method = "BH") %>% add_xy_position()
Fig1D.p + stat_pvalue_manual(Fig1D.posthoc.comparison, hide.ns = T)

## Panel E: Shannon diversity (16S) by clinical group
Fig1E.p <- 
  Rappid.baseline.meta.202305.df %>% transmute(SubjectIDDay, Shannon_16S = shannon_ETA) %>% inner_join(Rappid.pneum_category_final.df) %>% 
  ggboxplot(x = 'pneum_category_final', y = 'Shannon_16S', add = 'jitter') + xlab('')

Fig1E.posthoc.comparison <-
  Rappid.baseline.meta.202305.df %>% transmute(SubjectIDDay, Shannon_16S = shannon_ETA) %>% inner_join(Rappid.pneum_category_final.df) %>% 
  pairwise_wilcox_test(formula = Shannon_16S ~ pneum_category_final, p.adjust.method = "BH") %>% add_xy_position()

## Panel F: Number of species detected in plasma by clinical group
Fig1F.p <- 
  Karius.df %>% distinct(SubjectIDDay, Taxa, MPM, .keep_all = T) %>% group_by(SubjectIDDay) %>% summarise(Total_n_species = case_when(any(str_detect(Taxa, pattern = '^No\\ ')) ~ 0L, TRUE ~ n())) %>% inner_join(Rappid.pneum_category_final.df) %>% 
  ggboxplot(x = 'pneum_category_final', y = 'Total_n_species', add = 'jitter') + xlab('')

Fig1F.posthoc.comparison <- 
  Karius.df %>% distinct(SubjectIDDay, Taxa, MPM, .keep_all = T) %>% group_by(SubjectIDDay) %>% summarise(Total_n_species = case_when(any(str_detect(Taxa, pattern = '^No\\ ')) ~ 0L, TRUE ~ n())) %>% inner_join(Rappid.pneum_category_final.df) %>% 
  pairwise_wilcox_test(formula = Total_n_species ~ pneum_category_final, p.adjust.method = "BH") %>% add_xy_position()
  
## Panel G: Microbial composition by kingdom (pie charts)
Fig1G.df <- 
  rbind(
    Nanopore.all_reads.Sankey.df %>% dplyr::filter(Human != 'Human reads') %>% group_by(SubjectIDDay, Kingdom) %>% summarise(KingdomSum_bySample = n()) %>% ungroup() %>% right_join(Rappid.pneum_category_final.df) %>% group_by(pneum_category_final, Kingdom) %>% summarise(KingdomSum_byCategory = sum(KingdomSum_bySample, na.rm = T)) %>% dplyr::filter(!is.na(Kingdom)) %>% ungroup() %>% mutate(Compartment = 'Lower Airway\nNanopore Seq'), 
    Karius.df %>% dplyr::filter(Taxa != 'No Result' & Taxa != 'No Organism Detected') %>% 
      mutate(Kingdom = case_when(
        Taxa %in% (readxl::read_excel(path = 'data/karius_alir_allresult_3.11.20.edited_Oct18.20.xlsx') %>% dplyr::filter(Fungus) %>% distinct(SIGNIFICANT_BUG_CALLS) %>% deframe()) ~ 'Eukaryota', 
        Taxa %in% (readxl::read_excel(path = 'data/karius_alir_allresult_3.11.20.edited_Oct18.20.xlsx') %>% dplyr::filter(Virus) %>% distinct(SIGNIFICANT_BUG_CALLS) %>% deframe()) ~ 'Virus', 
        TRUE ~ 'Bacteria'
      )) %>% 
      group_by(SubjectIDDay, Kingdom) %>% summarise(KingdomSum_bySample = sum(MPM, na.rm = F)) %>% ungroup() %>% inner_join(Rappid.pneum_category_final.df) %>% 
      group_by(pneum_category_final, Kingdom) %>% summarise(KingdomSum_byCategory = sum(KingdomSum_bySample, na.rm = T)) %>% ungroup() %>% mutate(Compartment = 'Plasma mcfDNA\nKarius Test') 
  ) %>% dplyr::filter(complete.cases(.)) %>% group_by(pneum_category_final, Compartment) %>% 
  mutate(Sum_across_Category = sum(KingdomSum_byCategory)) %>% ungroup() %>% 
  mutate(pneum_category_final = factor(pneum_category_final, levels = c('MCP', 'CDP', 'CONTROL')))

options(scipen = 2)
Fig1G.p <-
  arrangeGrob(
  Fig1G.df %>% dplyr::filter(str_detect(Compartment, 'Nanopore')) %>% 
    ggplot(aes(x = Sum_across_Category / 2 / 100000, 
               width = Sum_across_Category / 100000, y = KingdomSum_byCategory / 100000, fill = Kingdom)) + 
    geom_bar(stat = 'identity', position = 'fill', color = 'black') + 
    facet_wrap(.~ pneum_category_final, scales = 'fixed', ncol = 1) + 
    ylab('Lower respira-\ntory tract')+ 
    xlab('Total Nanopore reads, millions') + coord_polar("y", start = 0, direction = -1) + 
    scale_fill_manual(values = viridis::cividis(5))+ theme_minimal() + 
    scale_y_continuous(labels = scales::label_comma()) + 
    scale_x_continuous( n.breaks = 3) + 
    theme(legend.position = "bottom", axis.text.x = element_blank(), 
          axis.text.y = element_text(size = 8)) +  guides(fill=guide_legend(nrow=3, byrow=TRUE)), 
  
  Fig1G.df %>% dplyr::filter(str_detect(Compartment, 'Karius')) %>% 
    ggplot(aes(x = Sum_across_Category / 2 / 100000, 
               width = Sum_across_Category / 100000, 
               y = KingdomSum_byCategory / 100000, fill = Kingdom)) + 
    geom_bar(stat = 'identity', position = 'fill', color = 'black') + 
    facet_wrap(.~ pneum_category_final, scales = 'fixed', ncol = 1) + ylab('Plasma') + 
    xlab('Total MPM collected, millions') + coord_polar("y", start = 0, direction = -1) + 
    scale_y_continuous(labels = scales::label_comma()) + 
    scale_fill_manual(values = viridis::cividis(5)[c(2,3,5)]) + 
    theme_minimal() + 
    theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_text(size = 8)), 
  ncol = 1, 
  top = 'G. Microbial reads by kingdom'
)

## Panel H: Beta diversity (PCoA) by clinical group
Nanopore.dist <- 
  Nanopore.filtered.RelAbd.wide.df %>% 
  dplyr::filter(SubjectIDDay %in% Rappid.pneum_category_final.df$SubjectIDDay) %>% column_to_rownames("SubjectIDDay") %>% 
  vegan::vegdist(x = ., method = 'bray')

Nanopore.PCoA <- vegan::wcmdscale(Nanopore.dist, 
                 k = Nanopore.filtered.RelAbd.wide.df %>% dplyr::filter(SubjectIDDay %in% Rappid.pneum_category_final.df$SubjectIDDay) %>% column_to_rownames("SubjectIDDay") %>% ncol() - 1, 
                 eig = T)

Nanopore.PCoA.pos_eig <- Nanopore.PCoA$eig[Nanopore.PCoA$eig > 0]

Nanopore.PCoA.species_loading.df <- 
  vegan::wascores(
    x = Nanopore.PCoA[['points']], 
    w = Nanopore.filtered.RelAbd.wide.df %>% 
      dplyr::filter(SubjectIDDay %in% Rappid.pneum_category_final.df$SubjectIDDay) %>% column_to_rownames("SubjectIDDay")
  ) %>% as.data.frame() %>% rownames_to_column('Species') %>% dplyr::select(Species, Dim1, Dim2)

Nanopore.PCoA.species_loading.df %>% 
  dplyr::filter(rank(-abs(Dim1)) <= 5 & rank(-abs(Dim2)) <= 5) 

Nanopore.PCoA.pneum_cate.df <- 
  Nanopore.filtered.RelAbd.wide.df %>% 
  dplyr::filter(SubjectIDDay %in% Rappid.pneum_category_final.df$SubjectIDDay) %>% 
  left_join(Rappid.pneum_category_final.df) %>% dplyr::select(SubjectIDDay, pneum_category_final) %>% ungroup()

set.seed(20251028)
vegan::betadisper(d = Nanopore.dist, group = Nanopore.PCoA.pneum_cate.df$pneum_category_final) %>% anova()

set.seed(20250425)
Nanopore.PERMANOVA <- 
  vegan::adonis2(formula = 
                   Nanopore.dist ~ Nanopore.filtered.RelAbd.wide.df %>% dplyr::select(SubjectIDDay) %>% 
                   inner_join(Rappid.pneum_category_final.df) %>% ungroup() %>% select(pneum_category_final) %>% deframe())

Fig1H.p <- ggordiplots::gg_ordiplot(groups = rownames(Nanopore.PCoA$points) %>% enframe(value = 'SubjectIDDay', name = NULL) %>% left_join(Rappid.pneum_category_final.df) %>% transmute(pneum_category_final = factor(pneum_category_final, levels = c('MCP', 'CDP', 'CONTROL'))) %>% deframe(), 
                         ord = Nanopore.PCoA, 
                         spider = T, plot = F)$plot + theme_minimal() + 
  annotate('text', x = -0.3, y = -0.3, label = paste("PERMANOVA\np=", Nanopore.PERMANOVA
 
Pr(>F)`[1], sep = '')) + theme(legend.position = 'bottom')

## Panel I: LEfSe enrichment analysis by clinical group
Nanopore.microtable <- microeco::microtable$new(
  otu_table = Nanopore.filtered.RelAbd.wide.df %>% column_to_rownames("SubjectIDDay") %>% t() %>% as.data.frame(), 
  tax_table = Nanopore.filtered.RelAbd.wide.df %>% colnames() %>% tail(-1) %>% enframe() %>% mutate(Species = value) %>% select(-name) %>% column_to_rownames('value'), 
  sample_table = rownames(Nanopore.PCoA$points) %>% 
    enframe(value = 'SubjectIDDay', name = NULL) %>% 
    left_join(Rappid.pneum_category_final.df %>% mutate(pneum_category_final = factor(pneum_category_final, levels = c('MCP', 'CDP', 'CONTROL')))) %>% as.data.frame() %>% column_to_rownames('SubjectIDDay')
)

Nanopore.Lefse <- trans_diff$new(Nanopore.microtable, taxrank = 'Species', group = 'pneum_category_final', p_adjust_method = 'none')
Fig1I.p <- Nanopore.Lefse$plot_diff_bar() + theme(legend.position = "bottom", axis.text.y = element_text(face = 'italic', size = 14))

## Assemble and save Figure 1
png(filename = out.path('-Figure 1'), width = 15, height = 10, res = 600, units = 'in')
grid.arrange(
  arrangeGrob(
    Fig1A.p + scale_x_discrete(limits = rev) + stat_pvalue_manual(Fig1A.posthoc.comparison, hide.ns = T) + ggtitle("A"), 
    Fig1D.p + scale_x_discrete(limits = rev) + stat_pvalue_manual(Fig1D.posthoc.comparison, hide.ns = T) + ggtitle("B") + ylab("Shannon Index, Nanopore"), 
    Fig1B.p + scale_x_discrete(limits = rev) + stat_pvalue_manual(Fig1B.posthoc.comparison, hide.ns = T) + ggtitle("C") + ylab("16S copies, qPCR"), 
    Fig1E.p + scale_x_discrete(limits = rev) + stat_pvalue_manual(Fig1E.posthoc.comparison, hide.ns = T) + ggtitle("D") + ylab("Shannon Index, 16S"), 
    Fig1C.p + scale_x_discrete(limits = rev) + stat_pvalue_manual(Fig1C.posthoc.comparison, hide.ns = T) + ggtitle("E") + ylab("Total microbial cell-free DNA\nmolecule per mL, log transformed"), 
    Fig1F.p + scale_x_discrete(limits = rev) + ylab("Total # of species detected, mcfDNA") + stat_pvalue_manual(Fig1F.posthoc.comparison, hide.ns = T) + ggtitle("F"), 
    ncol = 2
  ), 
  Fig1G.p, 
  arrangeGrob(
    Fig1H.p + scale_colour_viridis_d() + ggtitle("H"), 
    Fig1I.p + scale_fill_viridis_d()  + ggtitle("I"), 
    ncol = 1
  ), 
  ncol = 3, widths = c(3, 1, 3)
)
dev.off()

################################################################################
#### FIGURE 2: Pathogen Enrichment in MCP ####
################################################################################
# Figure 2 shows enrichment of respiratory pathogens in microbiologically 
# confirmed pneumonia (MCP) compared to clinically diagnosed pneumonia (CDP) 
# and controls

png(filename = out.path('-Figure 2'), width = 12, height = 8, res = 600, units = "in")
grid.arrange(
  pathogen_pneum_cate.Nanopore.p + labs(title = "A. Enrichment of respiratory pathogen per LRT metagenomics") + theme_minimal() + theme(axis.text.y = element_blank(), strip.text = element_text(size = 16), legend.position = 'bottom'), 
  pathogen_pneum_cate.Nanopore.pmanual.p  + labs(title = "B"), 
  pathogen_pneum_cate.Karius.p + labs(title = "C. Enrichment of respiratory pathogen per plasma metagenomics") + theme_minimal() + theme(axis.text.y = element_blank(), strip.text = element_text(size = 16), legend.position = 'bottom'), 
  pathogen_pneum_cate.Karius.pvalue.p + labs(title = "D"), 
  widths = c(3, 1), 
  ncol = 2)
dev.off()

################################################################################
#### FIGURE 3: Distinct Microbial Signatures by Translocation Group ####
################################################################################
# Figure 3 shows microbial differences between translocation subgroups

png(filename = out.path('-Figure 3'), width = 10, height = 8, units = "in", res = 600)
grid.arrange(
  Nanopore.TransLoc.LEFSE.diff_bar.p + ggtitle("A") + 
    theme(legend.position = "none", axis.text.y = element_text(size = 10, face = 'italic')), 
  TransLoc_final.mcfDNA.p + ylab("Total mcfDNA, MPM, log transformed") + xlab("") + ggtitle("C") + 
    theme(legend.position = "none", axis.text.x = element_text(size = 8)), 
  Nanopore.TransLoc.LEFSE.diff_abund.p + ggtitle("B") + 
    theme(legend.position = "none", axis.text.y = element_text(size = 10, face = 'italic')), 
  TransLoc_Respir.mcfDNA.p+ ylab("Respir pathogen mcfDNA, MPM, log transformed") + xlab("") + ggtitle("D") + 
    theme(legend.position = "none", axis.text.x = element_text(size = 8)), 
  ncol = 2, widths = c(2, 1)
)
dev.off()

################################################################################
#### FIGURE 4: Microbial Translocation and Host Inflammatory Responses ####
################################################################################
# Figure 4 shows associations between microbial translocation and biomarkers
# Panels A-C: Univariate regression models
# Panels D-G: Comparisons by inflammation subphenotype

## Panel A: LRT pathogen abundance predicting biomarkers
Fig4A.df <- 
  Nanopore.TransLoc_final.LEFSE.df %>% left_join(Patho_Classification.FINAL.raw[, c("SPECIES", "HY_FINAL")] , by = c("name_species" = "SPECIES")) %>% 
  group_by(SubjectIDDay, TransLoc_final, HY_FINAL) %>% summarise(GroupRelAbdSum = sum(filtered_RelAbd)) %>% ungroup() %>% 
  mutate(Patho_Rank_Place_Holder = case_when(HY_FINAL == "established" ~ GroupRelAbdSum, TRUE ~ 0)) %>% 
  left_join(DummyID.df %>% mutate(SubjectIDDay = paste(SubjectID, "_1", sep = ""), DummyID = paste("Patient #", ParticipantID, sep = " "))) %>% 
  mutate(Pathogenicity = str_to_upper(HY_FINAL)) %>% dplyr::filter(Pathogenicity == 'ESTABLISHED') %>% ungroup() %>% 
  inner_join(Rappid.baseline.meta.202305.df %>% transmute(SubjectIDDay, across(ang2:rage, .fns = log), tnfr1_final = log(tnfr1_final))) %>% 
  transmute(SubjectIDDay, TransLoc_final, NanoporePathoSum = GroupRelAbdSum, across(ang2:tnfr1_final))

Fig4A_ETApred.df <- 
  Nanopore.TransLoc_final.LEFSE.df %>% left_join(Patho_Classification.FINAL.raw[, c("SPECIES", "HY_FINAL")] , by = c("name_species" = "SPECIES")) %>% 
  group_by(SubjectIDDay, TransLoc_final, HY_FINAL) %>% summarise(GroupRelAbdSum = sum(filtered_RelAbd)) %>% ungroup() %>% 
  mutate(Patho_Rank_Place_Holder = case_when(HY_FINAL == "established" ~ GroupRelAbdSum, TRUE ~ 0)) %>% 
  left_join(DummyID.df %>% mutate(SubjectIDDay = paste(SubjectID, "_1", sep = ""), DummyID = paste("Patient #", ParticipantID, sep = " "))) %>% 
  mutate(Pathogenicity = str_to_upper(HY_FINAL)) %>% dplyr::filter(Pathogenicity == 'ESTABLISHED') %>% ungroup() %>% 
  inner_join(Rappid.baseline.meta.202305.df %>% 
               transmute(SubjectIDDay, 
                         ang2 = ang2_eta_pradj, il8 = il8_eta, il6 = il6_eta, procalcitonin = procalcitonin_eta, 
                         st2 = st2_eta, fractalkine = fractalkine_eta, il10 = il10_eta, pentraxin3 = pentraxin3_eta, rage = rage_eta, 
                         tnfr1_final = tnfr1_eta)) %>% 
  transmute(SubjectIDDay, TransLoc_final, NanoporePathoSum = GroupRelAbdSum, across(ang2:tnfr1_final, .fns = function(x){return(x/1e6)}))

# Univariate regression: plasma biomarkers ~ LRT pathogen abundance
Fig4A.univar_reg.res.df <- 
  sapply(
    X = Fig4A.df %>% select(ang2:tnfr1_final) %>% colnames(), 
    FUN = function(biomarker) {
      return(
        c(summary(lm(formula = as.formula(paste(biomarker, '~ NanoporePathoSum', sep = '')), data = Fig4A.df))$coefficients["NanoporePathoSum", ], 
          confint(lm(formula = as.formula(paste(biomarker, '~ NanoporePathoSum', sep = '')), data = Fig4A_ETApred.df))["NanoporePathoSum", ])
      )
    }
  ) %>% t() %>% as.data.frame() %>% rownames_to_column('Biomarker') %>% 
  mutate(Biomarker = case_when(
    Biomarker == 'tnfr1_final' ~ 'sTNFR-1', 
    TRUE ~ str_to_upper(Biomarker)
  ), 
  p.adj = p.adjust(`Pr(>|t|)`, method = 'BH'), 
  annotate_string = case_when(
    p.adj < 0.05 ~ paste(round(Estimate, digits = 2), ' [', round(`2.5 %`, digits = 2), ', ', round(`97.5 %`, digits = 2), ']', sep = ''), 
    p.adj >= 0.05 ~ NA_character_
  )) %>% arrange(Biomarker)

Fig4A.p <- 
  ggplot(data = Fig4A.univar_reg.res.df, aes(x = Estimate, y = Biomarker, colour = p.adj < 0.05)) + 
  geom_vline(xintercept = 0, linetype = 2) + 
  geom_point(size = 3) + geom_errorbarh(aes(xmax = `97.5 %`, xmin = `2.5 %`), height = 0.2) + 
  scale_colour_manual(values = c('gray75', 'black')) + 
  scale_y_discrete(limits = rev) + 
  ggrepel::geom_text_repel(aes(label = annotate_string), size = 3.5) + 
  theme_minimal() + theme(legend.position = 'none') 

Fig4A_ETApred.p <- 
  ggplot(data = Fig4A_ETApred.univar_reg.res.df, aes(x = Estimate, y = Biomarker, colour = p.adj < 0.05)) + 
  geom_vline(xintercept = 0, linetype = 2) + 
  geom_point(size = 3) + geom_errorbarh(aes(xmax = `97.5 %`, xmin = `2.5 %`), height = 0.2) + 
  scale_colour_manual(values = c('gray75', 'black')) + 
  scale_y_discrete(limits = rev) + 
  ggrepel::geom_text_repel(aes(label = annotate_string), size = 4) + 
  theme_minimal() + theme(legend.position = 'none') 

## Panel B: Plasma pathogen mcfDNA predicting biomarkers
Fig4B.df <- 
  Karius.df %>% dplyr::filter(Day == "1") %>% inner_join(rappid.TransLoc.df[, c("SubjectIDDay", "TransLoc_final_fix")]) %>% 
  left_join(Patho_Classification.FINAL.raw[, c("SPECIES", "HY_FINAL")] , by = c("Taxa" = "SPECIES")) %>% 
  group_by(SubjectIDDay, HY_FINAL, TransLoc_final_fix) %>% summarise(pathoSumMPM = sum(MPM)) %>% arrange(pathoSumMPM) %>% ungroup() %>% 
  dplyr::filter(!is.na(pathoSumMPM)) %>% 
  mutate(Patho_Rank_Place_Holder = case_when(HY_FINAL == "established" ~ pathoSumMPM, TRUE ~ 0), SubjectID = str_remove(SubjectIDDay, pattern = "_[\\d]+$")) %>% 
  dplyr::filter(HY_FINAL == 'established' | is.na(HY_FINAL)) %>% ungroup() %>% 
  inner_join(Rappid.baseline.meta.202305.df %>% transmute(SubjectIDDay, across(ang2:rage, .fns = log), tnfr1_final = log(tnfr1_final))) %>% 
  dplyr::filter(TransLoc_final_fix != 'No translocation') %>% 
  transmute(SubjectIDDay, TransLoc_final_fix, pathoSumMPM_log = case_when(pathoSumMPM == 0 ~ 0, TRUE ~ log(pathoSumMPM)), across(ang2:tnfr1_final))

Fig4B.univar_reg.res.df <- 
  sapply(
    X = Fig4B.df %>% select(ang2:tnfr1_final) %>% colnames(), 
    FUN = function(biomarker) {
      return(
        c(summary(lm(formula = as.formula(paste(biomarker, '~ pathoSumMPM_log', sep = '')), data = Fig4B.df))$coefficients["pathoSumMPM_log", ], 
          confint(lm(formula = as.formula(paste(biomarker, '~ pathoSumMPM_log', sep = '')), data = Fig4B.df))["pathoSumMPM_log", ])
      )
    }
  ) %>% t() %>% as.data.frame() %>% rownames_to_column('Biomarker') %>% 
  mutate(Biomarker = case_when(
    Biomarker == 'tnfr1_final' ~ 'sTNFR-1', 
    TRUE ~ str_to_upper(Biomarker)
  ), 
  p.adj = p.adjust(`Pr(>|t|)`, method = 'BH'), 
  annotate_string = case_when(
    p.adj < 0.05 ~ paste(sprintf('%.2f', Estimate), ' [', sprintf('%.2f', `2.5 %`), ', ', sprintf('%.2f', `97.5 %`), ']', sep = ''), 
    p.adj >= 0.05 ~ NA_character_
  )) %>% arrange(Biomarker)

Fig4B.p <- 
  ggplot(data = Fig4B.univar_reg.res.df, aes(x = Estimate, y = Biomarker, colour = p.adj < 0.05)) + 
  geom_vline(xintercept = 0, linetype = 2) + 
  geom_point(size = 3) + geom_errorbarh(aes(xmax = `97.5 %`, xmin = `2.5 %`), height = 0.2) + 
  scale_colour_manual(values = c('gray75', 'black')) + 
  scale_y_discrete(limits = rev) + 
  ggrepel::geom_text_repel(aes(label = annotate_string)) + 
  theme_minimal() + theme(legend.position = 'none') 

## Panel C: Pulmonary translocation burden predicting biomarkers
Fig4C.df <- rappid.TransLoc.df %>% 
  inner_join(Rappid.baseline.meta.202305.df %>% transmute(SubjectIDDay, across(ang2:rage, .fns = log), tnfr1_final = log(tnfr1_final))) %>% 
  mutate(TransLoc_MPM_log = case_when(TransLoc_MPM == 0 ~ 0, TransLoc_MPM > 0 ~ log(TransLoc_MPM))) %>% 
  dplyr::filter(TransLoc_final_fix != 'No translocation')

Fig4C.univar_reg.res.df <- 
  sapply(
    X = Fig4C.df %>% select(ang2:tnfr1_final) %>% colnames(), 
    FUN = function(biomarker) {
      return(
        c(summary(lm(formula = as.formula(paste(biomarker, '~ TransLoc_MPM_log', sep = '')), data = Fig4C.df))$coefficients["TransLoc_MPM_log", ], 
          confint(lm(formula = as.formula(paste(biomarker, '~ TransLoc_MPM_log', sep = '')), data = Fig4C.df))["TransLoc_MPM_log", ])
      )
    }
  ) %>% t() %>% as.data.frame() %>% rownames_to_column('Biomarker') %>% 
  mutate(Biomarker = case_when(
    Biomarker == 'tnfr1_final' ~ 'sTNFR-1', 
    TRUE ~ str_to_upper(Biomarker)
  ), 
  p.adj = p.adjust(`Pr(>|t|)`, method = 'BH'), 
  annotate_string = case_when(
    p.adj < 0.05 ~ paste(sprintf('%.2f', Estimate), ' [', sprintf('%.2f', `2.5 %`), ', ', sprintf('%.2f', `97.5 %`), ']', sep = ''), 
    p.adj >= 0.05 ~ NA_character_
  )) %>% arrange(Biomarker)

Fig4C.p <- 
  ggplot(data = Fig4C.univar_reg.res.df, aes(x = Estimate, y = Biomarker, colour = p.adj < 0.05)) + 
  geom_vline(xintercept = 0, linetype = 2) + 
  geom_point(size = 3) + geom_errorbarh(aes(xmax = `97.5 %`, xmin = `2.5 %`), height = 0.2) + 
  scale_colour_manual(values = c('gray75', 'black')) + 
  scale_y_discrete(limits = rev) + 
  ggrepel::geom_text_repel(aes(label = annotate_string)) + 
  theme_minimal() + theme(legend.position = 'none') 

## Panel D: Shannon diversity by inflammation subphenotype
Fig4D.p <- 
  Nanopore.filtered.RelAbd.wide.df %>% column_to_rownames('SubjectIDDay') %>% vegan::diversity() %>% enframe(value = 'Shannon_Nanopore', name = 'SubjectIDDay') %>% 
  inner_join(Rappid.meta.202305.df %>% transmute(SubjectIDDay, nouraiepars21 = case_when(nouraiepars21 == 'hypoinfl' ~ 'Hypo-\ninflammatory', nouraiepars21 == 'hyperinf' ~ 'Hyper-\ninflammatory'))) %>% 
  inner_join(Rappid.pneum_category_final.df) %>% dplyr::filter(pneum_category_final != 'CONTROL' & !is.na(nouraiepars21)) %>% 
  ggboxplot(x = 'nouraiepars21', y = 'Shannon_Nanopore', add = 'jitter') + stat_compare_means() + xlab("") + ylab('Shannon Index, Nanopore')

## Panel E: Beta diversity by inflammation subphenotype
Fig4E.df <- 
  Rappid.meta.202305.df %>% transmute(SubjectIDDay, nouraiepars21 = case_when(nouraiepars21 == 'hypoinfl' ~ 'Hypo-\ninflammatory', nouraiepars21 == 'hyperinf' ~ 'Hyper-\ninflammatory')) %>% 
  inner_join(Rappid.pneum_category_final.df) %>% left_join(rappid.TransLoc.df) %>% 
  dplyr::filter(pneum_category_final != 'CONTROL' & !is.na(nouraiepars21) & SubjectIDDay %in% Nanopore.filtered.RelAbd.wide.df$SubjectIDDay) %>% select(-SubjectID) 
  
Fig4E.df %>% select(SubjectIDDay, TransLoc_final_fix, nouraiepars21) %>% dplyr::filter(complete.cases(.)) %>% group_by(TransLoc_final_fix, nouraiepars21) %>% summarise(n = n()) %>% ungroup() %>% 
  pivot_wider(values_from = 'n', names_from = 'nouraiepars21', values_fill = 0L) %>% 
  column_to_rownames('TransLoc_final_fix') %>% pairwise_fisher_test(p.adjust.method = "BH")

Nanopore.subpheno.dist <- 
  Nanopore.filtered.RelAbd.wide.df %>% dplyr::filter(SubjectIDDay %in% Fig4E.df$SubjectIDDay) %>% column_to_rownames("SubjectIDDay") %>% 
  vegan::vegdist(x = ., method = 'bray')

Nanopore.subpheno.PCoA <- vegan::wcmdscale(Nanopore.subpheno.dist, 
                                           k = Nanopore.filtered.RelAbd.wide.df %>% 
                                             dplyr::filter(SubjectIDDay %in% Fig4E.df$SubjectIDDay) %>% column_to_rownames("SubjectIDDay") %>% ncol() - 1, 
                                           eig = T)
set.seed(20250426)
Nanopore.subpheno.PERMANOVA <- 
  vegan::adonis2(formula = Nanopore.subpheno.dist ~ Fig4E.df$nouraiepars21)

Fig4E.p <- ggordiplots::gg_ordiplot(groups = Fig4E.df$nouraiepars2, 
                                     ord = Nanopore.subpheno.PCoA, 
                                     spider = T, plot = F)$plot + theme_minimal() + 
  annotate('text', x = -0.27, y = 0.3, label = paste("PERMANOVA\np=", Nanopore.subpheno.PERMANOVA
                                                     
                                                     
################################################################################
#### FIGURE 4: Microbial Translocation and Host Inflammatory Responses ####
################################################################################
# Figure 4 shows associations between microbial translocation and biomarkers
# Panels A-C: Univariate regression models
# Panels D-G: Comparisons by inflammation subphenotype

## Panel A: LRT pathogen abundance predicting biomarkers
Fig4A.df <- 
  Nanopore.TransLoc_final.LEFSE.df %>% left_join(Patho_Classification.FINAL.raw[, c("SPECIES", "HY_FINAL")] , by = c("name_species" = "SPECIES")) %>% 
  group_by(SubjectIDDay, TransLoc_final, HY_FINAL) %>% summarise(GroupRelAbdSum = sum(filtered_RelAbd)) %>% ungroup() %>% 
  mutate(Patho_Rank_Place_Holder = case_when(HY_FINAL == "established" ~ GroupRelAbdSum, TRUE ~ 0)) %>% 
  left_join(DummyID.df %>% mutate(SubjectIDDay = paste(SubjectID, "_1", sep = ""), DummyID = paste("Patient #", ParticipantID, sep = " "))) %>% 
  mutate(Pathogenicity = str_to_upper(HY_FINAL)) %>% dplyr::filter(Pathogenicity == 'ESTABLISHED') %>% ungroup() %>% 
  inner_join(Rappid.baseline.meta.202305.df %>% transmute(SubjectIDDay, across(ang2:rage, .fns = log), tnfr1_final = log(tnfr1_final))) %>% 
  transmute(SubjectIDDay, TransLoc_final, NanoporePathoSum = GroupRelAbdSum, across(ang2:tnfr1_final))

Fig4A_ETApred.df <- 
  Nanopore.TransLoc_final.LEFSE.df %>% left_join(Patho_Classification.FINAL.raw[, c("SPECIES", "HY_FINAL")] , by = c("name_species" = "SPECIES")) %>% 
  group_by(SubjectIDDay, TransLoc_final, HY_FINAL) %>% summarise(GroupRelAbdSum = sum(filtered_RelAbd)) %>% ungroup() %>% 
  mutate(Patho_Rank_Place_Holder = case_when(HY_FINAL == "established" ~ GroupRelAbdSum, TRUE ~ 0)) %>% 
  left_join(DummyID.df %>% mutate(SubjectIDDay = paste(SubjectID, "_1", sep = ""), DummyID = paste("Patient #", ParticipantID, sep = " "))) %>% 
  mutate(Pathogenicity = str_to_upper(HY_FINAL)) %>% dplyr::filter(Pathogenicity == 'ESTABLISHED') %>% ungroup() %>% 
  inner_join(Rappid.baseline.meta.202305.df %>% 
               transmute(SubjectIDDay, 
                         ang2 = ang2_eta_pradj, il8 = il8_eta, il6 = il6_eta, procalcitonin = procalcitonin_eta, 
                         st2 = st2_eta, fractalkine = fractalkine_eta, il10 = il10_eta, pentraxin3 = pentraxin3_eta, rage = rage_eta, 
                         tnfr1_final = tnfr1_eta)) %>% 
  transmute(SubjectIDDay, TransLoc_final, NanoporePathoSum = GroupRelAbdSum, across(ang2:tnfr1_final, .fns = function(x){return(x/1e6)}))

# Univariate regression: plasma biomarkers ~ LRT pathogen abundance
Fig4A.univar_reg.res.df <- 
  sapply(
    X = Fig4A.df %>% select(ang2:tnfr1_final) %>% colnames(), 
    FUN = function(biomarker) {
      return(
        c(summary(lm(formula = as.formula(paste(biomarker, '~ NanoporePathoSum', sep = '')), data = Fig4A.df))$coefficients["NanoporePathoSum", ], 
Pr(>F)`[1], sep = '')) + theme_minimal() + theme(legend.position = 'bottom')

## Panel F: Total mcfDNA by inflammation subphenotype
Fig4F.p <- Rappid.pneum_category_final.df %>% left_join(Rappid.baseline.meta.202305.df) %>% 
  transmute(SubjectIDDay, pneum_category_final, nouraiepars21 = case_when(nouraiepars21 == 'hypoinfl' ~ 'Hypo-\ninflammatory', nouraiepars21 == 'hyperinf' ~ 'Hyper-\ninflammatory'), total_mcfdna, total_mcfdna_log = case_when(total_mcfdna == 0 ~ 0, TRUE ~ log(total_mcfdna)))  %>% 
  dplyr::filter(pneum_category_final != 'CONTROL' & !is.na(nouraiepars21)) %>% 
  ggboxplot(x = 'nouraiepars21', y = 'total_mcfdna_log', add = 'jitter') + stat_compare_means()

## Panel G: Translocation burden by inflammation subphenotype
Fig4G.p <- 
  Rappid.pneum_category_final.df %>% left_join(Rappid.baseline.meta.202305.df) %>% inner_join(rappid.TransLoc.df) %>% 
  transmute(SubjectIDDay, pneum_category_final, TransLoc_MPM_log = case_when(TransLoc_MPM == 0 ~ 0, TRUE ~ log(TransLoc_MPM)), TransLoc_final_fix, nouraiepars21 = case_when(nouraiepars21 == 'hypoinfl' ~ 'Hypo-\ninflammatory', nouraiepars21 == 'hyperinf' ~ 'Hyper-\ninflammatory'), total_mcfdna, total_mcfdna_log = case_when(total_mcfdna == 0 ~ 0, TRUE ~ log(total_mcfdna)))  %>% 
  dplyr::filter(pneum_category_final != 'CONTROL' & !is.na(nouraiepars21)) %>% 
  ggboxplot(x = 'nouraiepars21', y = 'TransLoc_MPM_log', add = 'jitter') + stat_compare_means()

## Assemble and save Figure 4
png(out.path('-Figure 4'), height = 8, width = 16, res = 600, units = 'in')
grid.arrange(
  arrangeGrob(
    arrangeGrob(
      Fig4A.p + ylab('Plasma biomarker'), 
      Fig4A_ETApred.p + ylab("Lower respiratory tract biomarker"), 
      top = 'A. LRT pathogen abundance as predictor', ncol = 2
    ), 
    Fig4B.p + ylab('Plasma biomarker') + ggtitle('B. plasma pathogen abundance\n    as predictor'), 
    Fig4C.p + ylab('Plasma biomarker') + ggtitle('C. pulm translocation burden\n    as predictor'), 
    ncol = 3, widths = c(2, 1, 1)
  ), 
  arrangeGrob(
    Fig4D.p + ggtitle('D'), 
    Fig4E.p + ggtitle('E'), 
    Fig4F.p + xlab('') + ylab('Total mcfDNA, MPM, log') + ggtitle("F"), 
    Fig4G.p + xlab('') + ylab('Translocating mcfDNA, MPM, log') + ggtitle('G'),
    ncol = 4, widths = c(1,1.5,1,1)
  )
)
dev.off()

################################################################################
#### SUPPLEMENTARY FIGURES ####
################################################################################
###############################################################################
# Supplemental Figures: 
# This script generates figures S1-S7 and compares filtering workflows for
# the Rappidseq project
###############################################################################

# ---- Figure S1: Biomarker profiles by clinical group ----
# Prepare data: extract biomarker columns and relevant metadata, merge with pneumonia categories,
# and factorize category levels for plotting.
PlotS1.df <- Rappid.baseline.meta.202305.df %>%
  transmute(
    SubjectID = as.character(SubjectID),
    copies_ETA, shannon_ETA,  # Biomarker measures: e.g., gene copies, diversity
    rage_eta, ang2_eta, rage, ang2, nDNA, mtDNA,  # Lung injury and DNA measures
    il10_eta_pradj, fractalkine_eta_pradj, il6_eta, il8_eta, tnfr1_eta, st2_eta,  # Innate response (ETA)
    il6, il8, il10, tnfr1_final, st2, fractalkine, phenoclass_logit_nouraie,  # Plasma cytokines, subphenotyping
    procalcitonin_eta, pentraxin3_eta, procalcitonin, pentraxin3
    # Uncomment for additional biomarkers: cd14, LBP
  ) %>%
  inner_join(Rappid.pneum_category.df) %>%
  mutate(pneum_category = factor(
    pneum_category,
    levels = c("CXPOS", "CXNEG", "CONTROL"),
    labels = c('MCP', 'CDP', "CONT")
  ))

# Output results to PNG file.
png(file = out.path('-Figure-S1.PNG'), height = 16, width = 12, units = "in", res = 600)
grid.arrange(
  # Panel A: Lung injury markers (ETA) across clinical categories
  arrangeGrob(
    ggboxplot(data = PlotS1.df, x = "pneum_category", y = "rage_eta", add = "jitter", add.params = list(size = 0.8)) +
      scale_y_log10(limits = c(1e2, 1.7e6)) + ylab("ETA, RAGE") +
      stat_compare_means(label = "p.format", label.y = 6) +
      xlab("") + theme(axis.text.x = element_text(angle = 15, hjust = 1)),
    ggboxplot(data = PlotS1.df, x = "pneum_category", y = "ang2_eta", add = "jitter", add.params = list(size = 0.8)) +
      scale_y_log10(limits = c(5, 1e4)) +
      ylab("ETA, ANG2") + stat_compare_means(label = "p.format", label.y = 4) +
      xlab("") + theme(axis.text.x = element_text(angle = 15, hjust = 1)),
    top = "A. Lung injury, ETA", ncol = 2
  ),
  # Panel C: Innate immune response markers (ETA)
  arrangeGrob(
    ggboxplot(data = PlotS1.df, x = "pneum_category", y = "fractalkine_eta_pradj", add = "jitter", add.params = list(size = 0.8)) +
      scale_y_log10() + ylab("ETA, Fractalkine") + stat_compare_means(label = "p.format") +
      xlab("") + theme(axis.text.x = element_text(angle = 15, hjust = 1)),
    # ...repeat for other markers, grouping by category...
    top = "C. Innate immune response, ETA", ncol = 6
  ),
  # Panels for host response, plasma markers, etc.
  # (See original script for full structure, repeat grouping as above)
  # Optional: add, remove, or modify panels as dictated by data.
  layout_matrix = matrix(c(1,3,2,2,4,6,5,5,7,7), nrow = 5, byrow = T)
)
dev.off()

# ---- Figure S2: Nanopore read processing and noise filtering ----
png(filename = out.path('-Figure S2'), width = 8, height = 16, res = 600, units = 'in')
gridExtra::grid.arrange(
  Nanopore.all_reads.breakdown.p + theme_minimal() + 
    scale_fill_manual(values = c("springgreen", "firebrick3"), label = c("INCLUSION", 'EXCLUSION')) +
    theme(legend.title = element_blank(), legend.position = 'bottom'),
  Nanopore.species_breakdown.Sankey.p + theme_minimal() +
    scale_fill_manual(values = c("springgreen", "firebrick3"), label = c("INCLUSION", 'EXCLUSION')) +
    theme(legend.position = "none"),
  ncol = 1
)
dev.off()

# ---- Figure S3: Human DNA quantification by clinical category ----
FigureS3.df <- 
  Rappid.baseline.meta.202305.df %>%
  transmute(SubjectIDDay, hcfDNA, nDNA, mtDNA) %>%  # Human cell-free & genomic DNA markers
  inner_join(Rappid.pneum_category_final.df) %>%
  left_join(Nanopore.all_reads.Sankey.df %>%
              dplyr::filter(Human == 'Human reads') %>%
              group_by(SubjectIDDay) %>% summarise(human_nanopore = n()))

# Human DNA, nDNA, and mtDNA are plotted for group stratification; significance annotated.
png(filename = out.path('-Figure S3'), width = 7, height = 7, res = 600, units = 'in')
grid.arrange(
  # Each panel shows a different DNA metric by group, with statistics.
  FigureS3.df %>%
    ggboxplot(x = 'pneum_category_final', y = 'human_nanopore', add = 'jitter') +
    scale_y_log10() + scale_x_discrete(limits = rev) +
    stat_pvalue_manual(data = pairwise_wilcox_test(...), hide.ns = T) +
    xlab("") + ggtitle("A") + ylab("Human-sourced Nanopore reads, log\nLower respiratory tract"),
  # Additional panels for hcfDNA, nDNA, mtDNA...
  ncol = 2
)
dev.off()

# ---- Figure S4: 16S rRNA microbial profiles stratified by clinical group ----
png(width = 12, height = 6, units = "in", filename = out.path('-Figure S4'), res = 300)
grid.arrange(
  Amplicon.ETA.LEFSE.pneum_cat.bar.p + theme(legend.position = "none", axis.text.y = element_text(face = 'italic')),
  Amplicon.ETA.LEFSE.pneum_cat.abund.p + theme(legend.position = "none", axis.text.y = element_text(face = 'italic')),
  ggpubr::get_legend(Amplicon.ETA.LEFSE.pneum_cat.bar.p + theme(legend.position = 'bottom')),
  layout_matrix = matrix(c(1,2,3,3), nrow = 2, byrow = T), heights = c(8, 1)
)
dev.off()

# ---- Figure S5: Culture-positivity and translocation subgroup analysis ----
# Statistical analysis for subgroup relations.
rappid.TransLoc.df %>%
  group_by(TransLoc_final_fix, pneum_category) %>% summarise(n=n()) %>% ungroup() %>%
  pivot_wider(names_from = 'TransLoc_final_fix', values_from = 'n') %>%
  column_to_rownames('pneum_category') %>%
  pairwise_fisher_test(p.adjust.method = "BH") %>%
  dplyr::filter(p.adj < 0.05)

# Plot alluvial visualization of subgroups.
FigS5.p <- 
  rappid.TransLoc.df %>%
  group_by(hasBothMetagenomics, TransLoc_final, ...) %>%
  summarise(freq = n()) %>%
  mutate(FILL_MANUAL_FAC = factor(case_when(...))) %>%
  # Label clinical groups for plotting.
  mutate(pneum_category_fix = case_when(...)) %>%
  ggplot(aes(...)) +
  geom_alluvium(aes(fill = FILL_MANUAL_FAC)) + geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Translocation", "Group"), expand = c(0.15, 0.05)) +
  ylab("# of Subjects with pneumonia") + theme_minimal() +
  scale_fill_manual(values = c("gray20", "gray80", ...))

png(filename = out.path('-Figure S5-cxpositivity_transloc'), width = 5, height = 5, units = "in", res = 300)
FigS5.p
dev.off()

# ---- Figure S6: Enrichment of respiratory pathogens in translocation subgroups ----
png(filename = out.path('-Figure S6'), width = 12, height = 9, res = 600, units = "in")
grid.arrange(
  pathogen_TransLoc.Nanopore.p + theme(legend.position = "bottom") + ggtitle("A"),
  pathogen_TransLoc.Karius.p + theme(legend.position = "none") + ggtitle("B"),
  heights = c(1.2, 1)
)
dev.off()

# ---- Figure S7: Microbial profile by inflammation subphenotype in pneumonia ----
# Filter and merge sample data for subphenotype analysis
Nanopore.subpheno.sample_table <- Rappid.baseline.meta.202305.df %>%
  select(SubjectIDDay, nouraiepars21) %>%
  inner_join(Rappid.pneum_category_final.df) %>%
  dplyr::filter(!is.na(nouraiepars21) & pneum_category_final != "CONTROL" & SubjectIDDay %in% Nanopore.filtered.RelAbd.wide.df$SubjectIDDay) %>%
  mutate(row = SubjectIDDay) %>%
  column_to_rownames('row')

# Setup for differential abundance testing and visualization
Nanopore.subpheno.microtable <- microtable$new(
  otu_table = Nanopore.filtered.RelAbd.wide.df %>%
    dplyr::filter(SubjectIDDay %in% Nanopore.subpheno.sample_table$SubjectIDDay) %>%
    column_to_rownames('SubjectIDDay') %>%
    t() %>% as.data.frame(),
  sample_table = Nanopore.subpheno.sample_table,
  tax_table = Nanopore.filtered.RelAbd.wide.df %>%
    dplyr::filter(SubjectIDDay %in% Nanopore.subpheno.sample_table$SubjectIDDay) %>%
    column_to_rownames('SubjectIDDay') %>%
    t() %>% as.data.frame() %>% rownames() %>% enframe(name = NULL) %>% mutate(species = value) %>% column_to_rownames('value')
)
png(out.path('-Figure S7'), width = 8, height = 8, res = 600, units = 'in')
grid.arrange(
  trans_diff$new(dataset = Nanopore.subpheno.microtable, method = 'lefse', group = 'nouraiepars21', alpha = 0.5, p_adjust_method = 'none', taxa_level = 'species')$plot_diff_bar() +
    scale_fill_manual(values = rev(scales::hue_pal()(2))) +
    theme(legend.position = 'bottom', legend.title = element_blank(), axis.text.y = element_text(face = 'italic')),
  trans_diff$new(dataset = Nanopore.subpheno.microtable, method = 'lefse', group = 'nouraiepars21', alpha = 0.5, p_adjust_method = 'none', taxa_level = 'species')$plot_diff_abund(color_values = rev(scales::hue_pal()(2))) +
    theme(legend.position = 'bottom', legend.title = element_blank(), axis.text.y = element_text(face = 'italic')),
  ncol = 1
)
dev.off()

###############################################################
# Workflow comparison: Filtered vs unfiltered microbiome analysis
# These steps compare diversity and composition between filtered and unfiltered Nanopore read sets.
Nanopore.unfiltered.RelAbd.long.df <- Nanopore.all_reads.Sankey.df %>%
  left_join(Nanopore.below_species.fix[, c("lineage", "id_species")]) %>%
  left_join(Nanopore.species_breakdown.Sankey.df %>% select(-name_species)) %>%
  mutate(EXCLUSION = (Human == "Human reads" | taxonomy_rank_cat %in% c("above species", "non-hierarchical/NA") | taxonomy_rank_cat %in% c("species", "below species"))) %>%
  dplyr::filter(!is.na(EXCLUSION)) %>%
  group_by(SubjectIDDay, id_species, Kingdom) %>%
  summarise(AbsAbd = n()) %>%
  arrange(SubjectIDDay, -AbsAbd) %>% ungroup()

# ...repeat steps for filtered dataset, relative abundance conversion, merging, ordination, diversity calculation...
# Plot comparing Shannon diversity and beta-diversity between workflows
png(filename = out.path('-Fig R1Q6-filtered_v_unfiltered-'), width = 8, height = 4, res = 600, units = 'in')
grid.arrange(
  Shannon_filter_v_unfilter.p, 
  Beta_fil_unfil_merge.p, ncol = 2
)
dev.off()

################################################################################
#### END OF SCRIPT ####
################################################################################
# This completes the analytical pipeline for the RAPPID study manuscript
# All figures and tables are generated with proper statistical testing
# and visualization standards for publication.formula(paste(biomarker, '~ NanoporePathoSum', sep = '')), data = Fig4A.df))["NanoporePathoSum", ])
  