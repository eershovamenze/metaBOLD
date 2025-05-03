
library(dplyr)
library(tidyr)
library(readxl)
library(stringr)



# MOTHUR taxonomy ---------------------------------------------------------

tax_ranks <- read.csv('~/Documents/MetaBOLD_intercallibration/mothur_ranks.csv')
tax_table <- read.delim('~/Documents/MetaBOLD_intercallibration/data/reference/OTU_reference.MZGmothur_coi__MZGdbALL__o00__A.wang.txt', header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(tax_table) <- c("MIN", "taxonomy")
tax_table$taxonomy <- sub(";$", "", tax_table$taxonomy)
tax_table1 <- tax_table %>% 
  separate_wider_delim(taxonomy, ";", names=tax_ranks$rank) %>% 
  melt(id.vars='MIN', value.name='taxa', variable.name='rank') %>% 
  separate_wider_delim(taxa, delim="(", names=c('Taxa', 'Confidence'), too_few='align_start')%>% 
  mutate(Confidence = replace_na(Confidence, '0)'))%>% 
  mutate(Confidence = as.numeric(sub("\\)", "", Confidence)))%>% 
  filter(Confidence >= 95 & !str_detect(Taxa, "unclassified"))%>% 
  merge(tax_ranks, by='rank')%>% 
  group_by(MIN)%>% 
  mutate(Lowest_level = max(final_level))%>% 
  filter(Lowest_level==final_level)%>% 
  mutate(Taxa = gsub("_EXT", "", Taxa))%>%
  mutate(Taxa = gsub("_", " ", Taxa))%>%
  merge(tax_table, all.y=T)%>%
  mutate(Taxa = if_else(final_level < 13, "unclassified", Taxa))
Species_to_check_MZGdb <- tax_table1%>%
  filter(rank_simple == 'species')%>%
  distinct(Taxa)
#write.csv(as.data.frame(Species_to_check_MZGdb), 'WORMS_to_check_MZGdb.csv')
MZGdb_WORMS <- read.csv('WORMS_MZGdb_matched.csv')
mzgdb_taxonomy_final <- tax_table1 %>%
  left_join(MZGdb_WORMS, by = "Taxa") %>%
  mutate(Taxa = coalesce(ScientificName_accepted, Taxa)) %>%
  dplyr::select(
    MIN,
    MZGdb_final_assignment = Taxa,
    MZGdb_rank             = rank_simple,
    MZGdb_confidence       = Confidence) %>%
  mutate(MZGdb_final_assignment = replace_na(MZGdb_final_assignment, 'unclassified'),
         MZGdb_rank = replace_na(MZGdb_rank, 'no rank'),
         MZGdb_confidence = replace_na(MZGdb_confidence, 0))  %>%
  mutate(MIN = sub(" ", "", MIN))

write.csv(mzgdb_taxonomy_final, "mothur_taxonomy_converted.csv", row.names = FALSE)

# mkLTG taxonomy ----------------------------------------------------------
mkltg_tax <- read.delim('~/Documents/MetaBOLD_intercallibration/data/reference/Taxonomy_MIN_global_ltg.tsv', header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ltg_ranks <- read.csv('~/Documents/MetaBOLD_intercallibration/ltg_ranks.csv')
unique(mkltg_tax$ltg_rank)
mkltg_tax <- mkltg_tax %>% mutate(ltg_rank = na_if(ltg_rank, "")) %>%
  mutate(ltg_rank = replace_na(ltg_rank, 'no rank'))%>% 
  merge(ltg_ranks, by='ltg_rank')%>% 
  mutate(ltg_name = gsub(',', '', ltg_name)) 
  
mkltg_tax$Final_assignment <- mkltg_tax$ltg_name
mkltg_tax[mkltg_tax$final_level<13,]$Final_assignment <- 'unclassified'
#WORMS check
Species_to_check_mkltg <- mkltg_tax%>%
  filter(ltg_rank == 'species')%>%
  distinct(Final_assignment)
write.csv(as.data.frame(Species_to_check_mkltg), 'WORMS_to_check_ltg.csv')
WORMS_ltg <- read.csv('WORMS_ltg_matched.csv')
str(mkltg_tax)
mkltg_taxonomy_final <- mkltg_tax %>%
  left_join(WORMS_ltg, by = "Final_assignment") %>%
  mutate(ScientificName_accepted = coalesce(ScientificName_accepted, Final_assignment)) %>%
  dplyr::select(MIN = seqid,
        mkLTG_final_assignment = ScientificName_accepted,
         mkLTG_superkingdom = superkingdom,
         mkLTG_kingdom = kingdom,
         mkLTG_phylum = phylum,
         mkLTG_order = order,
         mkLTG_assigned = ltg_name,
        mkLTG_rank =ltg_rank, 
         mkLTG_confidence = pid) 
write.csv(mkltg_taxonomy_final, "mkltg_taxonomy_converted.csv", row.names = FALSE)

# BOLDigger taxonomy ----------------------------------------------------------
bold_tax <- read_excel('~/Documents/MetaBOLD_intercallibration/data/reference/OTU_reference_identification_result.xlsx')
bold_hits <- read_excel('~/Documents/MetaBOLD_intercallibration/data/reference/OTU_reference_bold_results_part_1.xlsx')
str(bold_tax)
bold_tax <- bold_tax %>% 
  mutate(selected_level = replace_na(selected_level, 'no rank'))
bold_tax$Final_assignment <- apply(bold_tax, 1, function(row) {
  sel <- row["selected_level"]
  if (sel %in% colnames(bold_tax)) {
    return(row[[sel]])
  } else {
    return("unclassified")
  }
})


#Dealing with flagged taxa
bold_tax_flagged <- bold_tax %>%
  filter(str_detect(flags, "\\b2\\b"))
#Species double hits
bold_hits_flagged_species <- bold_hits %>%
  filter(id %in% bold_tax_flagged[bold_tax_flagged$selected_level=='Species',]$id)%>%
  filter(pct_identity>97)%>%
  filter(!is.na(Species)) %>% 
  group_by(id, Order, Family, Genus, Species)%>%
  summarise(hits = length(database), pct_identity=mean(pct_identity))%>%
  group_by(id)%>%
  mutate(total_hits=sum(hits))%>%
  mutate(percent_hits=hits/total_hits)%>%
  group_by(id, total_hits) %>%
  summarise(best_id = Species[which.max(percent_hits)], percent_hits=max(percent_hits))
bold_hits_flagged_species$Final_assignment <- bold_hits_flagged_species$best_id
bold_hits_flagged_species[bold_hits_flagged_species$percent_hits < 0.8 | bold_hits_flagged_species$total_hits < 10,]$Final_assignment <- "_Manual_check"
bold_hits_flagged_species[bold_hits_flagged_species$percent_hits <= 0.5 & bold_hits_flagged_species$Final_assignment == "_Manual_check",]$Final_assignment <- "genus"
bold_hits_flagged_species$Level_check <- 'Species'
bold_hits_flagged_species[bold_hits_flagged_species$Final_assignment == "_Manual_check",]$Level_check <- '_Manual_check'
genus_ids <- c(bold_tax_flagged[bold_tax_flagged$selected_level=='Genus',]$id, bold_hits_flagged_species[bold_hits_flagged_species$Final_assignment=='genus',]$id)
bold_hits_flagged_species <- filter(bold_hits_flagged_species, Final_assignment != 'genus')

#Genus double hits
bold_hits_flagged_genus <- bold_hits %>%
  filter(id %in% genus_ids)%>%
  filter(pct_identity>95)%>%
  filter(!is.na(Genus)) %>% 
  group_by(id, Genus)%>%
  summarise(hits = length(database), pct_identity=mean(pct_identity))%>%
  group_by(id)%>%
  mutate(total_hits=sum(hits))%>%
  mutate(percent_hits=hits/total_hits)%>%
  group_by(id, total_hits) %>%
  summarise(best_id = Genus[which.max(percent_hits)], percent_hits=max(percent_hits))
bold_hits_flagged_genus$Final_assignment <- bold_hits_flagged_genus$best_id
bold_hits_flagged_genus[bold_hits_flagged_genus$percent_hits < 0.8 | bold_hits_flagged_genus$total_hits < 10,]$Final_assignment <- "_Manual_check"
bold_hits_flagged_genus[bold_hits_flagged_genus$percent_hits <= 0.5 & bold_hits_flagged_genus$Final_assignment == "_Manual_check",]$Final_assignment <- "family"
bold_hits_flagged_genus[bold_hits_flagged_genus$percent_hits == 1 | bold_hits_flagged_genus$total_hits == 1,]$Final_assignment <- bold_hits_flagged_genus[bold_hits_flagged_genus$percent_hits == 1,]$best_id
bold_hits_flagged_genus$Level_check <- 'Genus'
bold_hits_flagged_genus[bold_hits_flagged_genus$Final_assignment == "_Manual_check",]$Level_check <- '_Manual_check'
family_ids <- c(bold_tax_flagged[bold_tax_flagged$selected_level=='Family',]$id, bold_hits_flagged_genus[bold_hits_flagged_genus$Final_assignment=='family',]$id)
bold_hits_flagged_genus <- filter(bold_hits_flagged_genus, Final_assignment != 'family')

#Family double hits
bold_hits_flagged_family <- bold_hits %>%
  filter(id %in% family_ids)%>%
  filter(pct_identity>90)%>%
  filter(!is.na(Family)) %>% 
  group_by(id, Family)%>%
  summarise(hits = length(database), pct_identity=mean(pct_identity))%>%
  group_by(id)%>%
  mutate(total_hits=sum(hits))%>%
  mutate(percent_hits=hits/total_hits)%>%
  group_by(id, total_hits) %>%
  summarise(best_id = Family[which.max(percent_hits)], percent_hits=max(percent_hits))
bold_hits_flagged_family$Final_assignment <- bold_hits_flagged_family$best_id
bold_hits_flagged_family[bold_hits_flagged_family$percent_hits < 0.8 | bold_hits_flagged_family$total_hits < 10,]$Final_assignment <- "_Manual_check"
bold_hits_flagged_family[bold_hits_flagged_family$percent_hits <= 0.5 & bold_hits_flagged_family$Final_assignment == "_Manual_check",]$Final_assignment <- "order"
bold_hits_flagged_family[bold_hits_flagged_family$percent_hits == 1,]$Final_assignment <- bold_hits_flagged_family[bold_hits_flagged_family$percent_hits == 1,]$best_id
bold_hits_flagged_family$Level_check <- 'Family'
bold_hits_flagged_family[bold_hits_flagged_family$Final_assignment == "_Manual_check",]$Level_check <- '_Manual_check'
order_ids <- c(bold_tax_flagged[bold_tax_flagged$selected_level=='Order',]$id, bold_hits_flagged_family[bold_hits_flagged_family$Final_assignment=='order',]$id)
bold_hits_flagged_family <- filter(bold_hits_flagged_family, Final_assignment != 'order')

#Order double hits
bold_hits_flagged_order <- bold_hits %>%
  filter(id %in% order_ids)%>%
  filter(pct_identity>85)%>%
  filter(!is.na(Order)) %>% 
  group_by(id, Order)%>%
  summarise(hits = length(database), pct_identity=mean(pct_identity))%>%
  group_by(id)%>%
  mutate(total_hits=sum(hits))%>%
  mutate(percent_hits=hits/total_hits)%>%
  group_by(id, total_hits) %>%
  summarise(best_id = Order[which.max(percent_hits)], percent_hits=max(percent_hits))
bold_hits_flagged_order$Final_assignment <- bold_hits_flagged_order$best_id
bold_hits_flagged_order[bold_hits_flagged_order$percent_hits < 0.8,]$Final_assignment <- "unclassified"
bold_hits_flagged_order$Level_check <- 'Order'
bold_hits_flagged_order[bold_hits_flagged_order$Final_assignment == "unclassified",]$Level_check <- 'no rank'

flagged_all <- rbind(bold_hits_flagged_species, bold_hits_flagged_genus, bold_hits_flagged_family, bold_hits_flagged_order)
flagged_all$Manual_check <- 'automatic'
flagged_all[flagged_all$Final_assignment=='_Manual_check',]$Manual_check <- 'manual'
flagged_all <- flagged_all %>% dplyr::rename(Final_assignment_check = Final_assignment)

bold_tax_manual <- merge(bold_tax, flagged_all_subset[c('id', 'Final_assignment_check', 'Manual_check', 'Level_check')], by='id')
manual_checks <- bold_hits[bold_hits$id %in% bold_tax_manual[bold_tax_manual$Manual_check == 'manual',]$id,]
write.csv(bold_tax_manual, 'bold_tax_final_to_check.csv')
write.csv(manual_checks, 'bold_hits_for_manual_check.csv')

bd_taxonomy_final_checked <- bd_taxonomy_final  %>%
  filter(BD_flag_check != 'none')
  
bold_tax_manual <- merge(bold_tax, flagged_all[c('id', 'Final_assignment_check', 'Manual_check', 'Level_check')], by='id', all.x=T) %>%
  mutate(Manual_check = replace_na(Manual_check, 'none'))
manual_checks <- bold_hits[bold_hits$id %in% bold_tax_final[bold_tax_final$Manual_check == 'manual',]$id,]
write.csv(bold_tax_manual, 'bold_tax_final_to_check.csv')
write.csv(manual_checks, 'bold_hits_for_manual_check.csv')
bold_final_checked <- read.csv('bold_final_checked_upd.csv')

str(bold_final_checked)
Species_to_check <- bold_final_checked%>%
  filter(Level_check_BD == 'Species')%>%
  distinct(Assignment_checked_BD)
write.csv(as.data.frame(Species_to_check), 'Species_to_check.csv')
worms_check_bd <- read.csv('WORMS_matched_BD.csv')
colnames(worms_check_bd)
bd_taxonomy_final <- bold_final_checked %>%
  left_join(worms_check_bd, by = "Assignment_checked_BD") %>%
  mutate(ScientificName_accepted = coalesce(ScientificName_accepted, Assignment_checked_BD)) %>%
  mutate(Level_check_BD = tolower(Level_check_BD)) %>%
  dplyr::select(MIN = id,
                BD_final_assignment = ScientificName_accepted,
                BD_phylum = Phylum,
                BD_class = Class,
                BD_Order = Order,
                BD_Family = Family,
                BD_assigned_raw = Assignment_BD,    
                BD_flags = flags,
                BD_level = Level_check_BD,
                BD_flag_check = Manual_check_BD,
                BD_manual_check_outcome = Manual_outcome_BD,
                BD_reason_outcome = Reason_BD)


# Merged taxonomy ---------------------------------------------------------
lulu_results <- read.csv('LULU_results.csv')
colnames(merged_taxonomy)
merged_taxonomy$BD_flag_check
merged_taxonomy <- bd_taxonomy_final %>%
  left_join(mzgdb_taxonomy_final, by = "MIN") %>%
  left_join(mkltg_taxonomy_final, by = "MIN")%>%
  left_join(lulu_results, by = "MIN")%>%
  mutate(Agreement = ifelse(BD_final_assignment==MZGdb_final_assignment & BD_final_assignment==mkLTG_final_assignment & mkLTG_final_assignment==MZGdb_final_assignment, 'All', 'None'))%>%
  mutate(Agreement = ifelse(Agreement == 'None' & BD_final_assignment==MZGdb_final_assignment, 'BD+MZGbd', Agreement))%>%
  mutate(Agreement = ifelse(Agreement == 'None' & mkLTG_final_assignment==MZGdb_final_assignment, 'MZGbd+mkLTG', Agreement))%>%
  mutate(Agreement = ifelse(Agreement == 'None' & mkLTG_final_assignment==BD_final_assignment, 'BD+mkLTG', Agreement))%>%
  #If all three agree, assign to consensus
  mutate(Merged_final_assignment = ifelse(Agreement=='All', BD_final_assignment, '_check'))%>%
  #If two agree at species level or the third is 'unclassified', assign to consensus. If the consensus is "unclassified", then assign to third
  mutate(Merged_final_assignment = ifelse(Agreement=='BD+MZGbd' & mkLTG_final_assignment == 'unclassified', BD_final_assignment, Merged_final_assignment))%>%
  mutate(Merged_final_assignment = ifelse(Agreement=='BD+MZGbd' & BD_level == 'species', BD_final_assignment, Merged_final_assignment))%>%
  mutate(Merged_final_assignment = ifelse(Agreement=='BD+MZGbd' & BD_final_assignment == 'unclassified', mkLTG_final_assignment, Merged_final_assignment))%>%
  mutate(Merged_final_assignment = ifelse(Agreement=='MZGbd+mkLTG' & BD_final_assignment == 'unclassified', mkLTG_final_assignment, Merged_final_assignment))%>%
  mutate(Merged_final_assignment = ifelse(Agreement=='MZGbd+mkLTG' & MZGdb_rank == 'species', mkLTG_final_assignment, Merged_final_assignment))%>%
  mutate(Merged_final_assignment = ifelse(Agreement=='MZGbd+mkLTG' & mkLTG_final_assignment == 'unclassified', BD_final_assignment, Merged_final_assignment))%>%
  mutate(Merged_final_assignment = ifelse(Agreement=='BD+mkLTG' & MZGdb_final_assignment == 'unclassified', mkLTG_final_assignment, Merged_final_assignment))%>%
  mutate(Merged_final_assignment = ifelse(Agreement=='BD+mkLTG' & BD_level == 'species', mkLTG_final_assignment, Merged_final_assignment))%>%
  mutate(Merged_final_assignment = ifelse(Agreement=='BD+mkLTG' & mkLTG_final_assignment == 'unclassified', MZGdb_final_assignment, Merged_final_assignment))%>%
  #If ambiguity was manually validated with BOLDigger, assign to BOLDigger
  mutate(Merged_final_assignment = ifelse(BD_flag_check == "manual", BD_final_assignment, Merged_final_assignment))%>%
  #If there is no remaining agreement and ambiguity was automatically validated with BOLDigger, assign to BOLDigger
  mutate(Merged_final_assignment = ifelse(Merged_final_assignment=='_check' & BD_flag_check == "automatic", BD_final_assignment, Merged_final_assignment))%>%
  #If there is no remaining agreement and likely pseudogene, then assign 'unclassified'
  mutate(Merged_final_assignment = ifelse(Merged_final_assignment=='_check' & LULU_results == "likely_pseudogene", 'unclassified', Merged_final_assignment))%>%
  mutate(Merged_rank = ifelse(Merged_final_assignment==mkLTG_final_assignment, mkLTG_rank, '_check'))%>%
  mutate(Merged_rank = ifelse(Merged_final_assignment==MZGdb_final_assignment, MZGdb_rank, Merged_rank))%>%
  mutate(Merged_rank = ifelse(Merged_final_assignment==BD_final_assignment, BD_level, Merged_rank))
  
  
#Rest checked manually



  
  

write.csv(merged_taxonomy, 'Merged_taxonomy_final.csv')


