rm(list=ls())
library(pCIA)
library(sjmisc)
library(tidyverse)
library(readxl)

`%ni%` <- Negate(`%in%`)
#######################################
#metabolics 
#######################################
metabolics <- read.csv("C:/Users/zqr20/Documents/349/TJ.GDM800.HM700.Quant.csv")

metabolics_sampleinfo <- read_excel("C:/Users/zqr20/Documents/349/349对GDM病例对照的样本对应情况及检测情况20221101统计.xlsx",
                                   sheet = "代谢") %>%
  mutate(sample = paste0("X",tobetest2023)) %>%
  dplyr::rename(ID = 孕保手册号) %>%
  mutate(trimester = case_when(孕期 == "孕早期" ~ "early",
                               孕期 == "孕中期" ~ "middle",
                               孕期 == "孕晚期" ~ "late")) %>%
  dplyr::select(ID,sample,group,trimester) #%>%
  #filter(sample != "X17P0014065" & sample != "X17P0028143")



metabolics_matrix  <- metabolics  %>%
  #dplyr::select(-c(X17P0014065,X17P0028143)) %>%
  dplyr::filter(Name != "batch") %>%
  mutate(ID = paste0("M",c(1:413))) %>%
  dplyr::select(-Name) %>%
  column_to_rownames("ID") %>%
  rotate_df() %>%
  rownames_to_column("sample") %>%
  left_join(metabolics_sampleinfo,by ="sample")



save(metabolics_sampleinfo, file = "metabolics_sampleinfo.rda")
save(metabolics_matrix, file = "metabolics_matrix.rda")



#######################################
#######################################
lipids <- read.csv("C:/Users/zqr20/Documents/349/TJ.GDM.Lipids.Quant.csv")

lipids_matrix <- lipids  %>%
  dplyr::select(where(~ any(.x %ni% c("P9")))) %>%
  dplyr::filter(Name != "batch") %>%
  select_if(~sum(!is.na(.)) > 0) %>%
  mutate(ID = paste0("L",c(1:1349))) %>%
  dplyr::select(-Name) %>%
  column_to_rownames("ID") %>%
  rotate_df() %>%
   drop_na() %>%
  rownames_to_column("sample") %>%
  left_join(metabolics_sampleinfo,by ="sample")

lipids_sampleinfo <- metabolics_sampleinfo 
save(lipids_sampleinfo, file = "lipids_sampleinfo.rda")
save(lipids_matrix, file = "lipids_matrix.rda")
#######################################
#microbiome
#######################################

###load sample information 
sample_summary_meta <- read_excel("C:/Users/zqr20/Documents/349/349对GDM病例对照的样本对应情况及检测情况20221101统计.xlsx", 
                                  sheet = "孕妇粪便meta")


load("C:/Users/zqr20/Documents/349/species_filter_1.rda")


###clean sample information
sample_summary_meta$保健手册号 <- as.character(sample_summary_meta$保健手册号)
microbiome_sampleinfo <- sample_summary_meta %>%
  dplyr::distinct(样本编号,.keep_all = T) %>%
  dplyr::rename(status = `20221101数据下机`,
                ID = 保健手册号,
                sample = 样本编号,
                disease =病例对照) %>%
  filter(status=="已交付") %>%
  mutate(trimester = case_when(孕期 == "孕早期" ~ "early",
                               孕期 == "孕中期" ~ "middle",
                               孕期 == "孕晚期" ~ "late")) %>%
  mutate(disease = case_when(disease == "病例组" ~ "GDM",
                             disease == "GDM" ~ "GDM",
                             disease == "正常对照" ~ "Control",
                             disease == "对照组" ~ "Control")) %>%
  dplyr::select(ID,sample,trimester,disease)


microbiome_matrix  <- species_filter_1  %>%
  rotate_df() %>%
  mutate(across(where(is.numeric), ~log10(.+0.00001))) %>%
  rownames_to_column("sample") %>%
  left_join(microbiome_sampleinfo,by ="sample")


save(microbiome_sampleinfo, file = "microbiome_sampleinfo.rda")
save(microbiome_matrix, file = "microbiome_matrix.rda")



#######################################
#cfRNA
#######################################
cfRNA_sampleinfo <- read.csv("C:/Users/zqr20/Documents/349/cfRNA_coldata_1550.txt", sep="") %>%
  mutate(group = GDMgroup2023) %>%
  dplyr::select(ID,sample,trimester,group)  

cfRNA_sampleinfo$ID <- as.character(cfRNA_sampleinfo$ID)



cfRNA_tpm_1550 <- read.csv("C:/Users/zqr20/Documents/349/cfRNA_tpm_1550.csv")

cfRNA_matrix <- cfRNA_tpm_1550  %>%
  column_to_rownames("X") %>%
  mutate(across(everything(),~log10(.+0.00001)))  %>%
  rotate_df() %>%
  rownames_to_column("sample") %>%
  right_join(cfRNA_sampleinfo, by = "sample") 


save(cfRNA_sampleinfo, file = "cfRNA_sampleinfo.rda")
save(cfRNA_matrix, file = "cfRNA_matrix.rda")


