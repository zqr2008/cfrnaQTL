rm(list=ls())
library(tidyverse)
library(sjmisc)
library(R.utils)
library(readxl)
library(reticulate)
library(rgl)
library(knitr)
knit_hooks$set(webgl = hook_webgl)

options(scipen = 200)
`%ni%` <- Negate(`%in%`)
setwd("C:/Users/zqr20/Documents/349")
use_python("C:/Users/zqr20/AppData/Local/Programs/Python/Python311")



##################################################
#load each sampleinfo(3) and each matrix(3) 
##################################################
filedir <- list.files(pattern = "*sampleinfo.rda|*matrix.rda")
for (i in filedir){
  load(i)
}

##################################################
#load threshold results
#store_result for microbiome threshold
#metabolitereult for metabolite threshold
#rnaresult for rna threshold
##################################################
load("C:/Users/zqr20/Documents/349/store_result.rda")
metaboliteresult <- read_excel("~/349/metaboliteresult.xlsx", 
                               sheet = "metabolites_HMBD")
rnaresult  <- read_excel("deg_all.xlsx") %>%
  dplyr::mutate(period = case_when(period == "1"~"T1",
                            period == "2"~"T2",
                            period == "all" ~"T12"))

lipidresult <- read.csv("~/349/GDM.csv") %>%
  dplyr::rename("T12" = padj) %>%
  dplyr::select(lipid, T12)
lipidresult1 <- read.csv("C:/Users/zqr20/Documents/349/first trimester.csv") %>%
  dplyr::rename("T1" = padj)  %>%
  dplyr::select(lipid, T1)
lipidresult2 <- read.csv("C:/Users/zqr20/Documents/349/second trimester.csv") %>%
  dplyr::rename("T2" = padj) %>%
  dplyr::select(lipid, T2)


lipidresult <- lipidresult %>%
  left_join(lipidresult1, by = "lipid") %>%
  left_join(lipidresult2, by = "lipid")
##################################################
#load preassigned function
#filtermetabolics function: 
#filtermetabolics(trimester = c("T1","T2","T12",
#                               threshold.score = ,
#                               threshold.prna = ,
#                               threshold.pmicrobiome =)
#
##################################################
filtermetabolics <- function(trimester=NULL,
                             threshold.score = NULL,
                             threshold.prna = NULL,
                             threshold.pmicrobiome = NULL,
                             threshold.lipid = NULL){
  
  
  metaboliteresult <- metaboliteresult %>%
  dplyr::filter( get(trimester) >= as.numeric(threshold.score)) %>%
  dplyr::select(ID) %>%
  distinct(ID) 

  filterlist <- as.vector(metaboliteresult$ID)
  filterlist <- append(filterlist,c("ID","trimester"))
  filteredmetabolics <- metabolics_matrix[,(colnames(metabolics_matrix)%in% filterlist)]
  
  rnaresult <- rnaresult %>%
    filter(period == trimester) %>%
    filter(P.Value <= threshold.prna) %>%
    select(SYMBOL) %>%
    dplyr::rename(rna = SYMBOL)
  
  filterlist <- append(filterlist,as.vector(rnaresult$rna))
  filteredrna<- cfRNA_matrix[,(colnames(cfRNA_matrix) %in% filterlist)]
  
  microbioeresult <- store_result %>%
    dplyr::filter(get(trimester) <= as.numeric(threshold.pmicrobiome)) %>%
    dplyr::select(microbiomename) %>%
    distinct(microbiomename) 

  
  filterlist <- append(filterlist,as.vector(microbioeresult$microbiomename))
  filteredmicrobiome <- microbiome_matrix[,(colnames(microbiome_matrix)%in% filterlist)]
  
  
  lipidresult <- lipidresult %>%
    dplyr::filter(get(trimester) <= as.numeric(threshold.lipid)) %>%
    dplyr::select(lipid) %>%
    distinct(lipid) 
  
  filterlist <- append(filterlist,as.vector(lipidresult$lipid))
  filteredlipid <- lipids_matrix[,(colnames(lipids_matrix)%in% filterlist)]
  
  
  dataframelist <- list(filteredmetabolics,filteredmicrobiome,filteredrna,filteredlipid)
  
  return(dataframelist)
}

rownamematch <- function(x){
  newname <- x  %>%
    dplyr::mutate(newid = paste0(ID,trimester)) %>%
    column_to_rownames("newid") %>%
    select(-c(ID,trimester))

  return(newname)
}

getlabel <- function(datainput = NULL,
                     matchedsample = NULL,
                     option = NULL){
  
  sampletrimester <- as.data.frame(datainput)  %>%
    rownames_to_column("newid") %>%
    dplyr::mutate(ID = substr(newid,start = 1,stop= 14)) %>%
    dplyr::mutate(group = substr(newid,start = 15,stop= 30)) %>%
    left_join(matchedsample, by = "ID") %>%
    distinct(newid,.keep_all = T) %>%
    dplyr::mutate(newlabel = paste0(group,disease))
  
  if (option == "case_control"){
    label <- sampletrimester %>%
      dplyr::select(disease,newid) %>%
      column_to_rownames("newid")
  }
  if (option == "trimester"){
    label <- sampletrimester %>%
      dplyr::select(trimester,newid)  %>%
      column_to_rownames("newid")
  } 
  
  if (option == "both"){
    label <- sampletrimester %>%
      dplyr::select(newlabel,newid) %>%
      column_to_rownames("newid")
  }
  
  label <- label %>%
    dplyr::mutate_at(.vars = vars(1),.funs = as.factor) %>%
    dplyr::mutate_at(.vars = vars(1),.funs = as.numeric)
  return(label)
}

dataframelist <- list()

#trimestervector <- c("T1","T2","T12")
rnavector <- rev(seq(min(rnaresult$adj.P.Val),0.5, by = 0.04))
#scorevector <- seq(0.4,2.4, by = 0.3)
#microbiomevector <- rev(seq(0.002,
#                            max(as.numeric(store_result$T12)), by = 0.1))


trimestervector <- c("T1","T2","T12")
rnavector <- rev(seq(0.04,0.2, by = 0.05))
scorevector <- seq(0.4,2.4, by = 0.3)
microbiomevector <- rev(seq(0.9,1, by = 0.1))
lipidvector <-seq(0.4,2.4, by = 0.3)


trimestervector <- c("T1")
rnavector <- c(0.2)
scorevector <- c(1.4)
microbiomevector <- c(1)
lipidvector <-c(1.4)

trialn = length(trimestervector)*length(rnavector)*length(scorevector)*length(microbiomevector)*length(lipidvector)

biglist <- list()
biglistname <- list()
big1 <- list()
big2 <- list()
big3 <- list()
big4 <- list()
big5 <- list()



counter = 1 
for (x in trimestervector){
 for (j in rnavector){
   for (k in scorevector){
     for (z in microbiomevector){
       for (c in lipidvector){
       
     data <- filtermetabolics(trimester = x,
                         threshold.prna = j,
                         threshold.score = k,
                         threshold.pmicrobiome = z,
                         threshold.lipid =  c)
     
     coverteddata <- list()
     
     for (i in (1:length(data))){
       temp <- rownamematch(data[[i]])
       coverteddata[[i]] <- temp
       }
     
     matcheddata <- list()
     index <- c(1:length(coverteddata))

     for (i in (1:length(coverteddata))){
  
  temp <- coverteddata[[i]]
  leftindex <- index[index %ni% i]
  index1 <- leftindex[1]
  index2 <- leftindex[2]
  index3 <- leftindex[3]
  
  sub <- temp[(rownames(temp) %in% rownames(coverteddata[[index1]])) &
               (rownames(temp) %in% rownames(coverteddata[[index2]])) &
                (rownames(temp) %in% rownames(coverteddata[[index3]])),]

  sub <- sub[order(rownames(sub)), ]

  matcheddata[[i]] <- sub %>%
    dplyr::mutate(across(everything(),~as.numeric(.))) 
}


label <-  as.data.frame(getlabel(datainput = matcheddata[[1]],
                        matchedsample =microbiome_sampleinfo,
                        option = "case_control"))

matcheddata[[5]] <- label
nameslist <- list("subX","subY","subZ","subC","label")




########################################################
#only for SNF py
########################################################
for (i in (1:length(matcheddata))){
  write.table(as.data.frame(matcheddata[i]),
              file = paste0("C:/Users/zqr20/AppData/Local/Programs/Python/Python311/Lib/site-packages/snfpy-0.2.2-py3.11.egg/snf/tests/data/multi/",nameslist[i],".csv"),
              sep =",",
              fileEncoding="GBK",
              row.names = F,
              col.names = F)
}




py_run_file("C:/Users/zqr20/Documents/349/getoptimalSNF.py") 
biglist <- append(biglist,py$optima)
biglistname <- append(biglistname,paste0("trimestervector",x,
                                         "rnavector",j,
                                         "scorevector",k,
                                         "microbiomevector",z))

big1 <- append(big1,x)
big2 <- append(big2,j)
big3 <- append(big3,k)
big4 <- append(big4,z)
big5 <- append(big5,c)

print(paste0("performed ",counter," times iteration"))
counter = counter +1

if(py$optima>=0.1){
  break
  }
     }
   }
 }
 }
}



########################################################
#plot 3d for thresholds
########################################################


visual3d <- cbind(t(as.data.frame(big1)),t(as.data.frame(big2)),
            t(as.data.frame(big3)),t(as.data.frame(big4)),
            t(as.data.frame(big5)),t(as.data.frame(biglist)))
visual3d <- as.data.frame(visual3d) %>%
  mutate_at(.vars = vars(2:5),.funs = as.numeric) %>%
  dplyr::rename(RNA = V2,
                metablite = V3,
                microbiome = V4,
                lipid = V5,
                SNFvalue = V6)

mycolors <- c('royalblue1', 'darkcyan', 'oldlace')
visual3dcolor<−mycolors[as.numeric(as.factor(visual3dcolor <- mycolors[as.numeric(as.factor(visual3dV1))]



par(mar=c(0,0,0,0))
plot3d( 
  x=visual3dmicrobiome,y=visual3dmicrobiome, y=visual3dmetablite, z=visual3d$SNFvalue, 
  col = visual3d$color, 
  type = 's', 
  radius = .03,
  xlab="microbiome", ylab="metablite", zlab="SNFvalue")



########################################################
#only for mogonet
########################################################
relabel <- matcheddata[[5]]  %>%
  mutate(disease = case_when(disease == "1"~ 0,
                             disease == "2"~ 1))
matcheddata[[5]] <- relabel


for (i in (1:length(matcheddata))){
  write.table(as.data.frame(matcheddata[i]),
              file = paste0("C:/Users/zqr20/Documents/349/MOGONET/gdmstudy/multi/",nameslist[i],".csv"),
              sep =",",
              fileEncoding="GBK",
              row.names = F,
              col.names = F,)
}


for (i in (1:length(matcheddata))){
  write.table(as.data.frame(matcheddata[i]),
              file = paste0("C:/Users/zqr20/Documents/349/MOGONET/gdmstudy/labelled/",nameslist[i],".csv"),
              sep =",",
              fileEncoding="GBK",
              row.names = T,
              col.names = T)
}
