library(easyPubMed)
library(tidyverse)
library(data.table)
library(ggplot2)

#search performed April 13, 2020
pubmed_query_string="(NGS OR \"next generation sequencing\" OR RNA-Seq OR \"mRNA sequencing\" OR \"whole exome sequencing\" OR \"whole-exome sequencing\" OR \"high throughput sequencing\" OR \"high-throughput sequencing\" OR \"DNA sequencing\" OR \"RNA sequencing\" OR WXS OR WGS OR \"whole-genome sequencing\" OR \"whole genome sequencing\") AND (rheumatology OR \"rheumatologic disease\" OR \"rheumatologic disease\" ))"
#814 entries

pubmed_ids <- get_pubmed_ids(pubmed_query_string)
abstracts_xml <- fetch_pubmed_data(pubmed_ids,retmax = 5000)
abstracts_list <- articles_to_list(abstracts_xml)
final=c("")
for (key in abstracts_list){
keywords = custom_grep(key,"KeywordList","char")

if(!is.null(keywords)){
out=str_replace_all(keywords,"<Keyword MajorTopicYN=\"N\">","")
out=str_replace_all(out,"<Keyword MajorTopicYN=\"Y\">","")
out=strsplit(out,"</Keyword>")[[1]]
final=c(final,out)
}
}
final=final[final!=" "]
final=tolower(final)

finalICD11_filtered=c("")

icd11 <- read.delim("~/icd11.txt", stringsAsFactors=FALSE)
icd11$Title=tolower(icd11$Title)

for (key in final){
retrows=icd11[icd11$Title %like% key,]
if (length(retrows) >= 1){
  finalICD11_filtered=c(finalICD11_filtered,key)}}
finalTable_freq=table(finalICD11_filtered)  
View(finalTable_freq)  
#253 entries in finalTable_freq


#search performed April,13
pubmed_query_string="(NGS OR \"next generation sequencing\" OR RNA-Seq OR \"mRNA sequencing\" OR \"whole exome sequencing\" OR \"whole-exome sequencing\" OR \"high throughput sequencing\" OR \"high-throughput sequencing\" OR WXS OR WGS OR \"whole-genome sequencing\" OR \"whole genome sequencing\") AND (\"Deficiency of adenosine deaminase 2\" OR \"chilblain lupus\" OR \"giant cell arteritis\" OR \"neuropsychiatric lupus\" OR \"lupus nephritis\" OR \"Immunoglobulin A nephropathy\" OR \"sarcoidosis\" OR \"rheumatoid arthritis\" OR \"systemic sclerosis\" OR \"systemic lupus erythematosus\" OR \"juvenile idiopathic arthritis\" OR \"myositis\" OR \"osteoarthritis\" OR \"spondyloarthritis\" OR \"dermatomyositis\" OR \"familial mediterranean fever\" OR \"gout\" OR \"uveitis\" OR \"vasculitis\" OR \"granulomatosis with polyangiitis\" OR \"psoriatic arthritis\" OR \"sjögren syndrome\" OR \"sjögren's syndrome\")"

#613 hits
pubmed_ids <- get_pubmed_ids(pubmed_query_string)
abstracts_xml <- fetch_pubmed_data(pubmed_ids,retmax = 5000)
abstracts_list <- articles_to_list(abstracts_xml)
diseases=c("sacroiliitis","Deficiency of adenosine deaminase 2","chilblain lupus","giant cell arteritis","neuropsychiatric lupus","lupus nephritis","Immunoglobulin A nephropathy","sarcoidosis","rheumatoid arthritis","systemic sclerosis","systemic lupus erythematosus","juvenile idiopathic arthritis","myositis","osteoarthritis","spondyloarthritis","dermatomyositis","familial mediterranean fever","gout","vasculitis","granulomatosis with polyangiitis","psoriatic arthritis","sjögren syndrome","sjögren's syndrome")
RNA=c("rna sequencing", "transcriptome sequencing", "rna-seq","rna-based next-generation sequencing","mrna profiles")
WXS=c("whole exome sequencing","whole-exome sequencing","whole-exome-sequencing","amplicon sequencing")
WGS=c("whole-genome sequencing","whole genome sequencing","whole-genome shotgun sequencing")
Bacteria=c("16s","metagenomics")
Epigenomic=c("atac-seq","chip-seq","bisulfite sequencing","wgbs")


result=NULL
for (key in abstracts_list){
review=custom_grep(key,"PublicationType","char")
journal=str_replace(custom_grep(custom_grep(key,"Journal","char"),"Title","char"),"&amp;","&")
year=custom_grep(custom_grep(key,"PubDate","char"),"Year","char")
title=custom_grep(key,"ArticleTitle","char")
title=tolower(title)
abstract = custom_grep(key,"Abstract","char")
abstract=tolower(abstract)
keywords=custom_grep(key,"KeywordList","char")
keywords=tolower(keywords)
title_abstract_key=paste(title,abstract,sep=",")
title_abstract_key=paste(title_abstract_key,keywords,sep=",")

pubmedid=custom_grep(key,"PMID","char")

if(!is.null(title_abstract_key) & length(grep("Review",review))==0){
  assay=""  
    
    for (epi in Epigenomic){
      if (length(grep(epi,title_abstract_key))>0){assay=paste(assay,"Epigenomics",sep=",")}
    }
    for (bac in Bacteria){
      if (length(grep(bac,title_abstract_key))>0){assay=paste(assay,"Metagenomics",sep=",")}
    }
    
    for (rna in RNA){
      if (length(grep(rna,title_abstract_key))>0){
        
        if(!grepl("RNA",assay)){
          assay=paste(assay,"RNA-Seq",sep=",")}}
    }
    for (wxs in WXS){
      
      if (length(grep(wxs,title_abstract_key))>0){
        if(!grepl("WES",assay)){
        assay=paste(assay,"WES",sep=",")}}
    }
    for (wgs in WGS){
      if (length(grep(wgs,title_abstract_key))>0){assay=paste(assay,"WGS",sep=",")}
    }
    for (disease in diseases){
    #print(disease)
      if (length(title_abstract_key)>0){
      if(title_abstract_key %like% disease){
        
      result=rbind(result,c(disease,pubmedid[1],substr(assay,2,nchar(assay)),journal,year))
    }
}
}
}
}


result_df=as.data.frame(result)
names(result_df)=c("disease","pmid","assay","journal","year")
result_df$assay=as.character(result_df$assay)
result_df$disease=as.character(result_df$disease)
result_df$disease=stringr::str_to_title(result_df$disease)

#manually adding missing data by looking at each publication using the PMID
result_df$pmid=as.character(result_df$pmid)
result_df[result_df$pmid=="31428656",]$assay="RNA-Seq"
result_df[result_df$pmid=="31370803",]$assay="Metagenomics"
result_df[result_df$pmid=="31360262",]$assay="RNA-Seq"
result_df[result_df$pmid=="31344123",]$assay="RNA-Seq"
result_df[result_df$pmid=="31337345",]$assay="DNA gene panel"
result_df[result_df$pmid=="31101814",]$assay="WGS,WES"
result_df[result_df$pmid=="30783801",]$assay="DNA gene panel"
result_df[result_df$pmid=="30513227",]$assay="DNA gene panel"
result_df[result_df$pmid=="30319641",]$assay="RNA-Seq"
result_df[result_df$pmid=="29963250",]$assay="RNA-Seq"
result_df[result_df$pmid=="29892977",]$assay="RNA-Seq"
result_df[result_df$pmid=="29891556",]$assay="RNA-Seq"
result_df[result_df$pmid=="29735907",]$assay="DNA gene panel"
result_df[result_df$pmid=="29681619",]$assay="DNA gene panel"
result_df[result_df$pmid=="29659823",]$assay="WES"
result_df[result_df$pmid=="29657145",]$assay="RNA-Seq"
result_df[result_df$pmid=="29432186",]$assay="RNA-Seq"
result_df[result_df$pmid=="28791018",]$assay="DNA gene panel"
result_df[result_df$pmid=="28750028",]$assay="DNA gene panel"
result_df[result_df$pmid=="28728565",]$assay="RNA-Seq"
result_df[result_df$pmid=="28523199",]$assay="DNA gene panel"
result_df[result_df$pmid=="31856934",]$assay="WES"
result_df[result_df$pmid=="31856934",]$disease="DADA2"
result_df[result_df$pmid=="31848804",]$assay="DNA gene panel"
result_df[result_df$pmid=="31848804",]$disease="DADA2"
result_df[result_df$pmid=="32079724",]$assay="RNA-Seq"
result_df[result_df$pmid=="31993940",]$assay="WES"
result_df[result_df$pmid=="31988805",]$assay="RNA-Seq"
result_df[result_df$pmid=="31882654",]$assay="RNA-Seq"
result_df[result_df$pmid=="29681619",]$disease="DADA2"
result_df[result_df$pmid=="31647025",]$assay="RNA-Seq"
result_df[result_df$pmid=="31608065",]$assay="RNA-Seq"
result_df[result_df$pmid=="31598594",]$assay="DNA gene panel"
result_df[result_df$pmid=="31856934",]$year=2019
result_df[result_df$pmid=="29739689",]$year=2018
result_df[result_df$pmid=="31003835",]$year=2019
result_df[result_df$pmid=="31538826",]$year=2019
result_df[result_df$pmid=="31003835",]$assay="Metagenomics"
result_df[result_df$pmid=="23666743",]$assay="WES"
result_df[result_df$pmid=="25285625",]$assay="DNA gene panel"
result_df[result_df$pmid=="25091625",]$assay="Epigenomics,RNA-Seq"
result_df[result_df$pmid=="24795478",]$assay="WGS,RNA-Seq"
result_df[result_df$pmid=="23606709",]$assay="RNA-Seq"
result_df[result_df$pmid=="26245356",]$assay="RNA-Seq"
result_df[result_df$pmid=="32042094",]$assay="DNA gene panel"
result_df[result_df$pmid=="31926670",]$assay="Metagenomics"
result_df[result_df$pmid=="31921732",]$assay="Metagenomics"
result_df[result_df$pmid=="31898522",]$assay="RNA-Seq"
result_df[result_df$pmid=="31880128",]$assay="Metagenomics"
result_df[result_df$pmid=="31874111",]$assay="WGS,WES"
result_df[result_df$pmid=="31779271",]$assay="RNA-Seq"
result_df[result_df$pmid=="31832070",]$assay="DNA gene panel"
result_df[result_df$pmid=="31538826",]$assay="DNA gene panel"
result_df[result_df$pmid=="31523167",]$assay="RNA-Seq"
result_df[result_df$pmid=="31443670",]$assay="DNA gene panel"
result_df[result_df$pmid=="31412876",]$assay="DNA gene panel"
result_df[result_df$pmid=="31384939",]$assay="DNA gene panel"
result_df[result_df$pmid=="31237906",]$assay="RNA-Seq"
result_df[result_df$pmid=="30987788",]$assay="DNA gene panel"
result_df[result_df$pmid=="30850477",]$assay="RNA-Seq"
result_df[result_df$pmid=="30836987",]$assay="RNA-Seq"
result_df[result_df$pmid=="30724444",]$assay="RNA-Seq"
result_df[result_df$pmid=="30700427",]$assay="Metagenomics"
result_df[result_df$pmid=="30544699",]$assay="RNA-Seq"
result_df[result_df$pmid=="30306067",]$assay="RNA-Seq"
result_df[result_df$pmid=="30185675",]$assay="RNA-Seq"
result_df[result_df$pmid=="30123050",]$assay="RNA-Seq"
result_df[result_df$pmid=="29683194",]$assay="DNA gene panel"
result_df[result_df$pmid=="29614084",]$assay="DNA gene panel"
result_df[result_df$pmid=="29490685",]$assay="DNA gene panel"
result_df[result_df$pmid=="29465611",]$assay="DNA gene panel"
result_df[result_df$pmid=="29381936",]$assay="DNA gene panel"
result_df[result_df$pmid=="29371932",]$assay="RNA-Seq"
result_df[result_df$pmid=="29247798",]$assay="RNA-Seq"
result_df[result_df$pmid=="29137139",]$assay="RNA-Seq"
result_df[result_df$pmid=="29129473",]$assay="DNA gene panel"
result_df[result_df$pmid=="29123165",]$assay="RNA-Seq"
result_df[result_df$pmid=="29109423",]$assay="RNA-Seq"
result_df[result_df$pmid=="28808260",]$assay="RNA-Seq"
result_df[result_df$pmid=="28597968",]$assay="WGS,WES"
result_df[result_df$pmid=="30123050",]$assay="RNA-Seq"
result_df[result_df$pmid=="28495399",]$assay="DNA gene panel"
result_df[result_df$pmid=="28348750",]$assay="Metagenomics"
result_df[result_df$pmid=="28339495",]$assay="RNA-Seq"
result_df[result_df$pmid=="28179509",]$assay="Metagenomics"
result_df[result_df$pmid=="28115215",]$assay="WGS"
result_df[result_df$pmid=="27835701",]$assay="RNA-Seq"
result_df[result_df$pmid=="27821747",]$assay="PhIP-Seq"
result_df[result_df$pmid=="27571913",]$assay="Metagenomics"
result_df[result_df$pmid=="27542282",]$assay="RNA-Seq"
result_df[result_df$pmid=="27420474",]$assay="RNA-Seq"
result_df[result_df$pmid=="27390188",]$assay="DNA gene panel"
result_df[result_df$pmid=="27273876",]$assay="RNA-Seq"
result_df[result_df$pmid=="29635721",]$assay="DNA gene panel"
result_df[result_df$pmid=="29200130",]$assay="DNA gene panel"
result_df[result_df$pmid=="28043923",]$assay="WES"
result_df[result_df$pmid=="27183593",]$assay="RNA-Seq"
result_df[result_df$pmid=="26815131",]$assay="DNA gene panel"
result_df[result_df$pmid=="26713667",]$assay="RNA-Seq"
result_df[result_df$pmid=="26603474",]$assay="DNA gene panel"
result_df[result_df$pmid=="26227771",]$assay="DNA gene panel"
result_df[result_df$pmid=="26001779",]$assay="RNA-Seq"
result_df[result_df$pmid=="25946710",]$assay="RNA-Seq"
result_df[result_df$pmid=="25700344",]$assay="DNA gene panel,RNA-Seq"
result_df[result_df$pmid=="25638290",]$assay="DNA gene panel"
result_df[result_df$pmid=="25638290",]$assay="WES"
result_df[result_df$pmid=="25498120",]$assay="DNA gene panel"
result_df[result_df$pmid=="25335895",]$assay="RNA-Seq"
result_df[result_df$pmid=="30481710",]$assay="Metagenomics"
result_df[result_df$pmid=="29465611",]$assay="DNA gene panel"
result_df[result_df$pmid=="31465725",]$assay="RNA-Seq"
result_df[result_df$pmid=="29997562",]$assay="DNA gene panel"
result_df[result_df$pmid=="22472776",]$assay="DNA gene panel"
result_df[result_df$pmid=="22294635",]$assay="RNA-Seq"
result_df[result_df$pmid=="29163767",]$assay="DNA gene panel"
result_df[result_df$pmid=="30049532",]$assay="RNA-Seq"
result_df[result_df$pmid=="30746468",]$assay="RNA-Seq"
result_df[result_df$pmid=="31682074",]$assay="RNA-Seq"
result_df[result_df$pmid=="32265907",]$assay="RNA-Seq"
result_df[result_df$pmid=="32237059",]$assay="RNA-Seq"
result_df[result_df$pmid=="32199921",]$assay="DNA gene panel"
result_df[result_df$pmid=="32159782",]$assay="WGS"
result_df[result_df$pmid=="32115259",]$assay="RNA-Seq"

#delete entries
#32047518 Necrotizing Sarcoid Granulomatosis (NGS)
#31586650,31084224,30774014,30054207 no access to look into methods
#30679579,9607210,10750629,11317015,16227418,15022312,28417081,27421624,26984631,26275432,25828588,31443722,31347786,11396511 no NGS
#30247649 GWAS
#31858722,24070009,25598697,31652165 review
#29607230 methods missing
#29710329 meninigitis
#26110394,25990289 microarray
#31006511 not using Rheuma seq data
#30906878,29578063,29180447,28686088 no NGS
#30337930,31745718,27788994 no rheuma
#30049830
#26404574 commenmtray, no acces
#26284111	PERSPECTIVE ARTICLE
#29248936 not clear from the methods, if they did ngs
# 24065968 not clear from the methods what kimd of NGS they downloaded
#18037627,18390570,20842748, 22627728, 22402741, 22275833, 26705675 no NGS
#31033783,30100989 oncology
#31021403 periodontal diagnosis
#31008566 nominal groups (NGs)
# 30946936,25799352,25541309,25373782	no access
#30643044 retracted
#30528862,30441946, 31494232 no access and no signs of NGS in Abstract
#30463656 article in chinese
#30211458 erratum
#29920970 goose-origin astrovirus
#29769526 Neuromyelitis optica
#29674726 astrovirus in goslings, 28012117 in chicken
#29257335 Re-analysis of DNA microarray data
#28286571 Commentary about NGS and gut microbiome
#27325253,32264838 not clear from the case reports (no method section), what they did
#26182405 Aicardi-Goutières-Syndrom
#22933063,30632919 review
#31744152 chicken gout
deletepmids=c("32047518","32264838","31494232","27788994","11396511","31744152","30632919","22933063","25373782","25541309","25799352","26182405","26705675","27325253","28012117","28286571","28417081","28686088","29257335","29674726","29769526","29920970","30100989","30211458","30441946","30463656","30528862","30643044","30946936","31008566","31021403","31033783","31652165","31745718","26284111","26404574","22402741","25598697","24065968","24070009","22275833","22627728","18037627","18390570","20842748","16227418","15022312","11317015","10750629","9607210","30054207","30774014","31084224","31586650","29180447","29248936","29578063","30049830","30337930","30679579","30247649","31858722","29607230","29710329","27421624","26984631","26275432","26110394","25990289","25828588","31443722","31347786","31006511","30906878")
result_df=result_df[(!result_df$pmid %in% deletepmids),]
#Aicardi-Goutières syndrome–like
result_df[result_df$pmid=="31874111" & result_df$disease=="Gout",]=NA
result_df <- na.omit(result_df)

#Rename categories
result_df[result_df$assay=="RNA-Seq,WGS",]$assay="WGS,RNA-Seq"
result_df[result_df$assay=="WES,WGS",]$assay="WGS,WES"
result_df[result_df$assay=="DNA gene panel",]$assay="Targeted DNA Seq"
result_df[result_df$assay=="DNA gene panel,RNA-Seq",]$assay="Targeted DNA Seq,RNA-Seq"
result_df$assay_main=result_df$assay

#missed by this search
new_row=data.frame("Rheumatoid Arthritis","23497938","PhIP-Seq","Journal of Autoimmunity","2013","Other")

names(new_row)=names(result_df)
newdf=rbind(result_df,new_row)
result_df=newdf

result_df[result_df$assay=="Metagenomics",]$assay_main="Other"
result_df[result_df$assay=="Epigenomics",]$assay_main="Other"
result_df[result_df$assay=="Epigenomics,Epigenomics,RNA-Seq",]$assay="Epigenomics,RNA-Seq"
result_df[result_df$assay=="Epigenomics,RNA-Seq",]$assay_main="Other"
result_df[result_df$assay=="Epigenomics,RNA-Seq,WES",]$assay_main="Other"
result_df[result_df$assay=="Epigenomics,WGS",]$assay_main="Other"
result_df[result_df$assay=="Metagenomics,RNA-Seq",]$assay_main="Other"
result_df[result_df$assay=="Methylation,RNA-Seq",]$assay_main="Other"
result_df[result_df$assay=="PhIP-Seq",]$assay_main="Other"
result_df[result_df$disease=="Systemic Lupus Erythematosus",]$disease="SLE"
result_df[result_df$disease=="Familial Mediterranean Fever",]$disease="FMF"


#-------redo all figures from paper--------------------------------
ggplot(data=result_df, aes(x=fct_infreq(disease),fill=assay_main)) + geom_bar(position="stack",stat="count") + 
  labs(y="# publications on pubmed") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13, colour="black"), 
        axis.title.x = element_blank(), axis.text.y=element_text(colour="black", size = 13), axis.title.y=element_text(colour="black", size = 14)) +
  theme(legend.title = element_blank()) +
  scale_fill_brewer(palette = "Paired") + scale_y_continuous(breaks=seq(0,110,10)) +
  theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.grid.minor = element_blank(), panel.background = element_blank())

ggplot(data=subset(result_df,assay_main=="Other"), aes(x=fct_infreq(disease),fill=assay)) + geom_bar(position="stack",stat="count") + 
  labs(y="# publications on pubmed") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13, colour="black"), 
        axis.title.x = element_blank(), axis.text.y=element_text(colour="black", size = 13), axis.title.y=element_text(colour="black", size = 14)) +
  theme(legend.title = element_blank()) +
  scale_fill_brewer(palette = "Paired") +
  theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.grid.minor = element_blank(), panel.background = element_blank())

result_df_unique=result_df
result_df_unique$disease=NULL
result_df_unique=unique(result_df_unique)

ggplot(data=result_df_unique, aes(x=year,fill=assay_main)) + geom_bar(position="stack",stat="count") + 
  labs(y="# unique publications on pubmed") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13, colour="black"), 
        axis.title.x = element_blank(), axis.text.y=element_text(colour="black", size = 13), axis.title.y=element_text(colour="black", size = 14)) +
  theme(legend.title = element_blank()) +
  scale_fill_brewer(palette = "Paired") + scale_y_continuous(breaks=seq(0,160,10)) +
  theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.grid.minor = element_blank(), panel.background = element_blank())



ggplot(data=subset(result_df_unique,journal %in% journal_freq[journal_freq$freq>=5,]$journal), aes(x=fct_infreq(journal),fill=assay_main)) + geom_bar(position="stack",stat="count") + 
  labs(y="# unique publications on pubmed") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13, colour="black"), 
        axis.title.x = element_blank(), axis.text.y=element_text(colour="black", size = 13), axis.title.y=element_text(colour="black", size = 14)) +
  theme(legend.title = element_blank()) +
  scale_fill_brewer(palette = "Paired") + scale_y_continuous(breaks=seq(0,40,10)) +
  theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.grid.minor = element_blank(), panel.background = element_blank())




#interesting studies
#2509162 DNA methylation and mRNA and microRNA expression of SLE CD4+ T cells correlate with disease phenotype.
#32042094 HLA typing in RA patients, 29614084 HLA typing to identify repsonse to adalimumab treatment
#22933063 first mention of NGS
#31428656,31110316 scRNA-Seq
#31344123 re-use of public data, 
#30319641,30700427,mouse 27183593 (miRNA) 
#29892977,31779271,28339495 miRNA, 30836987 lnc + miRNA, 28808260 lncRNA inclusion body moysitis
#22294635,29891556,29657145,26245356,26815131,26227771,29163767,31682074 TCR, 27420474 (after nivolumab!), 26713667 (identification T cell epitope ), 26001779 (SLE longitudinale)), 25946710 (gout)
#29432186,31988805,31882654,23606709,28053320 BCR, 27273876 antibody repertiore
#29659823 WES and then zebrafish
#28750028 targeted DNA
#28728565 CD4+ CD14+ cells
#28523199,31538826 DNA gene panel case report
#31993940 case report WES
#31832070 case report gene panel
#31647025 RNA-Seq for RA and response to anti-TNF-a!!!
#31608065 circ_RNA Biomarker SLE, 29247798 circ_RNA mouse model osteoarthritis
#28532706 WES mutation, Blau Syndrome


TableForPub=NULL
for (key in abstracts_list){
  pubmedid=custom_grep(key,"PMID","char")
  if (pubmedid[1] %in% result_df$pmid){
    title=custom_grep(key,"ArticleTitle","char")
    abstract = custom_grep(key,"Abstract","char")
    TableForPub=rbind(TableForPub,c(pubmedid[1],title,abstract))
    }
}

#create table with title and abstract of all identified publications
TableForPub=as.data.frame(TableForPub)
names(TableForPub)=c("pmid","title","abstract")
setdiff(result_df$pmid,TableForPub$pmid)
new_row_1=data.frame("31465725","Novel Findings From Determination of Common Expressed Plasma Exosomal microRNAs in Patients With Psoriatic Arthritis, Psoriasis Vulgaris, Rheumatoid Arthritis, and Gouty Arthritis","Background: Circulating exosomal microRNAs modulate not only cancer cell metabolism but also the immune response, and therefore plasma exosomal microRNAs might have the potential to be the biomarkers for a number of immune disorders. Objective: This study was conducted to identify the common mechanisms among psoriatic arthritis (PsA), psoriasis vulgaris (PV), rheumatoid arthritis (RA), and gouty arthritis (GA). The common expressed plasma exosomal microRNAs in these diseases were determined. Methods: The expression of microRNAs derived from plasma exosome of patients with PsA (n=30), PV (n=15), RA (n=15), GA (n=15), and healthy controls (n=15) was evaluated via sequencing. Function analysis of common expressed microRNAs was conducted by the Gene ontology (GO) and Kyoto encyclopedia of genes and genomes (KEGG) enrichment analyses. Coexpression analysis was conducted to identify novel and significant genes and proteins by using the Search Tool for the Retrieval of Interacting Genes (STRING). A systematic literature review was conducted to uncover the role of the common microRNAs in the pathogenesis of PsA, PV, RA, and GA. Results: A total of 36 common expressed microRNAs were detected in patients with PsA, PV, RA, and GA. The most significantly enriched biological processes, cellular components, and molecular functions were homophilic cell adhesion via plasma membrane adhesion molecules, CCR4-NOT complex, and calcium ion binding, respectively. Antigen processing and presentation was the most significantly enriched pathway. A total of 91 validated coexpressed gene pairs were identified and 16 common expressed microRNAs and 85 potential target genes were screened based on Cytoscape. Of 36 common expressed microRNAs, 5 microRNAs, including hsa-miR-151a-3p, hsa-miR-199a-5p, hsa-miR-370-3p, hsa-miR-589-5p, and hsa-miR-769-5p, were considered to be connected with the common pathogenesis of PsA, PV, RA, and GA. Systemic review revealed that the roles of these 5 microRNAs are related to immune disorder and bone injury, which matches the conclusion from GO and KEGG analyses. Conclusion: (1) Both immune disorder and bone metabolic dysregulation could be the shared mechanism in the development of PsA, PV, RA, and GA. (2) Immune dysfunction is involved in GA. Our study may shed new light on the diagnosis and treatment strategy of these autoimmune diseases and GA, which warrants further studies.")
names(new_row_1)=names(TableForPub)
TableForPub=rbind(TableForPub,new_row_1)
new_row_2=data.frame("23497938","PhIP-Seq Characterization of Autoantibodies From Patients With Multiple Sclerosis, Type 1 Diabetes and Rheumatoid Arthritis","Autoimmune disease results from a loss of tolerance to self-antigens in genetically susceptible individuals. Completely understanding this process requires that targeted antigens be identified, and so a number of techniques have been developed to determine immune receptor specificities. We previously reported the construction of a phage-displayed synthetic human peptidome and a proof-of-principle analysis of antibodies from three patients with neurological autoimmunity. Here we present data from a large-scale screen of 298 independent antibody repertoires, including those from 73 healthy sera, using phage immunoprecipitation sequencing. The resulting database of peptide-antibody interactions characterizes each individual's unique autoantibody fingerprint, and includes specificities found to occur frequently in the general population as well as those associated with disease. Screening type 1 diabetes (T1D) patients revealed a prematurely polyautoreactive phenotype compared with their matched controls. A collection of cerebrospinal fluids and sera from 63 multiple sclerosis patients uncovered novel, as well as previously reported antibody-peptide interactions. Finally, a screen of synovial fluids and sera from 64 rheumatoid arthritis patients revealed novel disease-associated antibody specificities that were independent of seropositivity status. This work demonstrates the utility of performing PhIP-Seq screens on large numbers of individuals and is another step toward defining the full complement of autoimmunoreactivities in health and disease.")
names(new_row_2)=names(TableForPub)
TableForPub=rbind(TableForPub,new_row_2)

#export both result tables 
write.csv2(result_df,"~/Documents/Lupus_NGS/main_results.csv")
write.csv2(TableForPub,"~/Documents/Lupus_NGS/title_and_abstracts.csv")


#2nd analysis----------------------------------------------------
#SRA datasets
sra_datasets <- read.delim("~/Documents/Lupus_NGS/sra_datasets.txt")
sra_datasets=as.data.frame(sra_datasets)
sra_datasets$assay=as.character(sra_datasets$assay)
sra_datasets[sra_datasets$assay=="miRNA-Seq",]$assay="miRNA/ncRNA-Seq"
sra_datasets[sra_datasets$assay=="ncRNA-Seq",]$assay="miRNA/ncRNA-Seq"
#Methylation=Bisulfite-Seq, MEDIP-Seq/MRE-Seq
sra_datasets[sra_datasets$assay=="Bisulfite-Seq",]$assay="Methylation"
sra_datasets[sra_datasets$assay=="MEDIP-Seq/MRE-Seq",]$assay="Methylation"
sra_datasets$disease=as.character(sra_datasets$disease)
sra_datasets$samples=as.integer(sra_datasets$samples)
sra_datasets[sra_datasets$disease=="myositis, polymyositis, dermatomyositis",]$disease="(poly/derma)myositis"
sra_datasets$disease=stringr::str_to_title(sra_datasets$disease)
SRA_tissues <- read.delim("~/Documents/Lupus_NGS/SRA_tissues.txt")
#---plot figures
ggplot(data=sra_datasets, aes(x=reorder(disease,-samples),y=samples,fill=assay)) + geom_bar(position="stack",stat="identity") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13, colour="black"), 
        axis.title.x = element_blank(), axis.text.y=element_text(colour="black", size = 13), axis.title.y=element_text(colour="black", size = 14)) +
  theme(legend.title = element_blank()) +
  scale_fill_brewer(palette = "Paired") + scale_y_continuous(breaks=seq(0,7000,1000)) +
  theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.grid.minor = element_blank(), panel.background = element_blank())
ggplot(data=SRA_tissues, aes(x="",y=samples,fill=tissue)) + geom_bar(stat="identity",position = position_fill()) + 
  geom_text(aes(label = samples), position = position_fill(vjust = 0.5)) + coord_polar(theta = "y") + 
  facet_wrap(~ disease) + theme(axis.title.x = element_blank(),axis.title.y = element_blank())

