library(easyPubMed)
library(tidyverse)
library(data.table)
library(ggplot2)
library(dplyr)

#search performed August 15, 2020
pubmed_query_string="(methylomics OR epigenomics OR NGS OR \"next generation sequencing\" OR RNA-Seq OR \"mRNA sequencing\" OR \"RNA sequencing\" OR \"RNA-sequencing\" OR \"transcriptome sequencing\" OR \"whole exome sequencing\" OR \"whole-exome sequencing\" OR \"high throughput sequencing\" OR \"high-throughput sequencing\" OR \"DNA sequencing\" OR \"RNA sequencing\" OR \"RNA-sequencing\" OR \"DNA-sequencing\" OR WXS OR WGS OR \"whole-genome sequencing\" OR \"whole genome sequencing\") AND (rheumatology OR \"rheumatologic disease\" OR \"rheumatologic disease\" ))"


pubmed_ids <- get_pubmed_ids(pubmed_query_string)
#1097 entries
abstracts_xml <- fetch_pubmed_data(pubmed_ids,retmax = 5000)
abstracts_list <- articles_to_list(abstracts_xml)
saveRDS(abstracts_list,file="/Users/sebastian/pubmed_rheuma_HTS/abstracts_list_for_keywords.RDS")
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
write.csv2(finalTable_freq,"/Users/sebastian/pubmed_rheuma_HTS/finalTable_freq.csv")
#326 entries in finalTable_freq


#search performed Sep,4, 2020

pubmed_query_string="(methylomics OR epigenomics OR NGS OR \"next generation sequencing\" OR RNA-Seq OR \"mRNA sequencing\" OR \"RNA sequencing\" OR \"RNA-sequencing\" OR \"transcriptome sequencing\" OR \"whole exome sequencing\" OR \"whole-exome sequencing\" OR \"high throughput sequencing\" OR \"high-throughput sequencing\" OR WXS OR WGS OR \"whole-genome sequencing\" OR \"whole genome sequencing\") AND (\"autoinflammatory syndrome\" OR dermatomyositis OR enthesitis OR \"familial mediterranean fever\" OR \"granulomatosis with polyangiitis\" OR \"juvenile idiopathic arthritis\" OR myositis OR osteoarthritis OR polymyositis OR \"psoriatic arthritis\" OR \"rheumatoid arthritis\" OR sacroiliitis OR \"sjögren syndrome\" OR \"sjögren's syndrome\" OR spondyloarthritis OR synovitis OR \"systemic lupus erythematosus\" OR \"systemic sclerosis\" OR vasculitis OR uveitis OR gout OR polychondritis)"



pubmed_ids <- get_pubmed_ids(pubmed_query_string)
pubmed_ids$Count
#1162 hits
abstracts_xml <- fetch_pubmed_data(pubmed_ids,retmax = 5000)
abstracts_list <- articles_to_list(abstracts_xml)
diseases=c("autoinflammatory syndrome","dermatomyositis","enthesitis","familial mediterranean fever","granulomatosis with polyangiitis","juvenile idiopathic arthritis","myositis","osteoarthritis","polymyositis","psoriatic arthritis","rheumatoid arthritis","sacroiliitis","sjögren syndrome","sjögren's syndrome","spondyloarthritis","synovitis","systemic lupus erythematosus","systemic sclerosis","vasculitis","uveitis","gout","polychondritis")
RNA=c("rna sequencing", "transcriptome sequencing", "rna-seq","rna-based next-generation sequencing","mrna profiles","rna-sequencing")
WXS=c("whole exome sequencing","whole-exome sequencing","whole-exome-sequencing","amplicon sequencing")
WGS=c("whole-genome sequencing","whole genome sequencing","whole-genome shotgun sequencing")
Bacteria=c("16s","metagenomics")
Epigenomic=c("atac-seq","chip-seq","wgbs","methylomics","methylomic","epigenomics")
single=c("single cell","single-cell")

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
    for (sin in single){
      if (length(grep(sin,title_abstract_key))>0){assay=paste(assay,"single",sep=",")}
    }
    for (disease in diseases){
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
result_df$journal=as.character(result_df$journal)
result_df$disease=as.character(result_df$disease)
result_df$disease=stringr::str_to_title(result_df$disease)

#manually adding missing data by looking at each publication using the PMID

#995 entries, 847 unique Pubmed-IDs
result_df$pmid=as.character(result_df$pmid)
result_df[result_df$pmid=="30944248",]$assay="RNA-Seq,single"
result_df[result_df$pmid=="31428656",]$assay="RNA-Seq,single"
result_df[result_df$pmid=="31370803",]$assay="Metagenomics"
result_df[result_df$pmid=="31360262",]$assay="RNA-Seq"
result_df[result_df$pmid=="31344123",]$assay="RNA-Seq"
result_df[result_df$pmid=="31337345",]$assay="DNA gene panel"
result_df[result_df$pmid=="31101814",]$assay="WGS,WES"
result_df[result_df$pmid=="30783801",]$assay="DNA gene panel"
result_df[result_df$pmid=="30513227",]$assay="DNA gene panel"

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
result_df[result_df$pmid=="31848804",]$assay="DNA gene panel"
result_df[result_df$pmid=="31993940",]$assay="WES"
result_df[result_df$pmid=="31988805",]$assay="RNA-Seq"
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
result_df[result_df$pmid=="31237906",]$assay="RNA-Seq"
result_df[result_df$pmid=="30987788",]$assay="DNA gene panel"
result_df[result_df$pmid=="30850477",]$assay="RNA-Seq"
result_df[result_df$pmid=="30836987",]$assay="RNA-Seq"
result_df[result_df$pmid=="30724444",]$assay="RNA-Seq"
result_df[result_df$pmid=="30700427",]$assay="Metagenomics"
result_df[result_df$pmid=="30544699",]$assay="RNA-Seq"
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
result_df[result_df$pmid=="28808260",]$assay="RNA-Seq"
result_df[result_df$pmid=="28597968",]$assay="WGS,WES"
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
result_df[result_df$pmid=="29997562",]$assay="DNA gene panel"
result_df[result_df$pmid=="22472776",]$assay="DNA gene panel"
result_df[result_df$pmid=="22294635",]$assay="RNA-Seq"
result_df[result_df$pmid=="30746468",]$assay="RNA-Seq"
result_df[result_df$pmid=="31682074",]$assay="RNA-Seq"
result_df[result_df$pmid=="32265907",]$assay="RNA-Seq"
result_df[result_df$pmid=="32237059",]$assay="RNA-Seq"
result_df[result_df$pmid=="32199921",]$assay="DNA gene panel"
result_df[result_df$pmid=="32159782",]$assay="WGS"
result_df[result_df$pmid=="32115259",]$assay="RNA-Seq"
result_df[result_df$pmid=="24136464",]$assay="DNA gene panel"
result_df[result_df$pmid=="27402083",]$assay="RNA-Seq"
result_df[result_df$pmid=="28053302",]$assay="DNA gene panel"
result_df[result_df$pmid=="28053320",]$assay="RNA-Seq"
result_df[result_df$pmid=="28532706",]$assay="WES"
result_df[result_df$pmid=="29472286",]$assay="WES" #case report
result_df[result_df$pmid=="29500522",]$assay="DNA gene panel"
result_df[result_df$pmid=="29040051",]$assay="DNA gene panel" #case reports
result_df[result_df$pmid=="29774027",]$assay="Epigenomics" #bisulfite amplicon seq using NGS
result_df[result_df$pmid=="30962246",]$assay="RNA-Seq" #PRJNA427177
result_df[result_df$pmid=="31101603",]$assay="DNA gene panel" #TCR receptor seq
result_df[result_df$pmid=="31882654",]$assay="RNA-Seq" #BCR repertoire RA
result_df[result_df$pmid=="31926583",]$assay="RNA-Seq" #no making data available
result_df[result_df$pmid=="32191636",]$assay="RNA-Seq" #GEO but embargo
result_df[result_df$pmid=="32196497",]$assay="RNA-Seq" #re-used RNA-Seq data from e.g. SRA
result_df[result_df$pmid=="32332704",]$assay="RNA-Seq" #no making data available
result_df[result_df$pmid=="32365362",]$assay="RNA-Seq" #no making data available
result_df[result_df$pmid=="32403239",]$assay="RNA-Seq" #no making data available
result_df[result_df$pmid=="32471379",]$assay="DNA gene panel" #case reports
result_df[result_df$pmid=="32510848",]$assay="DNA gene panel"
result_df[result_df$pmid=="32518584",]$assay="RNA-Seq" #RA, datasets available on request
result_df[result_df$pmid=="32552384",]$assay="DNA gene panel" #case reports
result_df[result_df$pmid=="32560314",]$assay="RNA-Seq" #PsA, no making data available
result_df[result_df$pmid=="32565918",]$assay="Metagenomics"
result_df[result_df$pmid=="32607330",]$assay="DNA gene panel" 
result_df[result_df$pmid=="32659156",]$assay="WES" #case report
result_df[result_df$pmid=="32714036",]$assay="Metagenomics"
result_df[result_df$pmid=="32714036",]$assay="Metagenomics" #case report Virus identification
result_df[result_df$pmid=="32735477",]$assay="Metagenomics" #case report Virus identification
result_df[result_df$pmid=="32745911",]$assay="RNA-Seq" #OA, no making data available
result_df[result_df$pmid=="32807082",]$assay="Metagenomics" #case report Virus identification
result_df[result_df$pmid=="32818496",]$assay="RNA-Seq" #SLE (LN), https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155405 (Mouse)
result_df[result_df$pmid=="31651337",]$assay="Epigenomics" #whole genome bisulfite seq using HiSeq4000
result_df[result_df$pmid=="32438470",]$assay="Epigenomics" #bioinformatics whole genome bisulfite seq 
result_df[result_df$pmid=="31926583",]$pmid=31465725

#Epigenomics=bisulfite seq using NGS, ATAC-Seq
#DNA gene panel=targeted re-sequencing

#delete entries
#no information about NGS assay: 24065968,29607230,32264838,23942639
#no NGS: 21262369,22046267 (GWAS),20514598, 29257335 microarray, 28686088,22402741,22275833,22627728,18037627,18390570,20842748,16227418,15022312,11317015,10750629,26705675,30679579,27421624,
#26984631,25990289 microarray,25828588,28417081,31347786,30906878,12522564,21192791,21436623,21865261,22014533 ChIP microarray,22300536, 23340290,23353785,23950730 (Methylation400 assay),
#24812288, 25261579 (HumanMethylation450 BeadChip arrays), 25537884,25661834 (microarray), 25789080, 26079624, 26492816, 26556652(HumanMethylation450 BeadChips), 26628598(Microarray), 27325182, 27818202
#27876072, 27885849(HumanMethylation450 BeadChips),28081217,28430662,28475762(HumanMethylation450 BeadChips), 28785300(GWAS), 29018507(BS-Seq),26111028,26787370(BS-Seq),29422534,30385847,
#30598166,30779748,31467281(llumina MethylationEPIC BeadChip),31513634(Infinium MethylationEPIC BeadChip array),32509015(microarray),25414732,29854846,32641481
#32581359:GWAS,24445253, 24223582 GWAS,24013639 GWAS
#not disease defined in 1st round:
#29769526, 31021403 periodontal diagnosis, 30100989,31033783 oncology, 31745718, 30337930,#26182405 Aicardi-Goutières-Syndrom,29710329 meninigitis,31006511 not using Rheuma seq data, 29434599(keratinocyes)
#32489448 glioma,31091335,31081255,30974145,30313048(myocard infarct),29608426(zika virus infection),29184082(cervical cancer),28920417
#no access and no info of assay in abstract: #30441946,31494232,30946936,30054207,30774014,31586650,25799352,25541309,25373782,25791244,25791245,29113828,30523707,31693422,31755746,30304699,31883828,32157911
#29674726 astrovirus in goslings, 28012117, 31744152 in chicken, #29920970 goose-origin astrovirus, 32532177 (pig)
#not english: 30463656,23981988,30124204 article in chinese
#retracted: 30643044
#not primary research article:#31652165,25598697,24070009,22933063,30632919,31858722 review, 26284111 PERSPECTIVE ARTICLE,28286571 Commentary about NGS and gut microbiome,21273237,
#23025638 & 23244304 Editorial,23275983, 25167330 Commentary,25366188,25557657, 26231343 (editorial), 27028588(meeting report),27482953(editorial),27607472(review),28447857(editorial),
#25165988(Commetary),28770613(review),29644081(review&meta-analysis),31699217,31876195,31902241(review&meta-analysis),32203190(review),32475529(book chapter, review),32611772(editorial),32746644(review),32806879(review)
#30399325(review)
deletepmids=c("32264838","31494232","31744152","30632919","22933063","25373782","25541309","25799352","26182405","26705675","28012117","28286571","28417081",
              "28686088","29257335","29674726","29769526","29920970","30100989","30441946","30463656","30643044","30946936","31021403","31033783","31652165",
              "31745718","26284111","22402741","25598697","24065968","24070009","22275833","22627728","18037627","18390570","20842748","16227418","15022312",
              "11317015","10750629","30054207","30774014","31586650","30337930","30679579","31858722","29607230","29710329","27421624","26984631","25990289",
              "25828588","31347786","31006511","30906878","12522564","21192791","21273237","21436623","21865261","22014533","22300536","23025638","23244304",
              "23275983","23340290","23353785","23942639","23950730","23981988","24812288","25167330","25261579","25366188","25537884","25557657","25661834",
              "25789080","25791244","25791245","26079624","26231343","26492816","26556652","26628598","27028588","27325182","27482953","27607472","27818202",
              "27876072","27885849","28081217","28430662","28447857","25165988","28475762","28770613","28785300","29018507","26787370","26111028","29113828",
              "29422534","29644081","30523707","31693422","31755746","30124204","30304699","30385847","30598166","30779748","31467281","31513634","31699217",
              "31876195","31883828","31902241","32157911","32203190","32475529","32509015","32532177","32611772","32746644","32806879","25414732","29854846",
              "32641481","29434599","32581359","32489448","31091335","31081255","30974145","30399325","30313048","29608426","29184082","28920417","24445253",
              "24223582","24013639","20514598","21262369","22046267")
result_df=result_df[(!result_df$pmid %in% deletepmids),]
#Aicardi-Goutières syndrome–like
result_df[result_df$pmid=="31874111" & result_df$disease=="Gout",]=NA
result_df <- na.omit(result_df)
#no NGS on SLE samples
result_df[result_df$pmid=="29720240" & result_df$disease=="SLE",]=NA
result_df <- na.omit(result_df)
#manually add studies identified via SRA project search
result_df=rbind(result_df,c("SLE",24645875,"RNA-Seq","Connective tissue research",2014,"RNA-Seq"))
result_df=rbind(result_df,c("SLE",31263277,"Epigenomics,RNA-Seq","Nature Immunology",2019,"Other"))
result_df=rbind(result_df,c("SLE",31890206,"RNA-Seq","Clinical & translational immunology",2019,"RNA-Seq"))
result_df=rbind(result_df,c("SLE",29717110,"RNA-Seq","Nature communications",2018,"RNA-Seq"))
result_df=rbind(result_df,c("SLE",28420548,"RNA-Seq","Journal of autoimmunity",2017,"RNA-Seq"))
result_df=rbind(result_df,c("SLE",30130253,"RNA-Seq","The Journal of clinical investigation",2018,"RNA-Seq"))
result_df=rbind(result_df,c("SLE",31200750,"RNA-Seq","Arthritis research & therapy",2019,"RNA-Seq"))
result_df=rbind(result_df,c("SLE",30478422,"Epigenomics,RNA-Seq","Nature immunology",2019,"Other"))
result_df=rbind(result_df,c("SLE",31616406,"RNA-Seq","Frontiers in immunology",2019,"RNA-Seq"))
result_df=rbind(result_df,c("RA",31616406,"RNA-Seq","Frontiers in immunology",2019,"RNA-Seq"))
result_df=rbind(result_df,c("Sjögren's Syndrome",31616406,"RNA-Seq","Frontiers in immunology",2019,"RNA-Seq"))

#Rename categories
result_df[result_df$assay=="RNA-Seq,WGS",]$assay="WGS,RNA-Seq"
result_df[result_df$assay=="Epigenomics,Epigenomics",]$assay="Epigenomics"
result_df[result_df$assay=="Epigenomics,Epigenomics,RNA-Seq",]$assay="Epigenomics,RNA-Seq"
result_df[result_df$assay=="WES,WGS",]$assay="WGS,WES"
result_df[result_df$assay=="single",]$assay="scRNA-Seq"
result_df[result_df$assay=="single,single",]$assay="scRNA-Seq"
result_df[result_df$assay=="RNA-Seq,single",]$assay="scRNA-Seq"
result_df[result_df$assay=="WES,single",]$assay="WES,scRNA-Seq"
result_df[result_df$assay=="RNA-Seq,single,single",]$assay="scRNA-Seq"
result_df[result_df$assay=="DNA gene panel",]$assay="Targeted DNA Seq"
result_df[result_df$assay=="DNA gene panel,RNA-Seq",]$assay="Targeted DNA Seq,RNA-Seq"
result_df$assay_main=result_df$assay


result_df[result_df$assay=="Metagenomics",]$assay_main="Other"
result_df[result_df$assay=="Metagenomics,RNA-Seq,WES",]$assay_main="Other"
result_df[result_df$assay=="Epigenomics",]$assay_main="Other"
result_df[result_df$assay=="Epigenomics,RNA-Seq",]$assay_main="Other"
result_df[result_df$assay=="Epigenomics,RNA-Seq,WES",]$assay_main="Other"
result_df[result_df$assay=="Epigenomics,WGS",]$assay_main="Other"
result_df[result_df$assay=="Metagenomics,RNA-Seq",]$assay_main="Other"

result_df[result_df$assay=="PhIP-Seq",]$assay_main="Other"
result_df[result_df$disease=="Systemic Lupus Erythematosus",]$disease="SLE"
result_df[result_df$disease=="Familial Mediterranean Fever",]$disease="FMF"
result_df[result_df$disease=="Rheumatoid Arthritis",]$disease="RA"
result_df[result_df$disease=="Juvenile Idiopathic Arthritis",]$disease="JIA"
result_df[result_df$disease=="Autoinflammatory Syndrome",]$disease="AutoSyn"
result_df[result_df$disease=="Granulomatosis With Polyangiitis",]$disease="GPA"

saveRDS(result_df,file="/Users/sebastian/pubmed_rheuma_HTS/result_df_11Sep2020.RDS")

#-------redo all figures from paper--------------------------------
#Figure 3: diseases vs assay
df2=result_df%>%
  group_by(disease)%>%
  count(disease,name="number")
ggplot(data=result_df, aes(x=fct_infreq(disease),fill=assay_main)) + geom_bar(position="stack",stat="count") + 
  labs(y="# publications on pubmed") +
  geom_text(data = df2, aes(x=disease,y=200,label=number,angle=60),inherit.aes = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13, colour="black"), 
        axis.title.x = element_blank(), axis.text.y=element_text(colour="black", size = 13), axis.title.y=element_text(colour="black", size = 14)) +
  theme(legend.title = element_blank()) +
  scale_fill_brewer(palette = "Paired") + scale_y_continuous(breaks=seq(0,200,10)) +
  theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.grid.minor = element_blank(), panel.background = element_blank())

#Figure S3. disease vs other_assays
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

#S2: scRNA-Seq vs year
ggplot(data=subset(result_df_unique, result_df_unique$assay %in% c("scRNA-Seq","WES,scRNA-Seq")), aes(x=year,fill=assay)) + geom_bar(position="stack",stat="count") + 
  labs(y="# unique publications on pubmed") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13, colour="black"), 
        axis.title.x = element_blank(), axis.text.y=element_text(colour="black", size = 13), axis.title.y=element_text(colour="black", size = 14)) +
  theme(legend.title = element_blank()) +
  scale_fill_brewer(palette = "Paired") + scale_y_continuous(breaks=seq(0,25,5)) +
  theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.grid.minor = element_blank(), panel.background = element_blank())

#Figure 2: publications vs year
df2=result_df_unique%>%
  group_by(year)%>%
  count(year,name="number")
abundance_plot=ggplot(data=result_df_unique, aes(x=year,fill=assay_main)) + geom_bar(position="stack",stat="count") + 
  labs(y="# unique publications on pubmed") +
  geom_text(data = df2, aes(x=year,y=200,label=number,angle=60),inherit.aes = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13, colour="black"), 
        axis.title.x = element_blank(), axis.text.y=element_text(colour="black", size = 13), axis.title.y=element_text(colour="black", size = 14)) +
  theme(legend.title = element_blank()) +
  scale_fill_brewer(palette = "Paired") + scale_y_continuous(breaks=seq(0,400,10)) +
  theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.grid.minor = element_blank(), panel.background = element_blank())
abundance_plot

#FigureS1: exponential growth
occurences=table(unlist(result_df_unique$year))
year_abundance=as.data.frame(occurences[9:18])
names(year_abundance)=c("year","abundance")
df <- cbind(year_abundance, index = 1:nrow(year_abundance))
fit <- lm(log(abundance) ~ index, data = df)
abundance_plot + stat_function(fun = function(x) exp(fit$coefficients[1] + x*fit$coefficients[2]))

#------journals vs assays
journal_freq=table(result_df_unique$journal)
journal_freq=as.data.frame(journal_freq)
names(journal_freq)=c("journal","freq")
journal_freq$freq=as.numeric(journal_freq$freq)

ggplot(data=subset(result_df_unique,journal %in% journal_freq[journal_freq$freq>=8,]$journal), aes(x=fct_infreq(journal),fill=assay_main)) + geom_bar(position="stack",stat="count") + 
  labs(y="# unique publications on pubmed") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13, colour="black"), 
        axis.title.x = element_blank(), axis.text.y=element_text(colour="black", size = 13), axis.title.y=element_text(colour="black", size = 14)) +
  theme(legend.title = element_blank()) +
  scale_fill_brewer(palette = "Paired") + scale_y_continuous(breaks=seq(0,60,10)) +
  theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.grid.minor = element_blank(), panel.background = element_blank())



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


#export both result tables 
write.csv2(result_df,"~/pubmed_rheuma_HTS/main_results.csv")
write.csv2(TableForPub,"~/pubmed_rheuma_HTS/title_and_abstracts.csv")

#number of assays in unique publications
length(unique(result_df[result_df$assay %in% c("RNA-Seq","Epigenomics,RNA-Seq","Epigenomics,RNA-Seq,WES","Metagenomics,RNA-Seq","Metagenomics,RNA-Seq,WES","RNA-Seq,WES","scRNA-Seq","Targeted DNA Seq,RNA-Seq","WES,scRNA-Seq","WGS,RNA-Seq"),]$pmid))
length(unique(result_df[result_df$assay %in% c("scRNA-Seq","WES,scRNA-Seq"),]$pmid))
length(unique(result_df[result_df$assay %in% c("WES","Epigenomics,RNA-Seq,WES","Metagenomics,RNA-Seq,WES","RNA-Seq,WES","WES,scRNA-Seq","WGS,WES","Targeted DNA Seq","Targeted DNA Seq,RNA-Seq"),]$pmid))
length(unique(result_df[result_df$assay %in% c("Metagenomics","Metagenomics,RNA-Seq,WES","Metagenomics,RNA-Seq"),]$pmid))
length(unique(result_df[result_df$assay %in% c("WGS","Epigenomics,WGS","WGS,RNA-Seq","WGS,WES"),]$pmid))
length(unique(result_df[result_df$assay %in% c("Epigenomics","Epigenomics,WGS","Epigenomics,RNA-Seq","Epigenomics,RNA-Seq,WES"),]$pmid))
#2nd analysis----------------------------------------------------
#no datasets: FMF, enthesitis, polychondritis
#gout: only cancer datasets
#-----------Figure 4 -------SRA datasets
sra_datasets <- read.delim("/Users/sebastian/pubmed_rheuma_HTS/sra_assays.tsv")
sra_datasets=as.data.frame(sra_datasets)
sra_datasets$assay=as.character(sra_datasets$assay)
sra_datasets[sra_datasets$assay=="chip",]$assay="Epigenomics"
sra_datasets[sra_datasets$assay=="chip-seq",]$assay="Epigenomics"
sra_datasets[sra_datasets$assay=="mirna-seq",]$assay="miRNA/ncRNA-Seq"
sra_datasets[sra_datasets$assay=="ncrna-seq",]$assay="miRNA/ncRNA-Seq"
#Methylation=Bisulfite-Seq, MEDIP-Seq/MRE-Seq
sra_datasets[sra_datasets$assay=="bisulfite-seq",]$assay="Epigenomics"
sra_datasets[sra_datasets$assay=="atac-seq",]$assay="Epigenomics"
sra_datasets[sra_datasets$assay=="dnase-hypersensitivity",]$assay="Epigenomics"
sra_datasets[sra_datasets$assay=="medip-seq",]$assay="Epigenomics"
sra_datasets[sra_datasets$assay=="mre-seq",]$assay="Epigenomics"
sra_datasets[sra_datasets$assay=="amplicon",]$assay="Targeted-capture"
sra_datasets[sra_datasets$assay=="targeted-capture",]$assay="Targeted-capture"
sra_datasets[sra_datasets$assay=="rna-seq",]$assay="RNA-Seq"
sra_datasets[sra_datasets$assay=="wgs",]$assay="WGS"
sra_datasets[sra_datasets$assay=="wxs",]$assay="WXS"
sra_datasets[sra_datasets$assay=="tn-seq",]$assay="TN-Seq"
sra_datasets[sra_datasets$assay=="mbd-seq",]$assay="MBD-Seq"
sra_datasets[sra_datasets$assay=="hi-c",]$assay="Hi-C"
sra_datasets[sra_datasets$assay=="other",]$assay="Other"
sra_datasets$disease=as.character(sra_datasets$disease)
sra_datasets$samples=as.integer(sra_datasets$samples)
sra_datasets[sra_datasets$disease=="MyoPolyDerma",]$disease="(poly/derma)myositis"
sra_datasets[sra_datasets$disease=="SysSclerosis",]$disease="Systemic Sclerosis"
sra_datasets[sra_datasets$disease=="Sjoegren",]$disease="Sjögren's Syndrome"
df2=sra_datasets%>%
  group_by(disease)%>%
  summarise(samples=sum(samples))
sra_datasets$disease=factor(sra_datasets$disease,levels=df2[order(df2$samples,decreasing = TRUE),]$disease)
ggplot(data=sra_datasets, aes(x=disease,y=samples,fill=assay)) + geom_bar(position="stack",stat="identity") +
  geom_text(data = df2, aes(x=disease,y=9000,label=samples,angle=60),inherit.aes = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13, colour="black"), 
        axis.title.x = element_blank(), axis.text.y=element_text(colour="black", size = 13), axis.title.y=element_text(colour="black", size = 14)) +
  theme(legend.title = element_blank()) +
  scale_fill_brewer(palette = "Paired") + scale_y_continuous(breaks=seq(0,9000,1000)) +
  theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.grid.minor = element_blank(), panel.background = element_blank())
#Figure S5 ------------samples per study----------------------------
sra_studies <- read.delim("/Users/sebastian/pubmed_rheuma_HTS/sra_studies.tsv")
sra_studies=as.data.frame(sra_studies)
sra_studies$number=as.integer(sra_studies$number)
sra_studies$disease=as.character(sra_studies$disease)
sra_studies[sra_studies$disease=="MyoPolyDerma",]$disease="(poly/derma)myositis"
sra_studies[sra_studies$disease=="SysSclerosis",]$disease="Systemic Sclerosis"
sra_studies[sra_studies$disease=="Sjoegren",]$disease="Sjoegren's Syndrome"
df2=sra_studies%>%
  group_by(disease)%>%
  count(disease,name="number")
sra_studies$disease=factor(sra_studies$disease,levels=df2[order(df2$number,decreasing = TRUE),]$disease)
dataMedian <- summarise(group_by(sra_studies, disease), MD = median(number))
ggplot(data=sra_studies, aes(x=disease,y=number)) + geom_jitter() + geom_boxplot(alpha = 0.2,outlier.shape = NA) + scale_y_log10() +
  geom_text(data = df2, aes(x=disease,y=8000,label=number,angle=60),inherit.aes = FALSE) +
  geom_text(data = dataMedian, aes(disease, MD, label = MD), position = position_dodge(width = 0.8), size = 3, vjust = -0.5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13, colour="black"), 
        axis.title.x = element_blank(), axis.text.y=element_text(colour="black", size = 13), axis.title.y=element_text(colour="black", size = 14)) + ylab("number of samples per study") +
  theme(legend.title = element_blank()) +
  theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.grid.minor = element_blank(), panel.background = element_blank())
#FigureS9-------sequencing_instrument
sra_instruments <- read.delim("/Users/sebastian/pubmed_rheuma_HTS/sra_instrument.tsv")
sra_instruments=as.data.frame(sra_instruments)
sra_instruments$instrument=as.character(sra_instruments$instrument)
sra_instruments$disease=as.character(sra_instruments$disease)
sra_instruments$samples=as.integer(sra_instruments$samples)
sra_instruments[sra_instruments==0] <- NA
sra_instruments=sra_instruments[complete.cases(sra_instruments),]
sra_instruments[sra_instruments$disease=="MyoPolyDerma",]$disease="(poly/derma)myositis"
sra_instruments[sra_instruments$disease=="SysSclerosis",]$disease="Systemic Sclerosis"
sra_instruments[sra_instruments$disease=="Sjoegren",]$disease="Sjögren's Syndrome"
df_instrument=sra_instruments%>%
  group_by(disease)%>%
  summarise(samples=sum(samples))
sra_instruments$disease=factor(sra_instruments$disease,levels=df_instrument[order(df_instrument$samples,decreasing = TRUE),]$disease)
ggplot(data=sra_instruments, aes(x=disease,y=samples,fill=instrument)) + geom_bar(position="stack",stat="identity") +
  geom_text(data = df_instrument, aes(x=disease,y=9000,label=samples,angle=60),inherit.aes = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13, colour="black"), 
        axis.title.x = element_blank(), axis.text.y=element_text(colour="black", size = 13), axis.title.y=element_text(colour="black", size = 14)) +
  theme(legend.title = element_blank()) +
  scale_fill_brewer(palette = "Paired") + scale_y_continuous(breaks=seq(0,9000,1000)) +
  theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.grid.minor = element_blank(), panel.background = element_blank())


#Figure S10-------layout
sra_layout <- read.delim("/Users/sebastian/pubmed_rheuma_HTS/sra_layout.tsv")
sra_layout=as.data.frame(sra_layout)
sra_layout$layout=as.character(sra_layout$layout)
sra_layout$disease=as.character(sra_layout$disease)
sra_layout$samples=as.integer(sra_layout$samples)
#sra_layout[sra_layout==0] <- NA
#sra_layout=sra_layout[complete.cases(sra_layout),]
sra_layout[sra_layout$disease=="MyoPolyDerma",]$disease="(poly/derma)myositis"
sra_layout[sra_layout$disease=="SysSclerosis",]$disease="Systemic Sclerosis"
sra_layout[sra_layout$disease=="Sjoegren",]$disease="Sjögren's Syndrome"
df2=sra_layout%>%
  group_by(disease)%>%
  summarise(samples=sum(samples))
sra_layout$disease=factor(sra_layout$disease,levels=df2[order(sra_layout$samples,decreasing = TRUE),]$disease)

ggplot(data=sra_layout, aes(x=disease,y=samples,fill=layout)) + geom_bar(position="stack",stat="identity") +
  geom_text(data = df2, aes(x=disease,y=9000,label=samples, angle=60),inherit.aes = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13, colour="black"), 
        axis.title.x = element_blank(), axis.text.y=element_text(colour="black", size = 13), axis.title.y=element_text(colour="black", size = 14)) +
  theme(legend.title = element_blank()) +
  #scale_y_continuous(breaks=seq(0,9000,1000)) +
  scale_fill_brewer(palette = "Paired") + scale_y_continuous(breaks=seq(0,9000,1000)) +
  theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.grid.minor = element_blank(), panel.background = element_blank())

#----------source------
#SRA datasets
sra_source <- read.delim("/Users/sebastian/pubmed_rheuma_HTS/sra_source.tsv")
sra_source=as.data.frame(sra_source)
sra_source$source=as.character(sra_source$source)
sra_source$disease=as.character(sra_source$disease)
sra_source$samples=as.integer(sra_source$samples)
sra_source[sra_source$disease=="MyoPolyDerma",]$disease="(poly/derma)myositis"
sra_source[sra_source$disease=="SysSclerosis",]$disease="Systemic Sclerosis"
sra_source[sra_source$disease=="Sjoegren",]$disease="Sjögren's Syndrome"
df2=sra_source%>%
  group_by(disease)%>%
  summarise(samples=sum(samples))
sra_source$disease=factor(sra_source$disease,levels=df2[order(df2$samples,decreasing = TRUE),]$disease)
ggplot(data=sra_source, aes(x=disease,y=samples,fill=source)) + geom_bar(position="stack",stat="identity") +
  geom_text(data = df2, aes(x=disease,y=9000,label=samples,),inherit.aes = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13, colour="black"), 
        axis.title.x = element_blank(), axis.text.y=element_text(colour="black", size = 13), axis.title.y=element_text(colour="black", size = 14)) +
  theme(legend.title = element_blank()) +
  scale_fill_brewer(palette = "Paired") + scale_y_continuous(breaks=seq(0,9000,1000)) +
  theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.grid.minor = element_blank(), panel.background = element_blank())

#Figure 5-------tissues
sra_tissue <- read.delim("/Users/sebastian/pubmed_rheuma_HTS/sra_tissue.tsv")
sra_tissue=as.data.frame(sra_tissue)
sra_tissue$tissue=as.character(sra_tissue$tissue)
sra_tissue$disease=as.character(sra_tissue$disease)
sra_tissue$samples=as.integer(sra_tissue$samples)
sra_tissue[sra_tissue==0] <- NA
sra_tissue=sra_tissue[complete.cases(sra_tissue),]
sra_tissue[sra_tissue$disease=="MyoPolyDerma",]$disease="(poly/derma)myositis"
sra_tissue[sra_tissue$disease=="SysSclerosis",]$disease="Systemic Sclerosis"
sra_tissue[sra_tissue$disease=="Sjoegren",]$disease="Sjögren's Syndrome"
df2=sra_tissue%>%
  group_by(disease)%>%
  summarise(samples=sum(samples))
sra_tissue$disease=factor(sra_tissue$disease,levels=df2[order(df2$samples,decreasing = TRUE),]$disease)

sra_tissue$tissue_main=sra_tissue$tissue
sra_tissue[sra_tissue$tissue=="fibroblast",]$tissue_main="other"
sra_tissue[sra_tissue$tissue=="hip",]$tissue_main="other"
sra_tissue[sra_tissue$tissue=="inner ear",]$tissue_main="other"
sra_tissue[sra_tissue$tissue=="joint",]$tissue_main="other"
sra_tissue[sra_tissue$tissue=="knee",]$tissue_main="other"
sra_tissue[sra_tissue$tissue=="lymph node",]$tissue_main="other"
sra_tissue[sra_tissue$tissue=="muscle",]$tissue_main="other"
sra_tissue[sra_tissue$tissue=="peritoneal lavage",]$tissue_main="other"
sra_tissue[sra_tissue$tissue=="retina",]$tissue_main="other"
sra_tissue[sra_tissue$tissue=="salivary gland",]$tissue_main="other"
sra_tissue[sra_tissue$tissue=="spleen",]$tissue_main="other"




ggplot(data=sra_tissue, aes(x=disease,y=samples,fill=tissue_main)) + geom_bar(position="stack",stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13, colour="black"), 
        axis.title.x = element_blank(), axis.text.y=element_text(colour="black", size = 13), axis.title.y=element_text(colour="black", size = 14)) +
  theme(legend.title = element_blank()) +
  #scale_y_continuous(breaks=seq(0,9000,1000)) +
  scale_fill_brewer(palette = "Paired") + scale_y_continuous(breaks=seq(0,9000,1000)) +
  theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.grid.minor = element_blank(), panel.background = element_blank())

ggplot(data=subset(sra_tissue,tissue_main=="other"), aes(x=disease,y=samples,fill=tissue)) + geom_bar(position="stack",stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13, colour="black"), 
        axis.title.x = element_blank(), axis.text.y=element_text(colour="black", size = 13), axis.title.y=element_text(colour="black", size = 14)) +
  theme(legend.title = element_blank()) +
  #scale_y_continuous(breaks=seq(0,9000,1000)) +
  scale_fill_brewer(palette = "Paired") + scale_y_continuous(breaks=seq(0,250,10)) +
  theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.grid.minor = element_blank(), panel.background = element_blank())


df2=sra_tissue%>%
  group_by(tissue)%>%
  summarise(samples=sum(samples))
#Figure S6----------organism------
sra_organism <- read.delim("/Users/sebastian/pubmed_rheuma_HTS/sra_organism.tsv")
sra_organism=as.data.frame(sra_organism)
sra_organism$organism=as.character(sra_organism$organism)
sra_organism$disease=as.character(sra_organism$disease)
sra_organism$samples=as.integer(sra_organism$samples)
sra_organism[sra_organism$organism=="klebsiella pneumoniae",]$organism="bacteria"
sra_organism[sra_organism$organism=="pseudomonas aeruginosa",]$organism="bacteria"
sra_organism[sra_organism$organism=="staphylococcus aureus",]$organism="bacteria"
sra_organism[sra_organism$organism=="staphylococcus pseudintermedius",]$organism="bacteria"
sra_organism[sra_organism$organism=="streptococcus pyogenes",]$organism="bacteria"
sra_organism[sra_organism$organism=="uncultured bacterium",]$organism="bacteria"
sra_organism[sra_organism$organism=="flavobacterium psychrophilum",]$organism="bacteria"
sra_organism[sra_organism$disease=="MyoPolyDerma",]$disease="(poly/derma)myositis"
sra_organism[sra_organism$disease=="SysSclerosis",]$disease="Systemic Sclerosis"
sra_organism[sra_organism$disease=="Sjoegren",]$disease="Sjögren's Syndrome"
df2=sra_organism%>%
  group_by(disease)%>%
  summarise(samples=sum(samples))
sra_organism$disease=factor(sra_organism$disease,levels=df2[order(df2$samples,decreasing = TRUE),]$disease)
ggplot(data=sra_organism, aes(x=disease,y=samples,fill=organism)) + geom_bar(position="stack",stat="identity") +
  geom_text(data = df2, aes(x=disease,y=9000,label=samples,angle=60),inherit.aes = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13, colour="black"), 
        axis.title.x = element_blank(), axis.text.y=element_text(colour="black", size = 13), axis.title.y=element_text(colour="black", size = 14)) +
  theme(legend.title = element_blank()) +
  scale_fill_brewer(palette = "Paired") + scale_y_continuous(breaks=seq(0,9000,1000)) +
  theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.grid.minor = element_blank(), panel.background = element_blank())

#Figure S7--------phenotype--------
sra_phenotype <- read.delim("/Users/sebastian/pubmed_rheuma_HTS/sra_phenotype.tsv")
sra_phenotype=as.data.frame(sra_phenotype)
sra_phenotype$phenotype=as.character(sra_phenotype$phenotype)
sra_phenotype$disease=as.character(sra_phenotype$disease)
sra_phenotype$samples=as.integer(sra_phenotype$samples)
sra_phenotype[sra_phenotype$disease=="MyoPolyDerma",]$disease="(poly/derma)myositis"
sra_phenotype[sra_phenotype$disease=="SysSclerosis",]$disease="Systemic Sclerosis"
sra_phenotype[sra_phenotype$disease=="Sjoegren",]$disease="Sjögren's Syndrome"
df2=sra_phenotype%>%
  group_by(disease)%>%
  summarise(samples=sum(samples))
sra_phenotype$disease=factor(sra_phenotype$disease,levels=df2[order(df2$samples,decreasing = TRUE),]$disease)
ggplot(data=sra_phenotype, aes(x=disease,y=samples,fill=phenotype)) + geom_bar(position="stack",stat="identity") +
  geom_text(data = df2, aes(x=disease,y=9000,label=samples,angle=60),inherit.aes = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13, colour="black"), 
        axis.title.x = element_blank(), axis.text.y=element_text(colour="black", size = 13), axis.title.y=element_text(colour="black", size = 14)) +
  theme(legend.title = element_blank()) +
  scale_fill_brewer(palette = "Paired") + scale_y_continuous(breaks=seq(0,9000,1000)) +
  theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.grid.minor = element_blank(), panel.background = element_blank())

#------sra vs pubmed------
sra_pubmed <- read.delim("~/pubmed_rheuma_HTS/sra_pubmed_copy.tsv")
sra_pubmed=as.data.frame(sra_pubmed)
sra_pubmed$assay=as.character(sra_pubmed$assay)
sra_pubmed$disease=as.character(sra_pubmed$disease)
sra_pubmed$patient.information=as.character(sra_pubmed$patient.information)
SLE_pubmed=sra_pubmed[sra_pubmed$disease=="SLE",]
SLE_pubmed_rna_patient_info=SLE_pubmed[SLE_pubmed$assay %in% c("rna-seq","ncrna-seq","mirna-seq"),]
SLE_pubmed_rna=SLE_pubmed[SLE_pubmed$assay %in% c("rna-seq","ncrna-seq","mirna-seq"),]$pubmed
SLE_result_df=result_df[result_df$disease=="SLE",]
SLE_result_df_rna=SLE_result_df[SLE_result_df$assay %in% c("RNA-Seq","scRNA-Seq","Targeted DNA Seq,RNA-Seq","Epigenomics,RNA-Seq"),]$pmid
intersect(SLE_result_df_rna,SLE_pubmed_rna)
diff_pubmed_result_df=setdiff(SLE_result_df_rna,SLE_pubmed_rna)
sra_pubmed_missing <- read.delim("~/pubmed_rheuma_HTS/data_availability_sle.txt")
setdiff(diff_pubmed_result_df,sra_pubmed_missing$pmid)
#export results
write.csv2(SLE_pubmed_rna_patient_info,"~/pubmed_rheuma_HTS/SLE_pubmed_rna_patient_info.csv")
#Figure S12-----plot patient information---------------------------
df2=SLE_pubmed_rna_patient_info%>%
  count(patient.information,name="number")
SLE_pubmed_rna_patient_info$patient.information=factor(SLE_pubmed_rna_patient_info$patient.information,levels=df2[order(df2$number,decreasing = TRUE),]$patient.information)
ggplot(data=SLE_pubmed_rna_patient_info, aes(x=patient.information)) + geom_bar() +
  geom_text(data = df2, aes(x=patient.information,y=10,label=number),inherit.aes = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13, colour="black"), 
        axis.title.x = element_blank(), axis.text.y=element_text(colour="black", size = 13), axis.title.y=element_text(colour="black", size = 14)) +
  theme(legend.title = element_blank()) +
  scale_fill_brewer(palette = "Paired") + scale_y_continuous(breaks=seq(0,15,5)) +
  theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.grid.minor = element_blank(), panel.background = element_blank())

#Figure S8-----plot raw data availability---------------------------
df2=sra_pubmed_missing%>%
  #group_by(availability)%>%
  count(availability,name="number")
sra_pubmed_missing$availability=factor(sra_pubmed_missing$availability,levels=df2[order(df2$number,decreasing = TRUE),]$availability)
ggplot(data=sra_pubmed_missing, aes(x=availability)) + geom_bar() +
  geom_text(data = df2, aes(x=availability,y=35,label=number,angle=60),inherit.aes = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13, colour="black"), 
        axis.title.x = element_blank(), axis.text.y=element_text(colour="black", size = 13), axis.title.y=element_text(colour="black", size = 14)) +
  theme(legend.title = element_blank()) +
  scale_fill_brewer(palette = "Paired") + #scale_y_continuous(breaks=seq(0,9000,1000)) +
  theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.grid.minor = element_blank(), panel.background = element_blank())


