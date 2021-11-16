#Recreate DRAM output from liquor

#packages I want to use
library(tidyverse)
#library(dplyr)
#library(tidyr)
#require(devtools)
#install_version("ggplot2", version = "3.3.3", repos = "http://cran.us.r-project.org")
#library(ggplot2)
library(gridExtra)
theme_set(theme_classic());theme_update(axis.text = element_text(color="black"))

#Working in
#setwd("~/Downloads")

#Read recruitment data (note that if using Libre Calc, saving as a csv actually saves as tsv)
Liquor1=as.data.frame(read.delim(file="DRAM/genome_summaries/product.tsv",header=T)) %>% 
  mutate(origin="UGent",modifications="None")
Liquor2=as.data.frame(read.delim(file="Mod_MAGs/DRAM/genome_summaries/product.tsv",header=T)) %>% 
  mutate(origin="UGent",modifications="Filt")
Liquor3=as.data.frame(read.delim(file="Brocad_refs/DRAM/genome_summaries/product.tsv",header=T)) %>% 
  mutate(origin="GTDBBrocadfromNCBI",modifications="None")
Liquor=rbind(Liquor1,Liquor2,Liquor3)
Liquor_Key=as.data.frame(read.delim(file="DramKey.txt",header=T))
Liquor_Key2=Liquor_Key %>% 
  mutate(Path=factor(Path,levels=unique(Path))) %>% 
  mutate(Path2=factor(Path2,levels=unique(Path2))) %>% 
  mutate(Component=factor(Component,levels=unique(Component))) %>% 
  mutate(Component2=factor(Component2,levels=unique(Component2)))
Genomes1=as.data.frame(read.delim(file="DRAM/genome_summaries/genome_stats.tsv",header=T)) %>% 
  mutate(origin="UGent",modifications="None")
Genomes2=as.data.frame(read.delim(file="Mod_MAGs/DRAM/genome_summaries/genome_stats.tsv",header=T)) %>% 
  mutate(origin="UGent",modifications="Filt")
Genomes3=as.data.frame(read.delim(file="Brocad_refs/DRAM/genome_summaries/genome_stats.tsv",header=T)) %>% 
  mutate(origin="GTDBBrocadfromNCBI",modifications="None")
Genomes=rbind(Genomes1,Genomes2,Genomes3)

Fulldat=full_join(Genomes %>% 
                    mutate(taxonomy=ifelse(grepl("__",taxonomy),as.character(taxonomy),"None"),
                           taxonomy=factor(taxonomy,levels=unique(taxonomy))) %>% 
                    separate(taxonomy,c("GTDB_Dom","GTDB_Phy","GTDB_Cla","GTDB_Ord","GTDB_Fam","GTDB_Gen","GTDB_Spe"),sep=";",remove=F) %>% 
                    arrange(origin,GTDB_Dom,GTDB_Phy,GTDB_Cla,GTDB_Ord,GTDB_Fam,GTDB_Gen) %>% 
                    mutate(genome2=factor(genome,levels=unique(genome)),
                           taxonomy2=gsub("d__.*p__","p__",taxonomy),
                           taxonomy2=gsub(";s__.*","",taxonomy2),
                           taxonomy3=gsub(".*;o__","o__",taxonomy),
                           taxonomy3=gsub(";s__.*","",taxonomy3),
                           taxonomy4=gsub(".*f__","f__",taxonomy),
                           taxonomy=gsub(";s__.*","",taxonomy),
                           taxonomy2=factor(taxonomy2,levels=unique(taxonomy2)),
                           taxonomy3=factor(taxonomy3,levels=unique(taxonomy3)),
                           taxonomy4=factor(taxonomy4,levels=unique(taxonomy4))), 
                  Liquor) 

#Massage data a bit
Liquor_ETC=Fulldat %>% 
  select(genome,genome2,origin,modifications,taxonomy,taxonomy2,taxonomy3,taxonomy4,GTDB_Phy:GTDB_Spe,X3.Hydroxypropionate.bi.cycle:Complex.V..V.A.type.ATPase..prokaryotes)
Liquor_ETC2=Liquor_ETC %>% 
  group_by(genome) %>% 
  gather("Path.Component","Completion",X3.Hydroxypropionate.bi.cycle:Complex.V..V.A.type.ATPase..prokaryotes)
Liquor_Key2_ETC=Liquor_Key2 %>% 
  filter(Process=="ETC") %>% 
  select(-Path_Order:-Orig_Order)
Liquor_ETC3=merge(Liquor_Key2_ETC,Liquor_ETC2,all.y = T)

Liquor_Other=Fulldat %>% 
  select(genome,genome2,origin,modifications,taxonomy,taxonomy2,taxonomy3,taxonomy4,GTDB_Phy:GTDB_Spe,CAZy..Alpha.galactans:Sulfur.metabolism..thiosulfate....sulfite)
Liquor_Other2=Liquor_Other %>%  
  group_by(genome) %>% 
  gather("Path.Component","Presence",CAZy..Alpha.galactans:Sulfur.metabolism..thiosulfate....sulfite)
Liquor_Key2_Other=Liquor_Key2 %>% 
  filter(Process=="Other") %>% 
  select(-Path_Order:-Orig_Order)
Liquor_Other3=merge(Liquor_Key2_Other,Liquor_Other2,all.y = T)

#Visualize metabolisms, first recreate, then try to improve
DRAM_ETC_Recreation=Liquor_ETC3 %>% 
  filter(origin=="UGent",modifications=="None") %>% 
  ggplot()+geom_tile(aes(x=Component,y=genome2,fill=Completion))+scale_fill_distiller(palette="YlGnBu",direction=1)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0),aspect.ratio = 1)+facet_grid(~Path2,scales="free",space="free")+xlab("")+ylab("")+ggtitle("ETC Complexes")
DRAM_ETC_Improved=Liquor_ETC3 %>% 
  filter(origin=="UGent") %>% 
  ggplot()+
  geom_tile(aes(x=Component2,y=genome2,fill=ifelse(Completion>0,Completion,NA)),size=0.1,color="black")+
  scale_fill_distiller("Component completion",palette="YlGnBu",direction=1,na.value = "white")+
  scale_y_discrete(breaks=Fulldat$genome2,labels=ifelse(Fulldat$taxonomy3=="None",as.character(Fulldat$genome2),as.character(Fulldat$taxonomy3)))+
  facet_grid(modifications~Path2,scales="free",space="free")+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0),aspect.ratio = 1)+xlab("")+ylab("")+ggtitle("Energy conservation pathways and ETC Complexes")+guides(color=F)
#ggsave("DRAM_ETCimp.pdf",plot=DRAM_ETC_Improved,height=20,width=30,units="cm",useDingbats=F)

DRAM_Other_Recreation=Liquor_Other3 %>% 
  ggplot()+geom_tile(aes(x=Component,y=genome2,fill=Presence))+scale_fill_manual(values=c("lightcyan","mediumseagreen"))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0),aspect.ratio = 1)+facet_grid(~Path2,scales="free",space="free")+xlab("")+ylab("")
DRAM_Other_Improved=Liquor_Other3 %>%
  filter(origin=="UGent") %>% 
  ggplot()+
  geom_tile(aes(x=Component2,y=genome2,fill=ifelse(Presence=="True",as.character(Path2),NA)),size=0.1,color="black")+
  scale_fill_brewer("Function presence",palette="YlGnBu",direction=-1)+
  scale_y_discrete(breaks=Fulldat$genome2,labels=ifelse(Fulldat$taxonomy3=="None",as.character(Fulldat$genome2),as.character(Fulldat$taxonomy3)))+
  facet_grid(modifications~Path2,scales="free",space="free")+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0),aspect.ratio = 1)+xlab("")+ylab("")+ggtitle("Dissimilatory (and some assimilatory) pathways")+guides(fill=F)
#ggsave("DRAM_Otherimp.pdf",plot=DRAM_Other_Improved,height=20,width=40,units="cm",useDingbats=F)

DRAM_Qual=Fulldat %>% 
  filter(origin=="UGent") %>% 
  ggplot()+
  geom_bar(aes(x=completeness.score,y=genome2,fill=GTDB_Phy,group=modifications),stat="identity",color="black")+
  geom_bar(aes(x=contamination.score,y=genome2,group=modifications),stat="identity",fill="gray50",color="black")+
  geom_text(aes(x=80,y=genome2,label=ifelse(X16S.rRNA!="","16S",NA),group=modifications),vjust=0.5)+
  geom_text(aes(x=95,y=genome2,label=assembly.quality,group=modifications),vjust=0.5)+
  geom_path(aes(x=sqrt(number.of.scaffolds),y=genome2,group=42),color="black",size=1)+
  geom_point(aes(x=tRNA.count,y=genome2,group=modifications),color="black",fill="gray70",shape=21,size=2.5)+
  scale_fill_brewer("Phylum",palette="Set1",na.value="brown")+
  coord_cartesian(xlim=c(0,100))+
  scale_y_discrete(breaks=Fulldat$genome2,labels=ifelse(Fulldat$taxonomy3=="None",as.character(Fulldat$genome2),as.character(Fulldat$taxonomy3)))+
  facet_grid(modifications~.,scales="free_y",space="free")+
  xlab("Checkm (bars, %)\ntRNA count (points)\nContigs (line, square root)")+ylab("")+ggtitle("MAG qualities")
#ggsave("DRAM_qual.pdf",plot=DRAM_Qual,height=20,width=30,units="cm",useDingbats=F)

DRAM_sum=grid.arrange(DRAM_Qual+theme(axis.title.x=element_text(angle=90,hjust=1,vjust=0.5)),
                      DRAM_ETC_Improved+theme(axis.ticks.y=element_blank(),axis.text.y=element_blank(),aspect.ratio = 1),
                      DRAM_Other_Improved+theme(axis.ticks.y=element_blank(),axis.text.y=element_blank(),aspect.ratio = 1,legend.position="none"),
                      nrow=1,widths=c(1,1.3,1.7))
#ggsave("DramSmry.pdf",plot=DRAM_sum,height=20,width=90,units="cm",useDingbats=F)


Broc_etc=Liquor_ETC3 %>% 
  filter(grepl("c__Brocadiae",taxonomy)) %>%
  arrange(GTDB_Ord,GTDB_Fam,GTDB_Gen,GTDB_Spe,taxonomy,genome2) %>% 
  mutate(taxonomy4=factor(taxonomy4,levels=unique(taxonomy4)),
         genome2=factor(genome2,levels=unique(genome2))) %>% 
  ggplot()+
  geom_tile(aes(x=Component2,y=genome2,fill=ifelse(Completion>0,Completion,NA),color=interaction(origin,modifications)),size=0.1)+
  scale_fill_distiller("Component\ncompletion",palette="YlGnBu",direction=1,na.value = "white")+
  scale_color_manual(values=c("blue","gray50","red"))+
  scale_y_discrete(breaks=Fulldat$genome2,labels=Fulldat$taxonomy4)+
  facet_grid(~Path2,scales="free",space="free")+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0),aspect.ratio = 1)+xlab("")+ylab("")+ggtitle("Energy conservation pathways and ETC Complexes")+guides(color=F)
Broc_other=Liquor_Other3 %>%
  filter(grepl("c__Brocad",taxonomy)) %>%
  arrange(GTDB_Ord,GTDB_Fam,GTDB_Gen,GTDB_Spe,taxonomy,genome2) %>% 
  mutate(taxonomy4=factor(taxonomy4,levels=unique(taxonomy4)),
         genome2=factor(genome2,levels=unique(genome2))) %>% 
  ggplot()+
  geom_tile(aes(x=Component2,y=genome2,fill=ifelse(Presence=="True",as.character(Path2),NA),color=interaction(origin,modifications)),size=0.1)+
  geom_point(aes(x=ifelse(Component2=="nitrite=>nitric oxide",as.character(Component2),NA),y=genome2),size=0.5,color="gray50")+
  scale_fill_brewer("Function presence",palette="YlGnBu",direction=-1)+
  scale_color_manual("",values=c("blue","gray50","red"))+
  scale_y_discrete(breaks=Fulldat$genome2,labels=Fulldat$taxonomy4)+
  facet_grid(~Path2,scales="free",space="free")+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0),aspect.ratio = 1)+xlab("")+ylab("")+ggtitle("Dissimilatory (and some assimilatory) pathways")+guides(fill=F,color=F)
Broc_qual=Fulldat %>%
  filter(grepl("c__Brocadiae",taxonomy)) %>%
  arrange(GTDB_Ord,GTDB_Fam,GTDB_Gen,GTDB_Spe,taxonomy,genome2) %>% 
  mutate(taxonomy4=factor(taxonomy4,levels=unique(taxonomy4)),
         genome2=factor(genome2,levels=unique(genome2))) %>% 
  ggplot()+
  geom_bar(aes(x=completeness.score,y=genome2,fill=interaction(origin,modifications)),stat="identity",color="black")+
  geom_bar(aes(x=contamination.score,y=genome2),stat="identity",fill="gray30",color="black")+
  geom_text(aes(x=80,y=genome2,label=ifelse(X16S.rRNA!="","16S",NA)),vjust=0.5)+
  geom_text(aes(x=95,y=genome2,label=assembly.quality),vjust=0.5)+
  geom_path(aes(x=sqrt(number.of.scaffolds),y=genome2,group=42),color="black",size=1)+
  geom_point(aes(x=tRNA.count,y=genome2),color="black",fill="gray70",shape=21,size=2.5)+
  scale_fill_manual(values=c("blue","gray50","red"))+theme(legend.position="none")+
  coord_cartesian(xlim=c(0,100))+
  scale_y_discrete(breaks=Fulldat$genome2,labels=Fulldat$taxonomy4)+
  xlab("Checkm (bars, %)\ntRNA count (points)\nContigs (line, square root)")+ylab("")+ggtitle("MAG qualities")+
  facet_grid(~"Various metrics",scales="free",space="free")

Broc_sum=grid.arrange(Broc_qual+theme(axis.title.x=element_text(angle=90,hjust=1,vjust=0.5)),
                      Broc_etc+theme(axis.ticks.y=element_blank(),axis.text.y=element_blank(),aspect.ratio = 1),
                      Broc_other+theme(axis.ticks.y=element_blank(),axis.text.y=element_blank(),aspect.ratio = 1,legend.position="none"),
                      nrow=1,widths=c(1,0.65,0.85))
#ggsave("Brocadia_UG_ETC.pdf",plot=Broc_etc,height=40,width=30,units="cm",useDingbats=F)
#ggsave("Brocadia_UG_other.pdf",plot=Broc_other,height=40,width=40,units="cm",useDingbats=F)
#ggsave("Brocadia_UG_Qual.pdf",plot=Broc_qual,height=40,width=30,units="cm",useDingbats=F)
#ggsave("Brocadia_UG_sum.pdf",plot=Broc_sum,height=30,width=75,units="cm",useDingbats=F)

