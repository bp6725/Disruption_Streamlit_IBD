library(vegan)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(heatmap3)
library(caret)
library(ggpubr)
user="Urigo10"
source(paste0("C:/Users/",user,"/Dropbox/R/functions/filtCols.R"))
source(paste0("C:/Users/",user,"/Dropbox/R/functions/strip_names_silva.R"))
source(paste0("C:/Users/",user,"/Dropbox/R/functions/filt_lefse.R"))
source(paste0("C:/Users/",user,"/Dropbox/R/functions/cor.1dim.r"))
source(paste0("C:/Users/",user,"/Dropbox/R/functions/cor.2dim.r"))
source(paste0("C:/Users/",user,"/Dropbox/R/functions/multiplot.r"))
source(paste0("C:/Users/",user,"/Dropbox/R/functions/change.to.date.R"))

path=paste0("C:/Users/",user,"/Dropbox/IBD/lihi/")###user supplies pathway
path1=paste0("C:/Users/",user,"/Dropbox/IBD/lihi/METABOLON/")###user supplies pathway

#####nice colors to chhose from####
pal=brewer.pal(12,"Set3")
pie(1:12,col=pal)
####read key####
key=read.csv(paste0(path,"key_new.csv"),header=T,check.names=F,comment.char = "",stringsAsFactors = F) ###supply file name
key$Calprotectin=as.numeric(gsub(">","",key$Calprotectin))
key$CRP=as.numeric(gsub(">","",key$CRP))
#remove extra PN10 sampe:
key=key[-grep("S113",key$ID),]
key$date=as.Date(key$date,format="%d/%m/%Y")

met=read.csv(paste0(path,"Metabolon/ScaledImpdata.csv"),header=T,check.names=F,comment.char = "",stringsAsFactors = F) 
#### match metabolon sample to 16s sample: ####
View(met[1:30,1:20])
rownames(met)=met$PATHWAY_SORTORDER

met=met[-1,]
#separate all data decsribing the metabolites to a different object:
met.info=met[ ,1:14]
met=met[ ,-c(1:14)]
identical(rownames(met),rownames(met.info))

met=t(met)
met=as.matrix(met)
mode(met)="numeric"
#### arrange met by key ####
met=met[match(key$ID,rownames(met)),] 
identical(rownames(met),key$ID)
apply(met,2,min)
####set data kind ####
kind= "L6" # "L6" #set here L2 , etc. for otu tables or "weighted_unifrac" or "unweighted_unifrac"
depth="2200"
####This loop reads in dat, format it and matches it to the key, dpending on which 'kind' it is ####
if (kind %in% c("L2","L3","L4","L5","L6","L7","L8","L9","L10","L11","L12")) {
dat=read.table(paste0(path,depth,"/",depth,"_",kind,".txt"),sep="\t",quote="",header=T,check.names=F,skip=1,comment.char = "") ###supply file name
#disregard the # at #OTUID; skip te 1st line (contructed from biom file)

rownames(dat)=dat[ ,1]
dat=dat[ ,-1]
dat=as.matrix(dat)
dat=t(dat)
mode(dat)="numeric"

keep.key=which(key$ID %in% rownames(dat))
keep.dat=which(rownames(dat) %in% key$ID )
key=droplevels(key[keep.key,])

dat=dat[keep.dat,]
dat=dat[match(key$ID,rownames(dat)),]
print(identical(rownames(dat),as.character(key$ID)))   } else {
  dat=read.table(paste0(path,"beta_diversity/",kind,"_",depth,".txt"),sep="\t",quote="",header=T,check.names=F,comment.char = "",stringsAsFactors = F) 
  rownames(dat)=dat[ ,1]
  dat=dat[ ,-1]
  dat=as.matrix(dat)
  mode(dat)="numeric"
  keep.key=which(key$ID %in% rownames(dat))
  keep.datR=which(rownames(dat) %in% key$ID )
  keep.datC=which(colnames(dat) %in% key$ID )
  key=key[keep.key,]
  
  dat=dat[keep.datR,keep.datC]
  dat=dat[match(key$ID,rownames(dat)),match(key$ID,colnames(dat))]
  print(identical(rownames(dat),as.character(key$ID)))
  print(identical(colnames(dat),as.character(key$ID)))
  dis=NULL #(to be used for titles in plots)
}

#### work on met.info AND MAKE MET.AGG####

met.info=as.data.frame(met.info,stringsAsFactors=F)
#to aggrgate met by pathways, 1st match eac compoind to his pathway:
identical(as.character(rownames(met.info)),colnames(met))
temp<-met
temp=t(temp)
temp=as.data.frame(temp)
identical(as.character(rownames(met.info)),rownames(temp))
temp$path<-as.character(met.info$SUB_PATHWAY)
met.agg=aggregate(temp[ ,1:112], by=list(temp$path),FUN=sum,stringsAsFactor=F)
met.agg=t(met.agg)
colnames(met.agg)=met.agg[1,]
met.agg=met.agg[-1,]
mode(met.agg)="numeric"
met.agg=met.agg[match(rownames(dat),rownames(met.agg)),]
identical(rownames(met.agg),rownames(dat))



####Test correlation of SCD, Enterobacter and Ceramides: ####
#manual fix:
grep("Enterobacteriaceae;Other",colnames(dat),ignore.case = T)->x
colnames(dat)[x]=gsub("Other","D_5__Enterobacter",colnames(dat)[x])
####2.	Remove antibiotic  takers: PN04 and 16 ,which are chronic takers, and 26, who took just before washout:####
An=c("PN04","PN16","PN26")
#### All samples





####subset SCD and washout ####
set="Before_After SCD"
key.s=key[which(key$Diet=="Before_SCD"),]  # & !key$PN %in% An ),]#|key$Diet=="Washout"),]
if (kind=="weighted_unifrac"|kind=="unweighted_unifrac") {
  dat.s=dat[match(key.s$ID,rownames(dat)),match(key.s$ID,colnames(dat))] 
  identical(rownames(dat.s),key.s$ID)
  identical(colnames(dat.s),key.s$ID)
}  else {
  dat.s=dat[match(key.s$ID,rownames(dat)),]                               
  identical(rownames(dat.s),key.s$ID) 
}

#### enterobacter and CRP_SCD ####
name1="D_5__Enterobacter"
key.s=key[which(key$Diet=="Before_SCD" & !key$PN %in% An ),]  #|key$Diet=="Before_SCD"),]  # & !key$PN %in% An ),]#|key$Diet=="Washout"),]
if (kind=="weighted_unifrac"|kind=="unweighted_unifrac") {
  dat.s=dat[match(key.s$ID,rownames(dat)),match(key.s$ID,colnames(dat))] 
  identical(rownames(dat.s),key.s$ID)
  identical(colnames(dat.s),key.s$ID)
}  else {
  dat.s=dat[match(key.s$ID,rownames(dat)),]                               
  identical(rownames(dat.s),key.s$ID) 
}

t=grep(name1,colnames(dat.s),ignore.case = T)
colnames(dat.s)[t]
var1=dat.s[ ,t]
identical(names(var1),key.s$ID)
kruskal.test(var1~as.factor(key.s$Increased_CRP_SCD)) #0.19 for AFter SCD, 0.84 for washout
boxplot(var1~as.factor(key.s$Increased_CRP_SCD))
stripchart(var1~as.factor(key.s$Increased_CRP_SCD),add=TRUE,method="jitter",vertical=TRUE,pch=20,col='blue')
tapply(var1,as.factor(key.s$Increased_CRP_SCD),summary)

df=as.data.frame(cbind(var1,key.s$Increased_CRP_SCD))

#caluclate pre/post SCD ratio per PN for enterobacter:
pats=levels(as.factor(key.s$PN))

d=c()
r=c()
for ( i in 1:length(pats)){
  k=key.s[which(key.s$PN==pats[i]),]
  k=k[order(k$date),]
  z=var1[match(k$ID,names(var1))]
  d[i]=z[2]-z[1]
  r[i]=unique(k$Increased_CRP_SCD)
  m=met[match(k$ID,rowname)]
}
df=as.data.frame(cbind(pats,d,r))
df$d=as.numeric(as.character(df$d))
boxplot(df$d~df$r)
kruskal.test(df$d~df$r)
#### lookminto ceramides in SCD ####
name2="Ceramide"
t=grep(name2,met.info$SUB_PATHWAY,ignore.case = T)
met.info$BIOCHEMICAL[t]  
View(met.info[t,])
keep.mets=rownames(met.info)[t]
names.mets=met.info$BIOCHEMICAL[t]
keep.paths=unique(met.info$SUB_PATHWAY[t])
names.paths=unique(met.info$SUB_PATHWAY[t])

met.agg.s=met.agg[match(key.s$ID,rownames(met.agg)), which(colnames(met.agg) %in% keep.paths)]
identical(rownames(met.agg.s),key.s$ID)

met.s=met[match(key.s$ID,rownames(met)), which(colnames(met) %in% keep.mets)]
identical(rownames(met.s),key.s$ID)
allp=c()
for (i in 1:ncol(met.s)){
  var1=met.s[ ,i]
  kruskal.test(var1~as.factor(key.s$Increased_CRP_SCD))->k
  print(paste0("Round ", i, " p.val ",round(k$p.value,2)))
  boxplot(var1~as.factor(key.s$Increased_CRP_SCD),main=paste(set, "Round ",i,names.mets[i]," p ",round(k$p.value,2)),xlab="CRP increased after SCD")
  stripchart(var1~as.factor(key.s$Increased_CRP_SCD), vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue')
   #print(tapply(var1,as.factor(key.s$Increased_CRP_SCD),summary))
  allp[i]<-k$p.value
}
p.adjust(allp,"fdr")
for (i in 1:ncol(met.agg.s)){
  var1=met.agg.s[ ,i]
  
  kruskal.test(var1~as.factor(key.s$Increased_CRP_SCD))->k
  print(paste0("Round ", i, " p.val ",round(k$p.value,2)))
  boxplot(var1~as.factor(key.s$Increased_CRP_SCD),main=paste(set, "Round ",i,names.paths[i]," p ",round(k$p.value,2)),xlab="CRP increased after SCD")
  stripchart(var1~as.factor(key.s$Increased_CRP_SCD), vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue')
  #print(tapply(var1,as.factor(key.s$Increased_CRP_SCD),summary))
  allp[i]<-k$p.value
}

name1="D_5__Enterobacter"
t=grep(name1,colnames(dat.s),ignore.case = T)
colnames(dat.s)[t]
var.t=dat.s[ ,t]
identical(names(var.t),key.s$ID)
pats=levels(as.factor(key.s$PN))

for (j in 1:ncol(met.s)){

d.ent=c()
r=c()
d.met=c()
res=list()
d.CRP=c()
for ( i in 1:length(pats)){
  k=key.s[which(key.s$PN==pats[i]),]
  k=k[order(k$date),]
  var.m=met.s[ ,j]
  name.met=met.info$BIOCHEMICAL[which(rownames(met.info)==colnames(met.s)[j])]
  z=var.m[match(k$ID,names(var.m))]
  d.met[i]=z[2]=z[1]
  
  z=var.t[match(k$ID,names(var.t))]
  d.ent[i]=z[2]-z[1]
  d.CRP[i]=k$CRP[2]-k$CRP[1]
  r[i]=unique(k$Increased_CRP_SCD)
 
}
df=as.data.frame(cbind(pats,r,d.ent,d.met,d.CRP))
df$d.ent=as.numeric(as.character(df$d.ent))
df$d.met=as.numeric(as.character(df$d.met))
df$d.CRP=as.numeric(as.character(df$d.CRP))
kruskal.test(df$d.met~as.factor(df$r))->kp
print(ggplot(df,aes(x=r,y=d.met))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(size=3,position = position_jitter(width=0.05,height=NULL))+
  labs(x="Icreased CRP after SCD",y=name.met,caption=paste("p",round(kp$p.value,2))))
  res[[j]]=df

#cor.test(d.met[which(r==1)],d.ent[which(r==1)],method="spearman")->c1
#cor.test(d.met[which(r==0)],d.ent[which(r==0)],method="spearman")->c2

#cor.test(d.met[which(r==1)],d.CRP[which(r==1)],method="spearman")->c3
#cor.test(d.met[which(r==0)],d.CRP[which(r==0)],method="spearman")->c4

#cor.test(d.ent[which(r==1)],d.CRP[which(r==1)],method="spearman")->c5
#cor.test(d.ent[which(r==0)],d.CRP[which(r==0)],method="spearman")->c6
#cors=list(c1,c2,c3,c4,c5,c6)
#which(lapply(cors, function(x) which(x$p.value<0.1))!=0)->ind
#if (length(ind)>0) {
  #print(paste("Round", j, name.met))
 # print(cors[ind])
  #}
}

cors


for (j in 1:ncol(met.agg.s)){
  
  d.ent=c()
  r=c()
  d.met=c()
  res=list()
  d.CRP=c()
  for ( i in 1:length(pats)){
    k=key.s[which(key.s$PN==pats[i]),]
    k=k[order(k$date),]
    var.m=met.agg.s[ ,j]
    name.met=colnames(met.agg.s)[j]
    z=var.m[match(k$ID,names(var.m))]
    d.met[i]=z[2]=z[1]
    
    z=var.t[match(k$ID,names(var.t))]
    d.ent[i]=z[2]-z[1]
    d.CRP[i]=k$CRP[2]-k$CRP[1]
    if (is.na(d.CRP[i])) {print(paste("Round" , i,"CRP is NA"))} else {
      if (d.CRP[i]<=0) {r[i]<-"No increase"} else {r[i]<-"Increase"}} 
    
  }
  df=as.data.frame(cbind(pats,r,d.ent,d.met,d.CRP))
  df$d.ent=as.numeric(as.character(df$d.ent))
  df$d.met=as.numeric(as.character(df$d.met))
  df$d.CRP=as.numeric(as.character(df$d.CRP))
  kruskal.test(df$d.met~as.factor(df$r))->kp
  print(ggplot(df,aes(x=r,y=d.met))+
          geom_boxplot(outlier.colour = NA)+
          geom_point(size=3,position = position_jitter(width=0.05,height=NULL))+
          labs(x="Icreased CRP after SCD",y=name.met,caption=paste("Round", j, "p",round(kp$p.value,2))))
  
  print(cor.test(d.CRP,d.met,method="spearman"))
  
}

#### drafts ####
cor.test(d.met[which(r==1)],d.ent[which(r==1)],method="spearman")->c1
cor.test(d.met[which(r==0)],d.ent[which(r==0)],method="spearman")->c2

cor.test(d.met[which(r==1)],d.CRP[which(r==1)],method="spearman")->c3
cor.test(d.met[which(r==0)],d.CRP[which(r==0)],method="spearman")->c4

cor.test(d.ent[which(r==1)],d.CRP[which(r==1)],method="spearman")->c5
cor.test(d.ent[which(r==0)],d.CRP[which(r==0)],method="spearman")->c6
cors=list(c1,c2,c3,c4,c5,c6)
which(lapply(cors, function(x) which(x$p.value<0.1))!=0)->ind
if (length(ind)>0) {
  print(paste("Round", j, name.met))
  print(cors[ind])
}
#### end ####

temp.name2="key$Response_CRP_SCD"
name2=unlist(strsplit(temp.name2,"$",fixed=T))[2]
z=summary(as.factor(eval(parse(text = temp.name2)))) #create to label x-axis
var2=eval(parse(text = temp.name2))
kruskal.test(var1,as.factor(var2))->k
p=round(k$p.value,4)
boxplot(var1~as.factor(var2), ylab=name1,main=paste("Mann-Whitney p=",p),cex.axis=0.9,xlab=name2,names=paste0(names(z)," (",z,")"))
stripchart(var1~as.factor(var2), vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue')

tapply(var1,INDEX=var2,FUN=median)->meds
tapply(var1,INDEX=var2,FUN=mean)->means
tapply(var1,INDEX=var2,FUN=max)->maxs
meds
means

jpeg(paste0(path,"Boxplots/",depth,"/No Pouch/",name1,"_",name2,".jpeg"),width=800,height=600)
boxplot(var1~as.factor(var2), ylab=paste(name1,"(Relative Abundance)"),outline=F,main="p (Mann-Whitney) for each group compared to 'Healthy'",cex.main=1.3,xlab=name2,names=paste0(names(z)," (",z,")"),cex.axis=0.8)
stripchart(var1~as.factor(var2), vertical = TRUE,method = "jitter", add = TRUE, pch = 20, col = 'blue')

dev.off()
####boxplot with ggplot (need  dataframe)####
ggplot(key.s,aes(Diagnosis,Shannon,fill=Diagnosis))+
  stat_boxplot(geom ='errorbar') +
  geom_point(size=2)+
  geom_jitter(width=0.05,height=NULL)+
  geom_boxplot(alpha=0.4)+
  scale_fill_manual(values=c(pal[7],pal[5],pal[6],pal[4],pal[8],pal[3]))    +
  labs(y=paste(name1,"(Relative abundance)"),caption="p (all groups):0.1")+
  theme(axis.title.y=element_text(face="italic"),axis.text.x = element_text(angle=70))
ggsave(paste0(path,"Boxplots/",depth,"/",kind,"_",dis," ",name1,"all_groups.jpeg"))



#### run correlations separaly among responders and non-responder####
path2=paste0(path1,"MED/Correlation analysis reponders and non_responders/")
#remove teh antibiotic users - 1 and 2:
key.s=key[which(key$Antibiotics_free_Interval!="1"&key$Antibiotics_free_Interval!="2"),]


dat.s=dat[match(key.s$ID,rownames(dat)),]
identical(rownames(dat.s),key.s$ID)




####divide samples randomly to 2 groups:####
x=key$ID

sample(x,56)->x1
x2=x[-which(x %in% x1)]

setdiff(x1,x2)
key.1=key[which(key$ID %in% x1),]
key.2=key[which(key$ID %in% x2),]

dat.1=dat[match(key.1$ID,rownames(dat)),]
dat.2=dat[match(key.2$ID,rownames(dat)),]
identical(rownames(dat.1),key.1$ID)
identical(rownames(dat.2),key.2$ID)

dat.1=filtCols(dat.1,6)
dat.2=filtCols(dat.2,6)

met.1=met[match(key.1$ID,rownames(met)),]
identical(rownames(dat.1),rownames(met.1))

met.2=met[match(key.2$ID,rownames(met)),]
identical(rownames(dat.2),rownames(met.2))
setdiff(colnames(dat.1),colnames(dat.2))
keep=intersect(colnames(dat.1),colnames(dat.2))
dat.1=dat.1[,which(colnames(dat.1) %in% keep)]
dat.2=dat.2[,which(colnames(dat.2) %in%  keep)]

#### correlate teh subsets of dat and met, tehn compare correlations of subsets####

res1=cor.2dim(met.1,dat.1)
res2=cor.2dim(met.2,dat.2)

t=0.2

length(which(res1$q<t))
res1f=res1[which(res1$q<t),]
#res1f$Var2=met.info$BIOCHEMICAL[match(res1f$Var2,rownames(met.info))]
res1f$interaction=paste(res1f$Var1,"|",res1f$Var2)
res1f$name=met.info$BIOCHEMICAL[match(res1f$Var2,rownames(met.info))]
res1f$Pathway=met.info$SUB_PATHWAY[match(res1f$Var2,rownames(met.info))]
res1f$Super_pathway=met.info$SUPER_PATHWAY[match(res1f$Var2,rownames(met.info))]
hist(res1f$Rho,breaks=60)
hist(res1$q,breaks=60)
hist(res1$p,breaks=60)


length(which(res2$q<t))
res2f=res2[which(res2$q<t),]

hist(res1f$Rho,breaks=60)

res2f$interaction=paste(res2f$Var1,"|",res2f$Var2)
length(setdiff(res1f$interaction,res2f$interaction))

res2f$name=met.info$BIOCHEMICAL[match(res2f$Var2,rownames(met.info))]
res2f$Pathway=met.info$SUB_PATHWAY[match(res2f$Var2,rownames(met.info))]
res2f$Super_pathway=met.info$SUPER_PATHWAY[match(res2f$Var2,rownames(met.info))]

hist(res2f$Rho, col=rgb(0,0,1,0.5),breaks=222,main="Rho of set1 (pink) and set2 (blue)",xlab="Rho")
hist(res1f$Rho, col=rgb(1,0,0,0.5), add=T,breaks=222)


summary(as.factor(res1f$Pathway),maxsum = 1000)->s1
summary(as.factor(res2f$Pathway),maxsum = 1000)->s2
s1=as.data.frame(s1)
s2=as.data.frame(s2)
ms=merge(s1,s2,by="row.names",all=T)
colnames(ms)=c("Pathways","Responders","Non-Responders")
write.csv(ms,paste0(path2,"Metabolites/map_sig_metabolites_to_pathways.csv"))





summary(as.factor(res1f$Super_pathway),maxsum = 1000)->sup1
summary(as.factor(res2f$Super_pathway),maxsum = 1000)->sup2

 
#### make distanmce matrix for both bact and met data: ####
#take only pre/post MED diet:
keep=c()
for (i in 1:length(pats)) {
  x=key[which(key$PN==pats[i]),]
  x=x[order(x$date),]
  k1=which(x$Diet=="MED")
  ids.to.keep=c(x$ID[k1],x$ID[k1-1]) #to account for the 3 pstients with different order
  keep=c(keep,ids.to.keep)
  
}
key.s=key[which(key$ID %in% keep),]
dat.s=dat[match(key.s$ID,rownames(dat)),match(key.s$ID,colnames(dat))]
identical(rownames(dat.s),key.s$ID)
identical(colnames(dat.s),key.s$ID)
met.s=met[match(key.s$ID,rownames(met)),]
identical(rownames(met.s),key.s$ID)
#make bary matrix on met.s:
metdis=as.matrix(vegdist(met.s,"bray"))
bacdis=dat.s
pats=unique(key.s$PN)
#grab metabolomic and bacterial distances per patient before/after diet:
d.bac=c()
d.met=c()

for ( i in 1:length(pats)){
  x=key.s[which(key.s$PN==pats[i]),]
  x=x[order(x$date),]
  d.bac[i]=bacdis[grep(paste0(x$ID[1],"$"),rownames(bacdis)),grep(paste0(x$ID[2],"$"),colnames(bacdis))]
  d.met[i]=metdis[grep(paste0(x$ID[1],"$"),rownames(metdis)),grep(paste0(x$ID[2],"$"),colnames(metdis))]
  
}
plot(d.bac,d.met)

df=as.data.frame(cbind(pats,d.bac,d.met),stringsAsFactors=F)
df$d.bac=as.numeric(df$d.bac)
df$d.met=as.numeric(df$d.met)
resp=c()
for (i in 1:nrow(df)){
  resp[i]=unique(key$Resp_calp_MED[which(key$PN==df$pats[i])])
}
df$resp=resp
df$d.calp=d.calp[match(names(d.calp),df$pats)]
ggplot(df,aes(d.bac,d.met,color=resp))+
  geom_point(size=5)+
  geom_text(aes(label=d.calp),size=3,hjust=-0.4)

cor.test(d.bac,d.met,method="spearman")
cor.test(df$d.met,df$d.bac,method="spearman")


df.f=df[which(df$d.calp>0),]
cor.test(df.f$d.bac,df.f$d.met,method="spearman")


####calculate change in MDI ####
#take only pre/post MED diet:

tax.name="Faecalibacterium"
colnames(dat.s)[grep(tax.name,colnames(dat.s))]
if (length(grep(tax.name,colnames(dat.s)))>1) tax=apply(dat.s[ ,grep(tax.name,colnames(dat.s))],1,sum) else {tax=dat[ ,grep(tax.name,colnames(dat.s))]}


pats=unique(key.s$PN)
d.mdi=c()
d.div=c()
d.calp=c()
d.faecali=c()
for ( i in 1:length(pats)){
  x=key.s[which(key.s$PN==pats[i]),]
  x=x[order(x$date),]
  d.mdi[i]=x$MDI[1]-x$MDI[2]
  d.div[i]=x$Shannon[1]-x$Shannon[2] 
  d.calp[i]=x$Calprotectin[1]-x$Calprotectin[2] 
  d.faecali[i]=tax[which(names(tax)==x$ID[1])]-tax[which(names(tax)==x$ID[2])]
}
names(d.calp)=pats
hist(d.calp,breaks=50)
hist(d.mdi,breaks=50) 
hist(d.div,breaks=50)
sort(d.calp)
cor.test(d.div,d.calp,method="spearman")
plot(d.calp,d.faecali)
#make key for plotting:
dfw=key.s[which(key.s$Diet=="Baseline"),] #to get 1 samp for patient. Note antibiotics ifor refers to baseline sample!

identical(dfw$PN,pats)
dfw$d.bac=d.bac
dfw$d.met=d.met
dfw$d.mdi=d.mdi
dfw$d.div=d.div
# Remove Antibiotic users on same day or previous week:

dfw=dfw[which(dfw$Antibiotics_free_Interval!=1 & dfw$Antibiotics_free_Interval!=2),]
cor.test(dfw$d.bac,dfw$d.met,method="pearson")->ct
p=round(ct$p.value,3)
r=round(ct$estimate,3)

ggsave(paste0(path1,set,"/all_pats show_calprotetin_response Weighted_unifrac.jpeg"))

name="dfw$Patients_Preference"
choose="MED"
name1=unlist(lapply(strsplit(name,"$",fixed=T),"[",2))

dfs=dfw[which(eval(parse(text=name))==choose),]
cor.test(dfs$d.bac,dfs$d.met,method="pearson")->ct
p=round(ct$p.value,3)
r=round(ct$estimate,3)
ggplot(dfs,aes(x=d.bac,y=d.met))+
  geom_point(size=4,color=pal[5])+
  #geom_text(size=4,vjust=1,hjust=1,aes(label=Antibiotics_free_Interval))+
  ggtitle(paste(set, "effects: ",name1,"=",choose,"only"))+
  labs(x="Microbiome shift",y="Metabolome shift",caption=paste("Pearson correlation: p=",p,",R=",r))
  

ggsave(paste0(path1,set,"/",name,"_",choose,".jpeg"))




#make key for plotting:
dfw=key.s[which(key.s$Diet=="Baseline"),] #to get 1 samp for patient. Note antibiotics ifor refers to baseline sample!

identical(dfw$PN,pats)
dfw$d.=d.bac
dfw$d.met=d.met

# Remove Antibiotic users on same day or previous week:

dfw=dfw[which(dfw$Antibiotics_free_Interval!=1 & dfw$Antibiotics_free_Interval!=2),]
dfw$d.mdi


dat.s=filtCols(dat.s,10)
#transfer data to absolute reads:

res=cor.1dim(dat,as.numeric(key$CRP))

dat.s=dat.s*as.numeric(depth)
rowSums(dat.s)
#get rid of zeros:
dat.nz=dat.s+0.2

#find pre/post ratio per patient:
pats=unique(key.s$PN)
fin=matrix(data=NA,length(pats),ncol(dat.nz))
colnames(fin)=colnames(dat.nz)

rownames(fin)=pats
for (i in 1:length(pats)){
  x=key.s[which(key.s$PN==pats[i]),]
  x=x[order(x$date),]
  for (j in 1:ncol(dat.nz)){
    fin[i,j]=dat.nz[grep(x$ID[2],rownames(dat.nz)),j]/dat.nz[grep(paste0("^",eval(x$ID[1]),"$"),rownames(dat.nz)),j]
  }
}

temp=t(fin)






dat.b=ifelse(dat.s>0,1,0) #first "binarize", then find taxa apperaing lerss than half teh samples in each group:
thresholds=c(1:10)
for (q in 1:10) {
tx=thresholds[q]
bef=key.s$ID[grep("Before",key.s$Diet)]
aft=key.s$ID[grep("After",key.s$Diet)]

keep=which(colSums(dat.b[which(rownames(dat.b) %in% bef),])>tx|colSums(dat.b[which(rownames(dat.b) %in% aft),])>tx)    

dat.f=dat.s[ ,keep]
#remove "blank":(last collumn)
dat.f=dat.f[ ,-ncol(dat.f)]
colnames(dat.f)=strip.names.silva(colnames(dat.f),5)
#transfer data to absolute reads:
dat.f=dat.f*as.numeric(depth)
rowSums(dat.f)
#get rid of zeros:
dat.nz=dat.f+0.2
# write a loop calculating the pre/post per raxa per patient:
pats=unique(key.s$PN)
fin=matrix(data=NA,length(pats),ncol(dat.nz))
colnames(fin)=colnames(dat.nz)
rownames(fin)=pats
for (i in 1:length(pats)){
  x=key.s[which(key.s$PN==pats[i]),]
  x=x[order(x$date),]
  for (j in 1:ncol(dat.nz)){
    fin[i,j]=dat.nz[grep(x$ID[2],rownames(dat.nz)),j]/dat.nz[grep(paste0("^",eval(x$ID[1]),"$"),rownames(dat.nz)),j]
  }
}
dis="bray" #define distance method here
if (dis=="bray") {X=F} else {X=T}
  dista=vegdist(fin,method=dis,binary=X)
clus=hclust(dista)
plot(clus)
my_palette <- colorRampPalette(c("green", "black", "red"))(n = 1000)

jpeg(paste0(path,"heatmaps/",set,"_filt",tx,".jpeg"),width=1200,height=1200)
heatmap3(log2(t(fin)),Colv=as.dendrogram(clus),showRowDendro=F,scale="row",margins=c(5,15),col=my_palette,cexRow=1.1,cexCol = 1.1,main=paste(set,"Filter threshold=",tx)) #ColSideColors = x,ColSideLabs = "Group"
dev.off()
}
####make subset of dat and met based on key.s: ####
set="Pre_Post_Med"
key.s=key[grep("Baseline|MED",key$Diet),]






#key.s=key[grep("Washout|SCD",key$Diet),]
if (kind=="weighted_unifrac"|kind=="unweighted_unifrac") {
  dat.s=dat[match(key.s$ID,rownames(dat)),match(key.s$ID,colnames(dat))] 
  identical(rownames(dat.s),key.s$ID)
  identical(colnames(dat.s),key.s$ID)
}  else {
  dat.s=dat[match(key.s$ID,rownames(dat)),]                               
  identical(rownames(dat.s),key.s$ID) 
}

met.s=met[match(key.s$ID,rownames(met)),] 
identical(rownames(met.s),key.s$ID)

met.agg.s=met.agg[match(key.s$ID,rownames(met.agg)),] 
identical(rownames(met.agg.s),key.s$ID)
####correlation mets to bacts: ####
dat.f=filtCols(dat,10)
res=cor.2dim(dat,met)
length(which(res$q<0.05))



#### PCOAs on metabolome data ####


key.s=key[which(key$Antibiotics_free_Interval!="1" &key$Antibiotics_free_Interval!="2"),]
key.s=key.s[which(key.s$Diet=="Baseline"|key.s$Diet=="MED"),]
met.s=met[match(key.s$ID,rownames(met)),] 
identical(rownames(met.s),key.s$ID)


dista=vegdist(met.s,method = "bray")
 
  co=cmdscale(dista,k=5,eig=T)
  eigs=co$eig
  eigs2=eigenvals(co)
  exp=round(100*(abs(eigs))/(sum(abs(eigs))),1) #gives same % explained as "PAST" does
  dfw=as.data.frame(co[[1]])
  identical(rownames(dfw),key.s$ID)
  dfw=merge(dfw,key.s,by.x="row.names",by.y="ID")
  rownames(dfw)=dfw[ ,1]
  dfw=dfw[ ,-1]

  temp.name1=paste0("key.s$Resp_calp_MED")
  name1=unlist(strsplit(temp.name1,"$",fixed=T))[2]
  var1=eval(parse(text = temp.name1))
  var1[which(is.na(var1)==T|var1=="#NULL!")]="Unknown"
  anosim(dista,as.factor(var1),999)->a1
  p1=a1$signif
  r1=round(a1$statistic,2)
 
  
  i=1
  j=2
  ggplot(dfw, aes(x=dfw[ ,i],y=dfw[ ,j],color=eval(parse(text=temp.name1))))+
    #geom_text(size=(3),vjust=-0.5,hjust=-0.5)+
    geom_point(  size=5)+          #aes(size=dat.s[ ,t]))+
    #scale_colour_gradientn(
    #colours = brewer.pal(9,"BrBG")[c(8:9,1:3)])
    #scale_size_continuous(c(2,10))+
    #geom_point(size=dat.s[ ,t])+
    labs(x=paste0("Coordinate ",i," ",exp[i],"%"),y=paste0("Coordinate ",j," ",exp[j],"%"))+
    labs(caption=paste(kind," ", dis," ","p ",name1,"=",p1," R=",r1,sep=" "),size="Candida albicans (RA)",color=name1)+
    # ggtitle("text labels show number of days till ANY KIND of non-antibiotic intervention")+
    theme(plot.title=element_text(size=12,face="italic"), axis.title=element_text(size=12,face="bold"))
  
  ggsave(paste0(path1,set,"/PCOAs/",name1,".jpeg"))
  
####correlate all data and met.agg ####
 dat.f=filtCols(dat,10)
 res=cor.2dim(met.agg,dat)
 res=res[which(res$q<0.01),]
write.csv(res,paste0(path1,"cors1.csv"))
View(met.info[which(met.info$SUB_PATHWAY=="Drug - Antibiotic"),])
####cipro modifiecations ####
cip.mod=c("S88","S32","S16","S63")
dat[which(rownames(dat) %in% cip.mod),grep("Enteroba",colnames(dat ))]
key[which(key$ID %in% cip.mod),]
#### make ratio tables of pre/post MED diet for metabolites/pathways ####
####pre filtration: keep only featuers with high variation.####
thres=min(met.s) #doesnt work - all but one feature ar way above this threshold
filt=c()
min.val=c()
for (i in 1:ncol(met.s)){
  x=met.s[ ,i]
  qs=quantile(x,probs=c(0.25,0.75))
  y=which(x<=qs[1]|x>qs[2])
  y1=which(x==min(x))
  filt[i]=sd(x)/mean(x)
  min.val[i]=min(x)
} #filter metabolites
hist(filt)


filt=c()
for (i in 1:ncol(met.agg.s)){
  x=met.agg.s[ ,i]
  qs=quantile(x,probs=c(0.25,0.75))
  y=which(x<=qs[1]|x>qs[2])
  y1=which(x==min(x))
  filt[i]=sd(x)/mean(x)
  min.val[i]=min(x)
} #filter pathways
min(filt)
hist(filt)


#### Do actaul comparisons of pre/post per patient ratios between responders and non responders ####
#filter object according to filt vector calculated above
met.agg.f=met.agg.s[ ,which(filt>1.2)]
rat=matrix(NA, 28,ncol(met.agg.f))
colnames(rat)=colnames(met.agg.f)
pat=unique(key.s$PN)
for (i in 1:length(pat)){
  k=key.s[which(key.s$PN==pat[i]),]
  k=k[order(k$date),]
  temp=met.agg.f[match(k$ID, rownames(met.agg.f)),]
   for (j in 1:ncol(temp)) {
  rat[ i,j]=temp[ 2,j]/temp[1,j] 
   }
  
}#calculate pre/pot diet ratio per patient per pathway (which passed filtration)


identical(rownames(met.agg.f),key.s$ID)
k.res1=c()
for (i in 1:ncol(rat)){
  wilcox.test(rat[ ,i]~as.factor(response))->k
  k.res1[i]=k$p.value
}
k.q1=p.adjust(k.res1,"fdr")

names(k.q1)<-colnames(rat)
keep=names(k.q1[which(k.q1<0.2)])
df=as.data.frame(rat)
for (i in 1:length(keep)){
  j=grep(keep[i],colnames(df),fixed=T)
  print(ggplot(df,aes(x=as.factor(response),y=df[ ,j]))+
  geom_boxplot(outlier.colour = NA)+
    stat_boxplot(geom="errorbar",width=0.2)+
    labs(y=paste(keep[i],"(Ratio pre/post MED diet)"),x="Response (Calprotectin decrease after MED diet) ",caption=paste("Wilcoxon p:",round(k.res1[grep(keep[i],names(k.res1))],3),"FDR-corrected:",round(k.q1[grep(keep[i],names(k.q1))],2)))+
    geom_jitter(width=0.05))
   # stat_compare_means(method="wilcox.test"))
   ggsave(paste0(path1,"MED/Compare ratios calp_responders vs non_reponders/paths/",keep[i],".jpeg")) 
}



met.f=met.s[ ,which(filt>1.2)]
rat1=matrix(NA, 28,ncol(met.f))
colnames(rat1)=colnames(met.f)
pat=unique(key.s$PN)
for (i in 1:length(pat)){
  k=key.s[which(key.s$PN==pat[i]),]
  k=k[order(k$date),]
  temp=met.f[match(k$ID, rownames(met.f)),]
  for (j in 1:ncol(temp)) {
    rat1[ i,j]=temp[ 2,j]/temp[1,j] 
  }
  
}



#create a vector of response in Calprotectin to MED:
response=key.s$Resp_calp_MED
names(response)=key.s$PN
response=response[unique(names(response))]
k.res1=c()
for (i in 1:ncol(rat1)){
  kruskal.test(rat1[ ,i],as.factor(response))->k
  k.res1[i]=k$p.value
}
names(k.res1)=colnames(rat1)
k.q1=p.adjust(k.res1,"fdr")
sort(k.q1)
names(k.q1)<-colnames(rat1)
keep=names(k.q1[which(k.q1<0.2)])
df=as.data.frame(rat1)
for (i in 1:length(keep)){
  j=grep(keep[i],colnames(df),fixed=T)
  name=as.character(met.info$BIOCHEMICAL[which(rownames(met.info)==keep[i])])
  path=as.character(met.info$SUB_PATHWAY[which(rownames(met.info)==keep[i])])
  print(ggplot(df,aes(x=as.factor(response),y=df[ ,j]))+
          geom_boxplot(outlier.colour = NA)+
          stat_boxplot(geom="errorbar",width=0.2)+
          labs(y=paste(keep[i],"(Ratio pre/post MED diet)"),x="Response (Calprotectin decrease after MED diet) ",caption=paste("Wilcoxon p:",round(k.res1[grep(keep[i],names(k.res1))],3),"FDR-corrected:",round(k.q1[grep(keep[i],names(k.q1))],2)))+
          geom_jitter(width=0.05)+
          ggtitle(paste(name,"; Pathway:",path)))
  # stat_compare_means(method="wilcox.test"))
  ggsave(paste0(path1,"MED/Compare ratios calp_responders vs non_reponders/metabolites/",keep[i],".jpeg")) 
}

#### Test teh BASAL levels of these metabolites in teh 2 classes ####
key.ss=key.s[which(key.s$Diet=="Baseline"),]
base=matrix(NA,nrow=length(pat),ncol=length(keep))
colnames(base)=keep
rownames(base)=pat
for (i in 1:length(pat)){
    k=key.ss$ID[which(key.ss$PN==pat[i])]
     x=met.f[which(rownames(met.f)==k),]
          for ( j in 1:length(keep)){
               base[i,j]=x[which(names(x)==keep[j])]
                    }
  print(is.na(base[ ,1]))
  print(i)
} #a matrix of baseline values of teh significant metabolites
base=as.data.frame(base)
base$Response=key.ss$Resp_calp_MED
p=c()
 for (i in 1:(ncol(base)-1)){
   w=kruskal.test(base[ ,i]~as.factor(base$Response))
   p[i]=w$p.value
  
   name=as.character(met.info$BIOCHEMICAL[which(rownames(met.info)==colnames(base)[i])])
   path=as.character(met.info$SUB_PATHWAY[which(rownames(met.info)==colnames(base)[i])])
   print(ggplot(base,aes(x=as.factor(response),y=base[ ,i]))+
           geom_boxplot(outlier.colour = NA)+
           stat_boxplot(geom="errorbar",width=0.2)+
           labs(y=paste(colnames(base)[i],"(amount at baseline)"),x="Response (Calprotectin decrease after MED diet) ",caption=paste("Wilcoxon p:",round(p[i],3),"FDR-corrected:",round(q[i],2)))+
           geom_jitter(width=0.05)+
           ggtitle(paste(name,"; Pathway:",path)))
   # stat_compare_means(method="wilcox.test"))
   ggsave(paste0(path1,"MED/Compare ratios calp_responders vs non_reponders/metabolites/",keep[i]," at Baseline.jpeg")) 
 }
names(p)=colnames(base)[1:13] 
q=round(p.adjust(p,"fdr"),4)

#### loop over multiple indexes AND variables ####
  #in 1st order patients only:
  
  path=paste0("C:/Users/",user,"/Dropbox/IBD/lihi/response/First_order_only/")
  
distance.method=c("bray","jaccard")
#response=c("Patients_Preference","Response_Calpro","Response_CRP","Response_PDAI","composite_response","Increased_CRP_SCD")
response=c("Resp_calp_MED","Resp_calp_SCD","Response_CRP_MED","Response_CRP_SCD","Any_Calp_response","Patients_Preference","Response_Calpro","Response_CRP","Response_PDAI","composite_response","Increased_CRP_SCD")
list.sum=list()
for (q in 1:2){
####make distance matrix on subset####
dis=distance.method[q] #define distance method here
sum.p=matrix(NA,length(response),5)
colnames(sum.p)=c("Variable","distance","p","R","q")
if (kind %in% c("L2","L3","L4","L5","L6","L7","L8","L9","L10","L11","L12")){
    if (dis=="bray") {X=F} else {X=T}
  dista=vegdist(dat.s,method=dis,binary=X)
} else {dista=dat.s}
for (u in 1:length(response)){
  resp=response[u]
####plot PCOA on subset####
co=cmdscale(dista,k=5,eig=T)
eigs=co$eig
eigs2=eigenvals(co)
exp=round(100*(abs(eigs))/(sum(abs(eigs))),1) #gives same % explained as "PAST" does
dfw=as.data.frame(co[[1]])
identical(rownames(dfw),key.s$ID)

####single variable ####
temp.name1=paste0("key.s$",eval(resp))
name1=unlist(strsplit(temp.name1,"$",fixed=T))[2]
var1=eval(parse(text = temp.name1))
var1[which(is.na(var1)==T|var1=="#NULL!")]="Unknown"
anosim(dista,as.factor(var1),999)->a1
p1=a1$signif
r1=round(a1$statistic,2)
sum.p[u,1]<-name1
sum.p[u,2]<-dis
sum.p[u,3]<-p1
sum.p[u,4]<-r1

i=1
j=2
ggplot(dfw, aes(x=dfw[ ,i],y=dfw[ ,j],color=as.factor(var1)))+
  #geom_text(size=(3),vjust=-0.5,hjust=-0.5)+
  geom_point(  size=5)+          #aes(size=dat.s[ ,t]))+
  #scale_colour_gradientn(
  #colours = brewer.pal(9,"BrBG")[c(8:9,1:3)])
  #scale_size_continuous(c(2,10))+
  #geom_point(size=dat.s[ ,t])+
  labs(x=paste0("Coordinate ",i," ",exp[i],"%"),y=paste0("Coordinate ",j," ",exp[j],"%"),color=name1)+
  labs(caption=paste(kind," ", dis," ","p ",name1,"=",p1," R=",r1,sep=" "),size="Candida albicans (RA)")+
  # ggtitle("text labels show number of days till ANY KIND of non-antibiotic intervention")+
  theme(plot.title=element_text(size=12,face="italic"), axis.title=element_text(size=12,face="bold"))

ggsave(paste0(path,set,"/pcoas/",kind,"_",dis," ",name1,".jpeg"))

}
list.sum[[q]]=sum.p
}
res=do.call("rbind",list.sum)
res[ ,5]=round(p.adjust(res[ ,3],"fdr"),2)
write.csv(res,paste0(path,set,"/",kind,"_AnosimSummary.csv"),row.names = F)




####double variable####

temp.name1="key.s$No.bowel"
name1=unlist(strsplit(temp.name1,"$",fixed=T))[2]
var1=eval(parse(text = temp.name1))

temp.name2="key.s$Diagnosis"
name2=unlist(strsplit(temp.name2,"$",fixed=T))[2]
var2=eval(parse(text = temp.name2))
#turn NAs to uinknow before ANOSIM:
var1[which(is.na(var1)==T|var1=="#NULL!")]="Unknown"
var2[which(is.na(var2)==T|var2=="#NULL!")]="Unknown"
anosim(dista,as.factor(var1),999)->a1
p1=a1$signif
r1=round(a1$statistic,2)

anosim(dista,as.factor(var2),999)->a2
p2=a2$signif
r2=round(a2$statistic,2)
i=1
j=2
ggplot(dfw, aes(x=dfw[ ,i],y=dfw[ ,j],color=as.numeric(var1),shape=as.factor(var2)))+
  # geom_text(size=(4),vjust=-0.5,hjust=-0.5)+
  #geom_point(size=6,aes(color=key.s$CDAI))   +
  #scale_colour_gradientn(colours = brewer.pal(9,"BrBG")[c(8:9,1:3)])
  geom_point(size=5)+
  labs(x=paste0("Coordinate ",i," ",exp[i],"%"),y=paste0("Coordinate ",j," ",exp[j],"%"),color=name2,shape=name1)+
  #scale_color_manual(values=c(pal[5],pal[2]))+
  #labs(caption=paste(kind," ", dis," ","p ",name1,"=",p1," R=",r1,"| p",name2,"=",p2,"R=",r2,sep=" "))+
  # ggtitle("text labels show number of days till ANY KIND of non-antibiotic intervention")+
  theme(plot.title=element_text(size=12,face="italic"), axis.title=element_text(size=12,face="bold"))

ggsave(paste0(path,"PCOAs/CD_healthy_pouch_non_suspIBD",kind,"_",dis," ",name1," ",name2,".jpeg"))

####test associations between two variables, using bosplots and krsukal-wallis ####
temp.name1="key.s$BasToAsc"
name1=unlist(strsplit(temp.name1,"$",fixed=T))[2]
var1=eval(parse(text = temp.name1))

temp.name2="key.s$Diagnosis"
name2=unlist(strsplit(temp.name2,"$",fixed=T))[2]
var2=eval(parse(text = temp.name2))
kruskal.test(var1,as.factor(var2))->k
p=round(k$p.value,2)
boxplot(var1~as.factor(var2), ylab=name1,main=paste("Mann-Whitney p=",p),cex.main=1.3,xlab=name2,ylim=c(0,0))
stripchart(var1~as.factor(var2), vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue')

tapply(var1,INDEX=var2,FUN=median)
jpeg(paste0(path,"/Boxplots/Subsets/PC vs Controls over 50 no pancreatitis/",name1,"_",name2,".jpeg"))
boxplot(var1~as.factor(var2), ylab=name1,main=paste("Mann-Whitney p=",p),cex.main=1.3,xlab=name2)
stripchart(var1~as.factor(var2), vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, col = 'blue')

dev.off()
#### associations of taxon with variable ####

name1="Shannon"
t=grep(name1,colnames(dat.s),ignore.case = T)
var1=dat.s[ ,t]

temp.name2="key.s$Diagnosis"
name2=unlist(strsplit(temp.name2,"$",fixed=T))[2]
z=summary(as.factor(key.s$Diagnosis)) #create to label x-axis
var2=eval(parse(text = temp.name2))
kruskal.test(var1,as.factor(var2))->k
p=round(k$p.value,4)
boxplot(var1~as.factor(var2), ylab=name1,main=paste("Mann-Whitney p=",p),cex.axis=0.9,xlab=name2,names=paste0(names(z)," (",z,")"))
stripchart(var1~as.factor(var2), vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue')

tapply(var1,INDEX=var2,FUN=median)->meds
tapply(var1,INDEX=var2,FUN=mean)->means
tapply(var1,INDEX=var2,FUN=max)->maxs
meds
means

jpeg(paste0(path,"Boxplots/",depth,"/No Pouch/",name1,"_",name2,".jpeg"),width=800,height=600)
boxplot(var1~as.factor(var2), ylab=paste(name1,"(Relative Abundance)"),outline=F,main="p (Mann-Whitney) for each group compared to 'Healthy'",cex.main=1.3,xlab=name2,names=paste0(names(z)," (",z,")"),cex.axis=0.8)
stripchart(var1~as.factor(var2), vertical = TRUE,method = "jitter", add = TRUE, pch = 20, col = 'blue')

dev.off()
####boxplot with ggplot (need  dataframe)####
ggplot(key.s,aes(Diagnosis,Shannon,fill=Diagnosis))+
  stat_boxplot(geom ='errorbar') +
  geom_point(size=2)+
 geom_jitter(width=0.05,height=NULL)+
  geom_boxplot(alpha=0.4)+
  scale_fill_manual(values=c(pal[7],pal[5],pal[6],pal[4],pal[8],pal[3]))    +
  labs(y=paste(name1,"(Relative abundance)"),caption="p (all groups):0.1")+
  theme(axis.title.y=element_text(face="italic"),axis.text.x = element_text(angle=70))
ggsave(paste0(path,"Boxplots/",depth,"/",kind,"_",dis," ",name1,"all_groups.jpeg"))




#for multiple comparisons:
x="BasToAsc"

temp=key[which(key$Diagnosis=="CD"|key$Diagnosis=="Healthy"),] 
kruskal.test(temp[ ,eval(x)],as.factor(temp$Diagnosis))
temp=key[which(key$Diagnosis=="CD"|key$Diagnosis=="Non-IBD"),]
kruskal.test(temp[ ,eval(temp[ ,eval(x)])],as.factor(temp$Diagnosis))
temp=key[which(key$Diagnosis=="CD"|key$Diagnosis=="susp-IBD"),] 
kruskal.test(temp[ ,eval(x)],as.factor(temp$Diagnosis))
temp=key[which(key$Diagnosis=="Healthy"|key$Diagnosis=="susp-IBD"),] 
kruskal.test(temp[ ,eval(x)],as.factor(temp$Diagnosis))
temp=key[which(key$Diagnosis=="Healthy"|key$Diagnosis=="Non-IBD"),]
kruskal.test(temp[ ,eval(x)],as.factor(temp$Diagnosis))
temp=key[which(key$Diagnosis=="CD"|key$Diagnosis=="Pouch"),] 
kruskal.test(temp[ ,eval(x)],as.factor(temp$Diagnosis))
temp=key[which(key$Diagnosis=="Healthy"|key$Diagnosis=="Pouch"),] 
kruskal.test(temp[ ,eval(x)],as.factor(temp$Diagnosis))
temp=key[which(key$Diagnosis=="FAP"|key$Diagnosis=="Pouch"),] 
kruskal.test(temp[ ,eval(x)],as.factor(temp$Diagnosis))
temp=key[which(key$Diagnosis=="Pouch"|key$Diagnosis=="Healthy"),] 
kruskal.test(temp[ ,eval(x)],as.factor(temp$Diagnosis))
####for lfse:####


identical(rownames(dat.s),key.s$ID)
for.lefse=t(dat.s)
key.s$Increased_CRP_SCD=gsub("1","Yes",key.s$Increased_CRP_SCD)
colnames(for.lefse)=key.s$Increased_CRP_SCD

final=filt.lefse(for.lefse,"No","Yes")
final=final[ ,order(colnames(final))]
write.table(final,paste0(path,"response/4_group/",set,"/Response_CRP_SCD/",kind,".txt"),row.names=F,quote=F,sep="\t")
####read in LefSe output ####
lef=read.table(paste0(path,"response/4_group/",set,"/Increased_CRP_SCD/",kind,"1.res"),sep="\t",stringsAsFactors = F)
sigtax=lef$V1[which(is.na(lef$V4)==F)]
####plot lefse bargraph####

p=lef[which(is.na(lef$V4)==F),]
p$V1=gsub("Bacteria","",p$V1)
#remove v2 - NOT SURE WHY ITS THERE!
p=p[ ,-2]
for ( i in 1:nrow(p)){
  if (p$V3[i]=="Control") {p$V4[i]=p$V4[i]*(-1)}
}
p=p[order(p$V4,decreasing = T),]
colnames(p)[3]="LDA"
colnames(p)[2]="Group"
helper=p$V1 #a character vector of taxa names, IN THE ORDER O WANT THEM PLOTTEd
p$V1=factor(p$V1,levels=helper)
jpeg(paste0(path,"lefse.res/L7_PC_vs_Control_over50/bars.jpeg"),width=1000,height=800)
ggplot(p,aes(p$V1,LDA))+
  coord_flip()+
  geom_bar(stat="identity",aes(fill=Group))+
  scale_fill_manual(values=c(pal[7],pal[4]))+
  theme(axis.text.y = element_text(face="italic",size=11))+
  ylab(NULL)
  labs(colour="Group",ylab=NULL)
dev.off()

 
  


####read dat into dl, than change taxa names to fit lefse nomeclature####
dl=dat.s
lev=c(";D_1__", ";D_2__" ,";D_3__" ,";D_4__", ";D_5__" ,";D_6__") 
for (i in 1:6) {
  colnames(dl)=gsub(lev[i],".",colnames(dl),fixed=T)
}
colnames(dl)=gsub("D_0__","",colnames(dl),fixed=T)
colnames(dl)=gsub(" ","",colnames(dl),fixed=T)
dl.s=dl[match(key.s$ID,rownames(dl)),]
identical(rownames(dl.s),key.s$ID)

plotTax=sigtax[which(sigtax %in% colnames(dl))]
res=as.data.frame(dl.s[ ,which(colnames(dl.s) %in% plotTax)])
identical(rownames(res),key.s$ID)
res$ID=rownames(res)
res$response=key.s$Increased_CRP_SCD
res$Status=key.s$Status
res$origin.short=key.s$origin.short
res=res[order(res$response),]
rem=melt(res,id.vars = c("ID","response"))
name=unlist(lapply(strsplit(plotTax,".",fixed=T),"[",6))
#for species:
#name=paste0(name,"_Uncultured bacterium")
for (i in 1:length(plotTax)) {
  x=rem[which(rem$variable==plotTax[i]),]
 x$ID=factor(x$ID,levels=x$ID)
#  x$Type=factor(x$Type,levels=c("Control","Donors","FattyLiver","IPMN","PC"))
  ggplot(x,aes(ID,value,fill=response))+
    geom_bar(stat="identity")+
   # scale_x_discrete(breaks = x$ID, labels =x$origin.short)+
    scale_fill_manual(values=c(pal[12],pal[10]))#,pal[4],pal[5]))+
    labs(y=paste(name[i],"(Relative Abundance)"))+
    theme(axis.title.y=element_text(face="italic"),axis.text.x = element_text(angle=90,vjust=0.3))
    
   # geom_line()+
      #geom_hline(yintercept =median(x$value))+
  ggsave(paste0(path,"response/4_group/",set,"/Increased_CRP_SCD/",name[i],".png"))
}

#draw samw plots on ALL samples:
res=as.data.frame(dl[ ,which(colnames(dl) %in% lef$V1 )])
res$ID=rownames(res)
res$Type=key$Type
res$Status=key$Status
res$origin.short=key$origin.short
res=res[order(res$Status,res$Type,res$origin.short),]
rem=melt(res,id.vars = c("ID","Type","Status","origin.short"))
name=unlist(lapply(strsplit(tax,".",fixed=T),"[",5))
for (i in 1:length(tax)) {
  x=rem[which(rem$variable==tax[i]),]
  x$ID=factor(x$ID,levels=x$ID)
 ggplot(x,aes(ID,value,fill=Type))+
    geom_bar(stat="identity")+
    labs(y=paste(tax[i]))+
   scale_x_discrete(breaks = x$ID, labels =x$origin.short)+
    theme(axis.title.y=element_text(face="italic",size=8),axis.text.x = element_text(angle=70,size=8))+
   theme(axis.ticks.x=element_blank())
      
 ggsave(paste0(path,"lefse.res/L5_PC_vs_Cont/plots_all_samps/",name[i],".jpeg"))
}





ggplot(data=data, aes(x=Category, y=Value, fill=Group))
p + geom_bar(position = 'dodge',stat="identity") +
  geom_text(aes(label=paste(Value, "%")), position=position_dodge(width=0.9),   vjust=-0.25) +
  geom_text(colour="darkgray", aes(y=-3, label=Group),  position=position_dodge(width=0.9), col=gray) +
  theme(legend.position = "none", 
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = "white"),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  annotate("segment", x = 0, xend = Inf, y = 0, yend = 0)

####test associations between two variables, using bosplots and krsukal-wallis ####
temp.name1="key.s$Firmicutes/Bacteroidetes"
name1=unlist(strsplit(temp.name1,"$",fixed=T))[2]
var1=eval(parse(text = temp.name1))

temp.name2="key.s$Status"
name2=unlist(strsplit(temp.name2,"$",fixed=T))[2]
var2=eval(parse(text = temp.name2))
kruskal.test(var1,as.factor(var2))->k
p=round(k$p.value,2)
tapply(var1,INDEX=var2,FUN=median)
jpeg(paste0(path,"/Boxplots/ratio_by_",name2,".jpeg"))
boxplot(var1~as.factor(var2), ylab=name1,main=paste("Mann-Whitney p=",p),cex.main=1.3,xlab=name2)
stripchart(var1~as.factor(var2), vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, col = 'blue')

dev.off()





####Diveristy####
jpeg(paste0(path,"/Boxplots/L5_Div_all_samples.jpeg"))
boxplot(key$div_L5~key$Type, ylab="L5 Shannon index",main="Mann-Whitney p=0.03")
stripchart(key$div_L5~key$Type, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, col = 'blue')
dev.off()
kruskal.test(key$div_L5,as.factor(key$Type))
#key$chao=specpool(dat,pool=rownames(dat))$chao
#key$shannon=diversity(dat,index="shannon",MARGIN=1)

#### make heatmap ####
#remove taxa below 0.005 IN ALL SAMPLES:
m=apply(dat,2,max)
dat.f=dat[ , which(m>0.01)]
dat.f=filtCols(dat.f,30)
colnames(dat.f)=strip.names.silva(colnames(dat.f),5)
clus=hclust(dista)
breaks=c(0,0.005,0.0075,0.015,0.03,0.06,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)  
palette=brewer.pal(9,"Blues")[3:9]

temp=c("white",colorRampPalette(palette)(14))
x=c()
for (i in 1:nrow(dat.f)){
  if (key.s$Group.x[i]=="1") {x[i]="blue"} else {x[i]="violet"}
}
heatmap3(t(dat.f),Colv=as.dendrogram(clus),Rowv=NA,showRowDendro=F,scale="none",breaks=breaks,col=temp,cexRow=0.8,cexCol = 1.1) #ColSideColors = x,ColSideLabs = "Group"

####make distance matrixes (if dat is a unifrac matrix, than 'dista' will be set to dat)####
if (kind %in% c("L2","L3","L4","L5","L6","L7")){
  dis="bray" #define distance method here
  if (dis=="bray") {X=F} else {X=T}
  dista=vegdist(dat,method=dis,binary=X)
} else {dista=dat}
####plot PCOA####

co=cmdscale(dista,k=5,eig=T)
eigs=co$eig
eigs2=eigenvals(co)
exp=round(100*(abs(eigs))/(sum(abs(eigs))),1) #gives same % explained as "PAST" does
dfw=as.data.frame(co[[1]])
identical(rownames(dfw),key$ID)
#set here which coordinates to draw:
i=1
j=2
#order factor levels as wanted:
#key$Status=factor(key$Status,levels=c("PC","IPMN","FattyLiver","Control"))
#key$Type=factor(key$Type,levels=c("PC","IPMN","FattyLiver","Control","Donors"))
#define variables for anosim and plotting:
temp.name1="key$Type"
name1=unlist(strsplit(temp.name1,"$",fixed=T))[2]
var1=eval(parse(text = temp.name1))

temp.name2="key$Status"
name2=unlist(strsplit(temp.name2,"$",fixed=T))[2]
var2=eval(parse(text = temp.name2))
#turn NAs to uinknow before ANOSIM:
#var1[which(is.na(var1)==T|var1=="#NULL!")]="Unknown"
#var2[which(is.na(var2)==T|var2=="#NULL!")]="Unknown"
anosim(dista,as.factor(var2),999)->a
p=a$signif
r=round(a$statistic,2)
ggplot(dfw, aes(x=dfw[ ,i],y=dfw[ ,j],color=as.factor(var2)))+
  # geom_text(size=(4),vjust=-0.5,hjust=-0.5)+
  #geom_point(size=6,aes(color=key.s$CDAI))   +
  #scale_colour_gradientn(colours = brewer.pal(9,"BrBG")[c(8:9,1:3)])
  geom_point(size=5)+
  labs(x=paste0("Coordinate ",i," ",exp[i],"%"),y=paste0("Coordinate ",j," ",exp[j],"%"),color=name2)+
  #scale_color_manual(values=c(pal[5],pal[2]))+
  ggtitle(paste(kind," ", dis," ","p=",p," R=",r,sep=""))+
  # ggtitle("text labels show number of days till ANY KIND of non-antibiotic intervention")+
  theme(plot.title=element_text(size=12,face="italic"), axis.title=element_text(size=12,face="bold"))

ggsave(paste0(path,"PCOAs/",kind,"_",dis," ",name1," ","coord2_by_1.jpeg"))
