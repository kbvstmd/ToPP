rm(list = ls(all = TRUE))

suppressMessages(library(rms))
suppressMessages(library(survminer))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))

setwd("data")



args = commandArgs()
ca = args[(grep("--ca",args)+1)]
tm = args[(grep("--tm",args)+1)]
variable1 = args[(grep("--variable1",args)+1)]
variable2 = args[(grep("--variable2",args)+1)]
timic = args[(grep("--timic",args)+1)]
conf.int = args[(grep("--confint",args)+1)]
risk.table = args[(grep("--risktable",args)+1)]
groupc1 = args[(grep("--groupc1",args)+1)]
groupc2 = args[(grep("--groupc2",args)+1)]
groupc3 = args[(grep("--groupc3",args)+1)]
groupc4 = args[(grep("--groupc4",args)+1)]
surv.median.line = args[(grep("--survmedianline",args)+1)]
survivaltype = args[(grep("--survivaltype",args)+1)]
datatype1 = args[(grep("--datatype1",args)+1)]
datatype2 = args[(grep("--datatype2",args)+1)]
cutoff = args[(grep("--cutoff",args)+1)]




# ca = "TCGA-LGG;TCGA-LIHC"
# tm = '22020200202'
# variable1 = 'TP53'
# variable2 = 'IDH1'
# timic = 'years'
# conf.int = FALSE
# risk.table = TRUE
# groupc1 = 'red'
# groupc2='green'
# groupc3 = 'blue'
# groupc4='pink'
# surv.median.line = 'hv'
# survivaltype = 'OS'
# datatype1 = 'Gene_expression'
# datatype2 = 'Mutation'
# cutoff = 'median'

print(variable2)

title=ca
conf.int=eval(parse(text = conf.int))
risk.table=eval(parse(text = risk.table))


title=ca
cas=str_split(ca,';')[[1]]
time=paste(survivaltype,'.time',sep='')
event=survivaltype
select1=c('sample','group',variable1,event,time)
select2=c('sample','group',variable2,event,time)

print(c(select1,variable2))

dataall1=NULL
dataall2=NULL
for (i in 1:length(cas)) {
  ca=cas[i]
  tmpdata1=fread(paste('rdata/pancancer/',datatype1,'/',ca,'.txt',sep=''),stringsAsFactors = F,select = select1)
  tmpdata1=tmpdata1[which(tmpdata1$group==ca),]
  dataall1=rbind(dataall1,tmpdata1)
  tmpdata2=fread(paste('rdata/pancancer/',datatype2,'/',ca,'.txt',sep=''),stringsAsFactors = F,select = select2)
  tmpdata2=tmpdata2[which(tmpdata2$group==ca),]
  dataall2=rbind(dataall2,tmpdata2)
}

delna=function(data,time,event){
  deltime=which(is.na(data[time]))
  if (length(deltime)>0) {
    data=data[-deltime,]
  }
  delevent=which(is.na(data[event]))
  if (length(delevent)>0) {
    data=data[-delevent,]
  }
  data[is.na(data)]<-0
  return(data.frame(data))
}



dataall2=data.frame(dataall2,stringsAsFactors = F)
dataall1=data.frame(dataall1,stringsAsFactors = F)

dataall1=delna(dataall1,time,event)
dataall2=delna(dataall2,time,event)



if(datatype1!="Mutation"&datatype1!="CNA"&datatype1!="fusion"&datatype1!="clinical"&datatype1!="Methylation"){
  
  if(cutoff=="Median"){
    dataall1[,variable1]=ifelse(dataall1[,variable1]<=summary(dataall1[,variable1])[3],"low","high")
    dataall1[,variable1]=factor(dataall1[,variable1],levels =c("low","high"))
    
  }else if(cutoff=="bestcut"){
    res.cut=surv_cutpoint(data, time = time, event = event, variable1,minprop = 0.1, progressbar = TRUE)
    res.point=res.cut$cutpoint$cutpoint
    dataall1[,variable1]=ifelse(dataall1[,variable1]<=res.point,"low","high")
    dataall1[,variable1]=factor(dataall1[,variable1],levels =c("low","high"))
  }else{
    dataall1[,variable1]=ifelse(dataall1[,variable1]<=summary(dataall1[,variable1])[3],"low","high")
    dataall1[,variable1]=factor(dataall1[,variable1],levels =c("low","high"))
  }
  
}else if(datatype1=="CNA"){
  data=dataall1[-which(dataall1[,variable1]==0),]
  dataall1[,variable1]=ifelse(dataall1[,variable1]<0,"CNV loss","CNV gain")
  dataall1[,variable1]=factor(dataall1[,variable1],levels =c("CNV loss","CNV gain"))
}else if(datatype1=="fusion"){
  
  dataall1[,variable1]=ifelse(dataall1[,variable1]==1,"fusion","Wildtype")
  dataall1[,variable1]=factor(dataall1[,variable1],levels =c("fusion","Wildtype"))
  
}else if(datatype1=="Mutation"){
  dataall1[,variable1]=ifelse(dataall1[,variable1]==1,"Mutant","Wildtype")
  dataall1[,variable1]=factor(dataall1[,variable1],levels =c("Mutant","Wildtype"))
}else if(datatype1=="Methylation"){
  dataall1=dataall1[-which(dataall1[,variable1]>0.2&dataall1[,variable1]<0.8),]
  dataall1[,variable1]=ifelse(dataall1[,variable1]>=0.8,"Hypermethylation","Hypomethylation")
  dataall1[,variable1]=factor(dataall1[,variable1],levels =c("Hypermethylation","Hypomethylation"))
}


if(datatype2!="Mutation"&datatype2!="CNA"&datatype2!="fusion"&datatype2!="clinical"&datatype1!="Methylation"){
  
  if(cutoff=="Median"){
    dataall2[,variable2]=ifelse(dataall2[,variable2]<=summary(dataall2[,variable2])[3],"low","high")
    dataall2[,variable2]=factor(dataall2[,variable2],levels =c("low","high"))
    
  }else if(cutoff=="bestcut"){
    res.cut=surv_cutpoint(data, time = time, event = event, variable2,minprop = 0.1, progressbar = TRUE)
    res.point=res.cut$cutpoint$cutpoint
    dataall2[,variable2]=ifelse(dataall2[,variable2]<=res.point,"low","high")
    dataall2[,variable2]=factor(dataall2[,variable2],levels =c("low","high"))
  }else{
    dataall2[,variable2]=ifelse(dataall2[,variable2]<=summary(dataall2[,variable2])[3],"low","high")
    dataall2[,variable2]=factor(dataall2[,variable2],levels =c("low","high"))
  }
  
}else if(datatype2=="CNA"){
  data=dataall2[-which(dataall2[,variable2]==0),]
  dataall2[,variable2]=ifelse(dataall2[,variable2]<0,"loss","gain")
  dataall2[,variable2]=factor(dataall2[,variable2],levels =c("loss","gain"))
}else if(datatype2=="fusion"){
  
  dataall2[,variable2]=ifelse(dataall2[,variable2]==1,"fusion","Wildtype")
  dataall2[,variable2]=factor(dataall2[,variable2],levels =c("fusion","Wildtype"))
  
}else if(datatype2=="Mutation"){
  dataall2[,variable2]=ifelse(dataall2[,variable2]==1,"Mutant","Wildtype")
  dataall2[,variable2]=factor(dataall2[,variable2],levels =c("Mutant","Wildtype"))
}else if(datatype2=="Methylation"){
  dataall2=dataall2[-which(dataall2[,variable2]>0.2&dataall2[,variable2]<0.8),]
  dataall2[,variable2]=ifelse(dataall2[,variable2]>=0.8,"Hypermethylation","Hypomethylation")
  dataall2[,variable2]=factor(dataall2[,variable2],levels =c("Hypermethylation","Hypomethylation"))
}



print(c(select1,variable2))

dataall=merge(dataall1,dataall2[,c('sample',variable2)],by='sample')
if(variable1==variable2){
  variable2=paste(variable2,'.y',sep='')
  dataall=dataall[,c(select1,variable2)]
  data=delna(dataall,time,event)
  colnames(data)[2]='catype'
  data['group']=paste(variable1,data[,variable1],"+",variable1,data[,variable2],sep=' ')
}else{
  dataall=dataall[,c(select1,variable2)]
  data=delna(dataall,time,event)
  colnames(data)[2]='catype'
  data['group']=paste(variable1,data[,variable1],"+",variable2,data[,variable2],sep=' ')
}

if(timic=="years")
{data[,time]=data[,time]/365}else if(timic=="months")
{data[,time]=data[,time]/30}


print(data[1:3,1:2])
variable='group'


f.np = npsurv(formula =as.formula(paste("Surv(time = as.numeric(data[,time]),event= as.numeric(data[,event])) ~ factor(`",variable,"`)",sep="")),data=data)

gp=data.frame(table(data$group))

colors=c("#DC143C","#4169E1","#CD853F","#3CB371","#808080","#FFC0CB","#00008B","#BDB76B","#EE82EE","#48D1CC","#9932CC","#4682B4","#3CB371","#D2B48C","#B22222")

if(length(gp$Var1)==4){
  palettecol=c(groupc1,groupc2,groupc3,groupc4)
}else{
  palettecol=colors[1:length(gp$Var1)]
}


print(surv.median.line)
ggsurv=ggsurvplot( fit=f.np,
                   data = data,
                   main = "Survival curve",
                   title=title,
                   #break.time.by=0.2*max(data$stime),#XÖácut
                   size = 1.1, # change line size
                   xlab = paste("Time in ",'years',sep = ""),# customize X axis label.
                   conf.int = conf.int, # Add confidence interval
                   pval = TRUE, # Add p-value
                   risk.table = risk.table, # Add risk table
                   risk.table.col = "strata",# Risk table color by groups
                   surv.median.line = surv.median.line, # add the median survival pointer.
                   legend= c(0.8,0.92),
                   legend.title = "",
                   palette = palettecol,
                   legend.labs = as.character(gp$Var1), # Change legend labels
                   risk.table.y.text = FALSE, # hide risktabl y
                   risk.table.height = 0.25, # Useful to change when you have multiple groups
                   risk.table.fontsize=5,
                   font.tickslab = c(12, "bold"),#km×ø±êÖá
                   font.legend=c(12, "bold"),#legend´óÐ¡ÑÕÉ«
                   #font.x = c(5, "bold", "red"),xÖá±êÇ©´óÐ¡
                   conf.int.style = "step", # customize style of confidence intervals
                   #fun="event",#ÀÛ»ý·çÏÕ
                   #xscale="d_y",#ÄêÔÂÈÕ×ª»»
                   ggtheme = theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))+theme(legend.background = element_rect(fill="NA", size=0.5, linetype="solid",colour ="NA"))+theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),plot.title = element_text(hjust = 0.5, face = "bold")) # Change ggplot2 theme
                   
)

cph = coxph(formula =as.formula(paste("Surv(time = as.numeric(data[,time]),event= as.numeric(data[,event])) ~ factor(`",variable,"`)",sep="")),data=data)
tmp=summary(cph)

pval=format(tmp$coefficients[,5],scientific=TRUE,digit=3)
HR=tmp$conf.int[,c(1,3,4)]
HR=signif(HR, digits = 2)

if(is.null(nrow(HR))){
  tabHR = as.data.frame(t(HR))
  HR2 = as.data.frame(t(HR))
}else{
  tabHR=as.data.frame(HR)
  HR2=as.data.frame(HR)
}
  
  
tabHR['P']=data.frame(pval)
tabHR=data.frame(tabHR)
colnames(tabHR)=c("HR","LOW","HIGH","pval")
write.csv(tabHR,file = paste("out/restb",variable,tm,".csv",sep=""),quote = F)
# tabHR=paste("HR",":",HR[1],"(",HR[2],"-",HR[3],")",sep="")
# restb=c(variable,pval,tabHR)
# print(restb)
# write.csv(restb,file = paste("out/restb",variable,tm,".csv",sep=""),sep = "\t",quote = F,row.names = F,col.names =F)

for (i in 1:nrow(tabHR)) {
  HR13=paste("HR",":",HR2[i,1],"(",HR2[i,2],"-",HR2[i,3],") p=",pval[i],sep="")
  st=paste("atop(bold('",HR13,"'))",sep="")
  ct=0.05
  xtm=round(max(data[time])/5)*5
  ggsurv$plot=ggsurv$plot+annotate("segment", x =  0.66*xtm, xend =  0.68*xtm, y =(0.85-i*ct), yend =(0.85-i*ct),colour = palettecol[i+1],size=2)+
    annotate("text",x =0.70*xtm , y =(0.82-i*ct),parse = T,label = "atop(bold('Vs.'))",size = 4.2)+annotate("text",x =0.9*xtm, y =(0.82-i*ct),parse = T,label = st,size = 4.2)+
    annotate("segment", x = 0.72*xtm, xend = 0.74*xtm, y =(0.85-i*ct), yend =(0.85-i*ct),colour = palettecol[1],size=2)
}



fname=paste("out/",title,"_",variable,"-",tm,sep = "")

png(file=paste(fname,".png",sep = ""),width = 938*3,height = 752*3,res=300)
print(ggsurv)
dev.off()

pdf(file=paste(fname,".pdf",sep = ""),width = 10,height = 8,onefile = FALSE)
print(ggsurv)
dev.off()

tiff(file=paste(fname,".tiff",sep = ""),width = 938*3,height = 752*3,res = 300)
print(ggsurv)
dev.off()













