rm(list = ls(all = TRUE))

suppressMessages(library(rms))
suppressMessages(library(survminer))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))

setwd("data")

args = commandArgs()
ca = args[(grep("--ca",args)+1)]
timic = args[(grep("--timic",args)+1)]
conf.int = args[(grep("--confint",args)+1)]
risk.table = args[(grep("--risktable",args)+1)]
HRR = args[(grep("--HRR",args)+1)]
low = args[(grep("--low",args)+1)]
high = args[(grep("--high",args)+1)]
surv.median.line = args[(grep("--survmedianline",args)+1)]
cutoff = args[(grep("--cutoff",args)+1)]
customlow = args[(grep("--customlow",args)+1)]
customhigh = args[(grep("--customhigh",args)+1)]
datatype = args[(grep("--datatype",args)+1)]
survivaltype = args[(grep("--survivaltype",args)+1)]
conditions = args[(grep("--conditions",args)+1)]
conditiongene = args[(grep("--conditiongene",args)+1)]
cli = args[(grep("--cli",args)+1)]
panclinicalcut = args[(grep("--panclinicalcut",args)+1)]
tm = args[(grep("--tm",args)+1)]

genelist = args[(grep("--genelist",args)+1)]
direction = args[(grep("--direction",args)+1)]

print(Sys.time())

palette=c(low,high)
genes=unlist(strsplit(genelist, split="+",fixed=T))



conf.int=eval(parse(text = conf.int))
risk.table=eval(parse(text = risk.table))
HRR=eval(parse(text = HRR))

title=ca
time=paste(survivaltype,'.time',sep='')
event=survivaltype
select=c('sample','group',genes,event,time)
print(select)

kmplot=function(expcli,time,event,datatype,variable,cutoff,low,high,title,timic,palette,conf.int, risk.table,surv.median.line,HRR){

  data=as.data.frame(expcli)
  deltime=which(is.na(data[time]))
  if (length(deltime)>0) {
    data=data[-deltime,]
  }
  delevent=which(is.na(data[event]))
  if (length(delevent)>0) {
    data=data[-delevent,]
  }
  data[is.na(data)]<-0

  if(datatype!="Mutation"&datatype!="CNA"&datatype!="fusion"){

    if(cutoff=="Median"){
      legend.labs = c(paste("low",variable,sep = " "), paste("high",variable,sep = " "))
      data[,variable]=ifelse(data[,variable]<=summary(data[,variable])[3],"low","high")
      data[,variable]=factor(data[,variable],levels =c("low","high"))
      palette=c(low,high)

    }else if(cutoff=="Quartile3"){
      legend.labs = c(paste("low",variable,sep = " "),paste("median",variable,sep = " "), paste("high",variable,sep = " "))
      data[,variable]=ifelse(data[,variable]<=summary(data[,variable])[2],"low",ifelse(data[,variable]>=summary(data[,variable])[5],"high","median"))
      data[,variable]=factor(data[,variable],levels =c("low","median","high"))
      palette=c(low,"green",high)
    }else if(cutoff=="Quartile"){
      legend.labs = c(paste("low",variable,"<1st",sep = " "), paste("high",variable,">3rd",sep = " "))
      data[,variable]=ifelse(data[,variable]<=summary(data[,variable])[2],"low",ifelse(data[,variable]>=summary(data[,variable])[5],"high",NA))
      data=data[which(!is.na(data[,variable])),]
      data[,variable]=factor(data[,variable],levels =c("low","high"))
      palette=c(low,high)
    }else if(cutoff=="bestcut"){
      res.cut=surv_cutpoint(data, time = time, event = event, variable,minprop = 0.1, progressbar = TRUE)
      res.point=res.cut$cutpoint$cutpoint
      legend.labs = c(paste("low",variable,sep = " "), paste("high",variable,sep = " "))
      data[,variable]=ifelse(data[,variable]<=res.point,"low","high")
      data[,variable]=factor(data[,variable],levels =c("low","high"))
      palette=c(low,high)
    }else{
      legend.labs = c(paste("low",variable,sep = " "), paste("high",variable,sep = " "))
      data[,variable]=ifelse(data[,variable]<=summary(data[,variable])[3],"low","high")
      data[,variable]=factor(data[,variable],levels =c("low","high"))
      palette=c(low,high)
    }

  }else if(datatype=="CNA"){
    data=data[-which(data[,variable]==0),]
    data[,variable]=ifelse(data[,variable]<0,"loss","gain")
    data[,variable]=factor(data[,variable],levels =c("loss","gain"))
    legend.labs = c(paste(variable,"loss",sep = " "), paste(variable,"gain",sep = " "))
    palette=c(low,high)
    print(legend.labs)
  }else if(datatype=="fusion"){
    legend.labs = c(paste(variable,"with fusion",sep = " "), paste(variable,"without fusion",sep = " "))
    palette=c(low,high)
    print(legend.labs)
  }else{
    legend.labs = c(paste(variable,"wild type",sep = " "), paste(variable,"mutation",sep = " "))
    palette=c(low,high)
  }


if(timic=="years")
{data[,time]=data[,time]/365}else if(timic=="months")
{data[,time]=data[,time]/30}

  ddist = datadist(data)
  options(datadist = "ddist")
  S.OS = Surv(time = as.numeric(data[,time]),event= as.numeric(data[,event]))
  f.cox=coxph(as.formula(paste("S.OS~`",variable,"`",sep="")),data=data)
  tmp=summary(f.cox)
  pval=signif(tmp$sctest[3], digits = 3)
  HR=tmp$conf.int[,c(1,3,4)]
  HR=signif(HR, digits = 3)
  tabHR=paste("HR",":",HR[1],"(",HR[2],"-",HR[3],")",sep="")
  #f.np = npsurv(formula =as.formula(paste("S.OS ~ factor(",variable,")",sep="")),data=data)
  f.np = npsurv(formula =as.formula(paste("Surv(time = as.numeric(data[,time]),event= as.numeric(data[,event])) ~ factor(`",variable,"`)",sep="")),data=data)
  restb=c(variable,pval,tabHR)
  write.table(restb,file = paste("out/restb",variable,tm,".txt",sep=""),sep = "\t",quote = F,row.names = F,col.names =F)
  forestplot=c(title,pval,HR)


  ggsurv=ggsurvplot( fit=f.np,
                     data = data,
                     main = "Survival curve",
                     title=title,
                     xlim=c(0,max(data[,time])),
                     #break.time.by=0.2*max(data$stime),#XÖácut
                     size = 1.1, # change line size
                     xlab = paste("Time in ",timic,sep = ""),# customize X axis label.
                     palette = palette,# custom color palettes
                     conf.int = conf.int, # Add confidence interval
                     pval = TRUE, # Add p-value
                     risk.table = risk.table, # Add risk table
                     risk.table.col = "strata",# Risk table color by groups
                     surv.median.line = surv.median.line, # add the median survival pointer.
                     legend.title = "",
                     legend.labs = legend.labs, # Change legend labels
                     risk.table.y.text = FALSE, # hide risktabl y
                     risk.table.height = 0.25, # Useful to change when you have multiple groups
                     risk.table.fontsize=5,
                     font.tickslab = c(12, "bold"),#km×ø±êÖá
                     font.legend=c(12, "bold"),#legend´óÐ¡ÑÕÉ«
                     #font.x = c(5, "bold", "red"),xÖá±êÇ©´óÐ¡
                     conf.int.style = "step", # customize style of confidence intervals
                     #fun="event",#ÀÛ»ý·çÏÕ
                     ggtheme = theme_bw()+theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),plot.title = element_text(hjust = 0.5, face = "bold")) # Change ggplot2 theme

  )
  if(HRR){ggsurv$plot=ggsurv$plot+annotate("text",x =max(data[,time])/6.5 , y = 0.1,parse = T,label = tabHR,size = 5)}

  fname=paste("out/",title,"_",variable,"-",tm,sep = "")

  png(file=paste(fname,".png",sep = ""),width = 816*3,height = 641*3,res=300)
  print(ggsurv)
  dev.off()

  pdf(file=paste(fname,".pdf",sep = ""),width = 10,height = 8,onefile = FALSE)
  print(ggsurv)
  dev.off()

  tiff(file=paste(fname,".tiff",sep = ""),width = 816*3,height = 641*3,res = 300)
  print(ggsurv)
  dev.off()

  return(forestplot)
}


kmpval=function(expcli,time,event,datatype,variable,cutoff,low,high,title,timic,palette,conf.int, risk.table,surv.median.line,HRR){


  data=as.data.frame(expcli)
  deltime=which(is.na(data[time]))
  if (length(deltime)>0) {
    data=data[-deltime,]
  }
  delevent=which(is.na(data[event]))
  if (length(delevent)>0) {
    data=data[-delevent,]
  }
  data[is.na(data)]<-0

  if(datatype!="Mutation"&datatype!="CNA"&datatype!="fusion"){

    if(cutoff=="Median"){
      legend.labs = c(paste("low",variable,sep = " "), paste("high",variable,sep = " "))
      data[,variable]=ifelse(data[,variable]<=summary(data[,variable])[3],"low","high")
      data[,variable]=factor(data[,variable],levels =c("low","high"))
      palette=c(low,high)
    }else if(cutoff=="Quartile3"){
      legend.labs = c(paste("low",variable,sep = " "),paste("median",variable,sep = " "), paste("high",variable,sep = " "))
      data[,variable]=ifelse(data[,variable]<=summary(data[,variable])[2],"low",ifelse(data[,variable]>=summary(data[,variable])[5],"high","median"))
      data[,variable]=factor(data[,variable],levels =c("low","median","high"))
      palette=c(low,"green",high)
    }else if(cutoff=="Quartile"){
      legend.labs = c(paste("low",variable,"<1st",sep = " "), paste("high",variable,">3rd",sep = " "))
      data[,variable]=ifelse(data[,variable]<=summary(data[,variable])[2],"low",ifelse(data[,variable]>=summary(data[,variable])[5],"high",NA))
      data=data[which(!is.na(data[,variable])),]
      data[,variable]=factor(data[,variable],levels =c("low","high"))
      palette=c(low,high)
    }else if(cutoff=="bestcut"){
      res.cut=surv_cutpoint(data, time = time, event = event, variable,minprop = 0.1, progressbar = TRUE)
      res.point=res.cut$cutpoint$cutpoint
      legend.labs = c(paste("low",variable,sep = " "), paste("high",variable,sep = " "))
      data[,variable]=ifelse(data[,variable]<=res.point,"low","high")
      data[,variable]=factor(data[,variable],levels =c("low","high"))
      palette=c(low,high)
    }

  }else if(datatype=="CNA"){
    data=data[-which(data[,variable]==0),]
    data[,variable]=ifelse(data[,variable]<0,"loss","gain")
    data[,variable]=factor(data[,variable],levels =c("loss","gain"))
    legend.labs = c(paste(variable,"loss",sep = " "), paste(variable,"gain",sep = " "))
    palette=c(low,high)
    print(legend.labs)
  }else if(datatype=="fusion"){
    legend.labs = c(paste(variable,"with fusion",sep = " "), paste(variable,"without fusion",sep = " "))
    palette=c(low,high)
    print(legend.labs)
  }else{
    legend.labs = c(paste(variable,"wild type",sep = " "), paste(variable,"mutation",sep = " "))
    palette=c(low,high)
  }


if(timic=="years")
{data[,time]=data[,time]/365}else if(timic=="months")
{data[,time]=data[,time]/30}

  ddist = datadist(data)
  options(datadist = "ddist")
  S.OS = Surv(time = as.numeric(data[,time]),event= as.numeric(data[,event]))
  f.cox=coxph(as.formula(paste("S.OS~`",variable,"`",sep="")),data=data)
  tmp=summary(f.cox)
  pval=signif(tmp$sctest[3], digits = 3)
  HR=tmp$conf.int[,c(1,3,4)]
  HR=signif(HR, digits = 3)
  tabHR=paste("HR",":",HR[1],"(",HR[2],"-",HR[3],")",sep="")
  #f.np = npsurv(formula =as.formula(paste("S.OS ~ factor(",variable,")",sep="")),data=data)
  f.np = npsurv(formula =as.formula(paste("Surv(time = as.numeric(data[,time]),event= as.numeric(data[,event])) ~ factor(`",variable,"`)",sep="")),data=data)
  restb=c(variable,pval,tabHR)

  dangtb=c(variable,pval,HR)

  return(dangtb)
}



expcli=NULL
dangtb=NULL

tmp=fread(paste('rdata/',datatype,'/',ca,'.txt',sep=''),stringsAsFactors = F,select = select)
tmp=as.data.frame(tmp)
tmp1=tmp

matchgene=intersect(genes,colnames(tmp))
nomatchgene=setdiff(genes,colnames(tmp))

write.table(matchgene,file = paste("out/matchgene",tm,".txt",sep=""),sep = "\t",quote = F,row.names = F,col.names =F)
write.table(nomatchgene,file = paste("out/nomatchgene",tm,".txt",sep=""),sep = "\t",quote = F,row.names = F,col.names =F)

genes=matchgene
genelist=paste(genes, collapse = "+")
print(genelist)

for (i in 1:length(genes)) {
  variable=genes[i]
  print(variable)

  if(conditiongene!=''){

    datatype2=strsplit(conditions, "_")[[1]][1]
    cdit=strsplit(conditions, "_")[[1]][2]
    cditgenes=str_split(conditiongene,';')[[1]]

    if(cli!='no'){
      print(ca)
      cditselect=c('sample',cditgenes,cli)
      tmpmut=fread(paste('rdata/',datatype2,'/',ca,'.txt',sep=''),stringsAsFactors = F,select = cditselect)

      conditioncli=paste("which(tmpmut$",cli,"=='",panclinicalcut,"')",sep='')
      clipick=eval(parse(text = conditioncli))
      tmpmut=tmpmut[clipick,]


    }else{
      cditselect=c('sample',cditgenes)
      tmpmut=fread(paste('rdata/',datatype2,'/',ca,'.txt',sep=''),stringsAsFactors = F,select = cditselect)
    }

    tmpmut=as.data.frame(tmpmut)
  if(datatype2=="Methylation"){


    if(cdit=='high'){
       conditiongenestr=str_replace_all(conditiongene,";",paste(">=",0.8,"|tmpmut$",sep=""))
       conditiongenestr=paste("which(tmpmut$",conditiongenestr,">=",0.8,")",sep='')

      }else{
        conditiongenestr=str_replace_all(conditiongene,";",paste("<=",0.2,"|tmpmut$",sep=""))
        conditiongenestr=paste("which(tmpmut$",conditiongenestr,"<=",0.2,")",sep='')
      }

  }else if(datatype2=="CNA"){


      if(cdit=='gain'){
        conditiongenestr=str_replace_all(conditiongene,";",paste(">",0,"|tmpmut$",sep=""))
       conditiongenestr=paste("which(tmpmut$",conditiongenestr,">",0,")",sep='')

      }else{
         conditiongenestr=str_replace_all(conditiongene,";",paste("<",0,"|tmpmut$",sep=""))
        conditiongenestr=paste("which(tmpmut$",conditiongenestr,"<",0,")",sep='')
      }

  }else{
     conditiongenestr=str_replace_all(conditiongene,";",paste("==",cdit,"|tmpmut$",sep=""))
    conditiongenestr=paste("which(tmpmut$",conditiongenestr,"==",cdit,")",sep='')
    }
    pickmut=eval(parse(text = conditiongenestr))
    pickmutsample=tmpmut[pickmut,'sample']
    print(199)

    tmp=tmp[which(tmp[,'sample'] %in% pickmutsample),]
    tmp=tmp[which(tmp[,'group']==ca),]
    km=kmpval(expcli=tmp,time,event,datatype,variable,cutoff,low,high,title=ca,timic,palette,conf.int, risk.table,surv.median.line,HRR)
    dangtb=rbind(dangtb,km)

  }else{
    if(cli!='no'){

      cditselect=c('sample',cli)
      tmpmut=fread('rdata/pancliall.txt',stringsAsFactors = F,select = cditselect)
      tmpmut=as.data.frame(tmpmut)
      conditiongenestr=paste("which(tmpmut$",cli,"=='",panclinicalcut,"')",sep='')
      pickmut=eval(parse(text = conditiongenestr))
      pickmutsample=tmpmut[pickmut,'sample']
      tmp=tmp[which(tmp[,'sample'] %in% pickmutsample),]
      tmp=tmp[which(tmp[,'group']==ca),]
      km=kmpval(expcli=tmp,time,event,datatype,variable,cutoff,low,high,title=ca,timic,palette,conf.int, risk.table,surv.median.line,HRR)
      dangtb=rbind(dangtb,km)

    }else{
      tmp=tmp[which(tmp[,'group']==ca),]
      km=kmpval(expcli=tmp,time,event,datatype,variable,cutoff,low,high,title=ca,timic,palette,conf.int, risk.table,surv.median.line,HRR)
      dangtb=rbind(dangtb,km)
    }
  }

}

expcli=tmp



  deltime=which(is.na(expcli[time]))
  if (length(deltime)>0) {
    expcli=expcli[-deltime,]
  }
  delevent=which(is.na(expcli[event]))
  if (length(delevent)>0) {
    expcli=expcli[-delevent,]
  }


print(Sys.time())

rownames(dangtb)=NULL

dangtb1=as.data.frame(dangtb)

colnames(dangtb)=c("symbol","pvalue","HR","low.95.CI","high.05.CI")
dangtb1=dangtb[,c(1,3,4,5,2)]
rownames(dangtb1)=NULL
options(stringsAsFactors = F)
dangtb1=as.data.frame(dangtb1)
write.table(dangtb1,file = paste("out/dangtb",tm,".txt",sep=""),sep = "\t",quote = F,row.names =F)


data=as.data.frame(expcli)
head(data)
ddist=datadist(data)
options(datadist="ddist")
S.OS = Surv(time = as.numeric(data[,time]),event= as.numeric(data[,event]))
f.m1=cph(as.formula(paste("S.OS~",genelist )),data=data,x=T,y=T,surv = TRUE)
f.m2=coxph(as.formula(paste("S.OS~",genelist )),data=data)



#can be one of "both", "backward", or "forward",both逐步回归
if(direction!='none'){
  f.step2=step(f.m2,direction = direction)
}else{
  f.step2=f.m2
}
aa=summary(f.step2)
mut=data.frame(aa$conf.int[,c(1,3,4)],aa$coefficients[,c(5,1)])
colnames(mut)=c("HR","low.95.CI","high.95.CI","pvalue","coef")
mut=signif(mut, digits = 3)
mut["symbol"]=rownames(mut)
mut=mut[,c(6,1,2,3,4,5)]
write.table(mut,file = paste("out/duotb",tm,".txt",sep=""),sep = "\t",quote = F,row.names = F)
mut["HR.%95CI"]=paste(mut[,2],":(",mut[,3],"-",mut[,4],")",sep="")
mut=mut[,c("symbol","HR.%95CI","pvalue","coef")]

dangtb1["HR.%95CI"]=paste(dangtb1[,2],":(",dangtb1[,3],"-",dangtb1[,4],")",sep="")
dangtb1=dangtb1[,c("symbol","HR.%95CI","pvalue")]

all=merge(dangtb1,mut,by.x = "symbol",by.y = "symbol")
colnames(all)=c("symbol","uni.HR.%95CI","uni.pvalue","mut.HR.%95CI","mut.pvalue","coef")
write.table(all,file = paste("out/alltb",tm,".txt",sep=""),sep = "\t",quote = F,row.names =F)



if(direction!='none'){
  f.step=step(f.m1,direction = direction)
}else{
  f.step=f.m1
}


f=f.step
nos <- function(x) {
  pre=predict(f,newdata=x)
  return(pre)
}

kk=as.data.frame(f$coefficients)
pick=data[,rownames(kk)]
risk=apply(pick,1,nos)
print(Sys.time())
risk=risk+f$center

# ###################
newdata=cbind(risk,data)

data=newdata
variable="risk"
print(66666)
print(summary(data[,variable]))
if(datatype=='Mutation'|datatype=='CNA'|datatype=='fusion'){
	datatype='risk level'
	cutoff='risk cut'
}

km=kmplot(expcli=data,time,event,datatype,variable,cutoff,low,high,title=ca,timic,palette,conf.int, risk.table,surv.median.line,HRR)
print(Sys.time())