rm(list = ls(all = TRUE))

suppressMessages(library(rms))
suppressMessages(library(survminer))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(forestplot))

setwd("data")

args = commandArgs()
ca = args[(grep("--ca",args)+1)]
variable = args[(grep("--variable",args)+1)]
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
Methylationhigh = args[(grep("--Methylationhigh",args)+1)]
Methylationlow = args[(grep("--Methylationlow",args)+1)]
cpg = args[(grep("--cpg",args)+1)]
cli = args[(grep("--cli",args)+1)]
panclinicalcut = args[(grep("--panclinicalcut",args)+1)]
pathway = args[(grep("--pathway",args)+1)]
tm = args[(grep("--tm",args)+1)]

print(cli)

#conditiongene=eval(parse(text = eval(parse(text = conditiongene))))

print(ca)

conf.int=eval(parse(text = conf.int))
risk.table=eval(parse(text = risk.table))
HRR=eval(parse(text = HRR))

title=ca
cas=str_split(ca,';')[[1]]
time=paste(survivaltype,'.time',sep='')
event=survivaltype
select=c('sample','group',variable,event,time)

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


  expcli=NULL
  forestplot=NULL
  for (i in 1:length(cas)) {
    ca=cas[i]

if(conditiongene!=''){
  print(11111)

  datatype2=strsplit(conditions, "_")[[1]][1]
  cdit=strsplit(conditions, "_")[[1]][2]
  cditgenes=str_split(conditiongene,';')[[1]]

  if(cli!='no'){
    print(6666)
    print(ca)
    cditselect=c('sample',cditgenes,cli)
    print(cditselect)
    tmpmut=fread(paste('rdata/pancancer/',datatype2,'/',ca,'.txt',sep=''),stringsAsFactors = F,select = cditselect)

    conditioncli=paste("which(tmpmut$",cli,"=='",panclinicalcut,"')",sep='')
    clipick=eval(parse(text = conditioncli))
    tmpmut=tmpmut[clipick,]
    print(tmpmut[1:3,1:2])


  }else{
     print(2222)
    cditselect=c('sample',cditgenes)
    tmpmut=fread(paste('rdata/pancancer/',datatype2,'/',ca,'.txt',sep=''),stringsAsFactors = F,select = cditselect)
  }

  tmpmut=as.data.frame(tmpmut)
  conditiongenestr=str_replace_all(conditiongene,";",paste("==",cdit,"|tmpmut$",sep=""))
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
  print(pickmutsample)
  tmp=fread(paste('rdata/pancancer/',datatype,'/',ca,'.txt',sep=''),stringsAsFactors = F,select = select)
  tmp=as.data.frame(tmp)

  tmp=tmp[which(tmp[,'sample'] %in% pickmutsample),]
    print("_____________________________________")
  print(tmp[1:3,1:3])
   print("_____________________________________")
  km=kmplot(expcli=tmp,time,event,datatype,variable,cutoff,low,high,title=ca,timic,palette,conf.int, risk.table,surv.median.line,HRR)
  forestplot=rbind(forestplot,km)
  expcli=rbind(expcli,tmp)

}else{
  if(cli!='no'){
    print(333333)
    cditselect=c('sample',cli)
    tmpmut=fread('rdata/pancancer/pancliall.txt',stringsAsFactors = F,select = cditselect)
    tmpmut=as.data.frame(tmpmut)
    conditiongenestr=paste("which(tmpmut$",cli,"=='",panclinicalcut,"')",sep='')
    pickmut=eval(parse(text = conditiongenestr))
    pickmutsample=tmpmut[pickmut,'sample']

    #print(pickmutsample)

    tmp=fread(paste('rdata/pancancer/',datatype,'/',ca,'.txt',sep=''),stringsAsFactors = F,select = select)
    tmp=as.data.frame(tmp)
    print(tmp[1:10,1:4])
    tmp=tmp[which(tmp[,'sample'] %in% pickmutsample),]




    km=kmplot(expcli=tmp,time,event,datatype,variable,cutoff,low,high,title=ca,timic,palette,conf.int, risk.table,surv.median.line,HRR)
    forestplot=rbind(forestplot,km)
    expcli=rbind(expcli,tmp)

  }else{
    print(44444)
    tmp=fread(paste('rdata/pancancer/',datatype,'/',ca,'.txt',sep=''),stringsAsFactors = F,select = select)

    km=kmplot(expcli=tmp,time,event,datatype,variable,cutoff,low,high,title=ca,timic,palette,conf.int, risk.table,surv.median.line,HRR)
    forestplot=rbind(forestplot,km)
    expcli=rbind(expcli,tmp)
  }
}

  }
  km=kmplot(expcli,time,event,datatype,variable,cutoff,low,high,title='pan-cancer',timic,palette,conf.int, risk.table,surv.median.line,HRR)
  forestplot=rbind(forestplot,km)
  rownames(forestplot)=1:nrow(forestplot)


#数据处理

sample<-data.frame(forestplot,stringsAsFactors=FALSE)
sample[nrow(sample),1]='pan-cancer'

sample=sample[order(sample[,2]),]


tabletext1<-as.character(sample[,1])

tabletext2<-as.numeric(sample[,2])

tabletext3<-paste(round(as.numeric(sample[,3]),3),round(as.numeric(sample[,4]),2),sep="(")

tabletext4<-paste(tabletext3,round(as.numeric(sample[,5]),2),sep="-")

tabletext5<-paste0(tabletext4,sep=")")

tabletext<-cbind(tabletext1,tabletext2,tabletext5)

#tabletext=as.character(tabletext)

tabletext=rbind(c('ca','logrank-p','HR(95%CI)'),tabletext)


csize <- data.frame(mean=c(NA,round(as.numeric(sample[,3]),3)),
                    lower=c(NA,round(as.numeric(sample[,4]),2)),
                    upper=c(NA, round(as.numeric(sample[,5]),2)))






print(csize)

fname=paste("out/forestplot_",variable,"-",tm,sep = "")

png(file=paste(fname,".png",sep = ""),width = 816*3,height = 641*3,res=300)
forestplot(labeltext=tabletext, #文本信息

           csize,

           hrzl_lines=gpar(col="black",lty=1,lwd =2),

           is.summary=c(TRUE,rep(FALSE,5)),

           boxsize = 0.15,##大小

           graph.pos=3,#图在表中的列位置

           graphwidth = unit(0.3,"npc"),#图在表中的宽度比例

           fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石,可以更改fpDrawNormalCI；fpDrawCircleCI等

           col=fpColors(box="steelblue", lines="black", zero = "black"),#颜色设置

           lwd.ci=2,ci.vertices.height = 0.1,ci.vertices=TRUE,#置信区间用线宽、高、型

           zero=1,#zero线横坐标

           lwd.zero=2,#zero线宽

           grid=T,

           lwd.xaxis=2,#X轴线宽

           txt_gp = fpTxtGp(ticks = gpar(cex = 1)), # configure fontsize
           #txt_gp = fpTxtGp(ticks = gpar(cex = 1), xlab = gpar(cex = 1.5), cex = 1.2), # configure fontsize
           title="Hazard Ratio",

           xlab="",#X轴标题

           clip=c(-Inf,3),#边界

           colgap = unit(0.5,"cm"),#列间隙

           new_page = TRUE#是否新页

)

dev.off()


pdf(file=paste(fname,".pdf",sep = ""),width = 10,height = 8,onefile = FALSE)
aa=forestplot(labeltext=tabletext, #文本信息

           csize,

           hrzl_lines=gpar(col="black",lty=1,lwd =2),

           is.summary=c(TRUE,rep(FALSE,5)),

           boxsize = 0.15,##大小

           graph.pos=3,#图在表中的列位置

           graphwidth = unit(0.3,"npc"),#图在表中的宽度比例

           fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石,可以更改fpDrawNormalCI；fpDrawCircleCI等

           col=fpColors(box="steelblue", lines="black", zero = "black"),#颜色设置

           lwd.ci=2,ci.vertices.height = 0.1,ci.vertices=TRUE,#置信区间用线宽、高、型

           zero=1,#zero线横坐标

           lwd.zero=2,#zero线宽

           grid=T,

           lwd.xaxis=2,#X轴线宽

           txt_gp = fpTxtGp(ticks = gpar(cex = 1)), # configure fontsize
           #txt_gp = fpTxtGp(ticks = gpar(cex = 1), xlab = gpar(cex = 1.5), cex = 1.2), # configure fontsize
           title="Hazard Ratio",

           xlab="",#X轴标题

           clip=c(-Inf,3),#边界

           colgap = unit(0.5,"cm"),#列间隙

           new_page = TRUE#是否新页

)

dev.off()


tiff(file=paste(fname,".tiff",sep = ""),width = 816*3,height = 641*3,res = 300)
forestplot(labeltext=tabletext, #文本信息

           csize,

           hrzl_lines=gpar(col="black",lty=1,lwd =2),

           is.summary=c(TRUE,rep(FALSE,5)),

           boxsize = 0.15,##大小

           graph.pos=3,#图在表中的列位置

           graphwidth = unit(0.3,"npc"),#图在表中的宽度比例

           fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石,可以更改fpDrawNormalCI；fpDrawCircleCI等

           col=fpColors(box="steelblue", lines="black", zero = "black"),#颜色设置

           lwd.ci=2,ci.vertices.height = 0.1,ci.vertices=TRUE,#置信区间用线宽、高、型

           zero=1,#zero线横坐标

           lwd.zero=2,#zero线宽

           grid=T,

           lwd.xaxis=2,#X轴线宽

           txt_gp = fpTxtGp(ticks = gpar(cex = 1)), # configure fontsize
           #txt_gp = fpTxtGp(ticks = gpar(cex = 1), xlab = gpar(cex = 1.5), cex = 1.2), # configure fontsize
           title="Hazard Ratio",

           xlab="",#X轴标题

           clip=c(-Inf,3),#边界

           colgap = unit(0.5,"cm"),#列间隙

           new_page = TRUE#是否新页

)

dev.off()

res=as.data.frame(tabletext)

write.csv(res[2:nrow(res),],file = paste("out/restb",variable,tm,".csv",sep=""),quote = F,row.names = F)