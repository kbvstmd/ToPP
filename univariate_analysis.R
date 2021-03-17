start=Sys.time()
print(start)
rm(list = ls(all = TRUE))

suppressMessages(library(rms))
suppressMessages(library(survminer))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(forestplot))
suppressMessages(library(ggpubr))
suppressMessages(library(ggsignif))

setwd("data")

print(paste("setp 1 takes ",Sys.time(),sep=""))

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
cli = args[(grep("--cli",args)+1)]
panclinicalcut = args[(grep("--panclinicalcut",args)+1)]
pathway = args[(grep("--pathway",args)+1)]
tm = args[(grep("--tm",args)+1)]




#conditiongene=eval(parse(text = eval(parse(text = conditiongene))))

conf.int=eval(parse(text = conf.int))
risk.table=eval(parse(text = risk.table))
HRR=eval(parse(text = HRR))

title=ca
cas=str_split(ca,';')[[1]]
time=paste(survivaltype,'.time',sep='')
event=survivaltype
select=c('sample','group',variable,event,time)


print(paste("setp 2 takes ",Sys.time(),sep=""))

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

if(conditiongene!=''){

  print(111111)

  datatype2=strsplit(conditions, "_")[[1]][1]
  print(datatype2)
  cdit=strsplit(conditions, "_")[[1]][2]
  cditgenes=str_split(conditiongene,';')[[1]]

  if(cli!='no'){
    cditselect=c('sample','group',cditgenes,cli)
    tmpmut=fread(paste('rdata/',datatype2,'/',ca,'.txt',sep=''),stringsAsFactors = F,select = cditselect)

    conditioncli=paste("which(tmpmut$",cli,"=='",panclinicalcut,"')",sep='')
    clipick=eval(parse(text = conditioncli))
    tmpmut=tmpmut[clipick,]


  }else{
    cditselect=c('sample','group',cditgenes)
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
  print(conditiongenestr)
  pickmut=eval(parse(text = conditiongenestr))
  pickmutsample=tmpmut[pickmut,'sample']
  print(199)
  tmp=fread(paste('rdata/',datatype,'/',ca,'.txt',sep=''),stringsAsFactors = F,select = select)
  tmp=as.data.frame(tmp)
  tmp=tmp[which(tmp[,'sample'] %in% pickmutsample),]
  tmp1=tmp
  tmp=tmp[which(tmp[,'group']==ca),]

  km=kmplot(expcli=tmp,time,event,datatype,variable,cutoff,low,high,title=ca,timic,palette,conf.int, risk.table,surv.median.line,HRR)

}else{
  if(cli!='no'){
     print(222222)

    cditselect=c('sample','group',cli)
    tmpmut=fread('rdata/pancliall.txt',stringsAsFactors = F,select = cditselect)
    tmpmut=as.data.frame(tmpmut)
    conditiongenestr=paste("which(tmpmut$",cli,"=='",panclinicalcut,"')",sep='')
    pickmut=eval(parse(text = conditiongenestr))
    pickmutsample=tmpmut[pickmut,'sample']

    tmp=fread(paste('rdata/',datatype,'/',ca,'.txt',sep=''),stringsAsFactors = F,select = select)
    tmp=as.data.frame(tmp)
    tmp=tmp[which(tmp[,'sample'] %in% pickmutsample),]
    tmp1=tmp
    tmp=tmp[which(tmp[,'group']==ca),]
    km=kmplot(expcli=tmp,time,event,datatype,variable,cutoff,low,high,title=ca,timic,palette,conf.int, risk.table,surv.median.line,HRR)

  }else{

     print(33333)

     print(paste("setp 3 takes ",Sys.time(),sep=""))

    tmp=fread(paste('rdata/',datatype,'/',ca,'.txt',sep=''),stringsAsFactors = F,select = select)
    print(select)
    tmp1=tmp
    tmp=tmp[which(tmp[,'group']==ca),]
    print(tmp)
    km=kmplot(expcli=tmp,time,event,datatype,variable,cutoff,low,high,title=ca,timic,palette,conf.int, risk.table,surv.median.line,HRR)
    print(paste("setp 3 takes ",Sys.time(),sep=""))
  }
}




# boxplot ca vs normal

zh.sig.boxplot <- function(mdat=data.frame(),x="",y="",fill="",add="jitter",
                           title="",subtitle="",
                           xlab="",ylab="",sig =T,
                           test.method="wilcox.test",caption=""){#comparisons=list(),
  if(!(all(c(x,y) %in% colnames(mdat)))){"please check your input!!!";break}
  if(fill ==""){fill=x}
  if(xlab ==""){xlab=x}
  if(ylab ==""){ylab=y}
  if(title==""){title=paste(x,y,sep = " ")}
  mgrps <- levels(mdat[,fill])

  if(sig){
    sumx <- summary(mdat[,y])
    mrge <- sumx[6] - sumx[1]

    if(length(levels(mdat[,fill])) == 2){comparisons=list(c(mgrps[1],mgrps[2]));
    y_pos=c(mrge*1.02 + sumx[1])}

    if(length(levels(mdat[,fill])) == 3){comparisons=list(c(mgrps[1],mgrps[2]),c(mgrps[2],mgrps[3]),c(mgrps[1],mgrps[3]));
    y_pos=c(mrge*1.02 + sumx[1],mrge*1.08 + sumx[1],mrge*1.14 + sumx[1])}
    if(length(levels(mdat[,fill])) == 4){comparisons=list(c(mgrps[1],mgrps[2]),c(mgrps[2],mgrps[3]),c(mgrps[3],mgrps[4]),
                                                          c(mgrps[1],mgrps[3]),c(mgrps[2],mgrps[4]),c(mgrps[1],mgrps[4]));
    y_pos=c(mrge*1.02 + sumx[1],mrge*1.08 + sumx[1],mrge*1.14 + sumx[1],
            mrge*1.20 + sumx[1],mrge*1.26 + sumx[1],mrge*1.32 + sumx[1])}

    ggboxplot(mdat,
              x = x,
              y = y,
              fill = fill,
              bxp.errorbar = T,
              bxp.errorbar.width = 0.2,
              palette = "npg",
              add = add) +
      labs(title = title,
           subtitle = subtitle,
           caption =caption,
           x = xlab,
           y = ylab) +
      geom_signif(comparisons = comparisons,
                  y_position = y_pos,
                  tip_length = c(0),
                  map_signif_level = T,
                  test = test.method) +
      theme(plot.title    = element_text(color = "black", size   = 16, hjust = 0.5),
            plot.subtitle = element_text(color = "black", size   = 14,hjust = 0.5),
            plot.caption  = element_text(color = "black", size   = 16,face = "italic", hjust = 1),
            axis.text.x   = element_text(color = "black", size = 16, angle = 0),
            axis.text.y   = element_text(color = "black", size = 16, angle = 0),
            axis.title.x  = element_text(color = "black", size = 16, angle = 0),
            axis.title.y  = element_text(color = "black", size = 16, angle = 90),
            legend.title  = element_text(color = "black", size  = 16),
            legend.text   = element_text(color = "black", size   = 16),
            axis.line.y = element_line(color = "black", linetype = "solid"),
            axis.line.x = element_line (color = "black",linetype = "solid"),
            panel.border = element_rect(linetype = "solid", size = 1.2,fill = NA) )

  } else {
    ggboxplot(mdat,
              x = x,
              y = y,
              fill = fill,
              bxp.errorbar = T,
              bxp.errorbar.width = 0.2,
              palette = "npg",
              add = add) +
      labs(title = title,
           subtitle = subtitle,
           caption =caption,
           x = xlab,
           y = ylab) +
      theme(plot.title    = element_text(color = "black", size   = 16, hjust = 0.5),
            plot.subtitle = element_text(color = "black", size   = 14,hjust = 0.5),
            plot.caption  = element_text(color = "black", size   = 16,face = "italic", hjust = 1),
            axis.text.x   = element_text(color = "black", size = 16, angle = 0),
            axis.text.y   = element_text(color = "black", size = 16, angle = 0),
            axis.title.x  = element_text(color = "black", size = 16, angle = 0),
            axis.title.y  = element_text(color = "black", size = 16, angle = 90),
            legend.title  = element_text(color = "black", size  = 16),
            legend.text   = element_text(color = "black", size   = 16),
            axis.line.y = element_line(color = "black", linetype = "solid"),
            axis.line.x = element_line (color = "black",linetype = "solid"),
            panel.border = element_rect(linetype = "solid", size = 1.2,fill = NA) )

  }

}




if(datatype=="Gene_expression"|datatype=="miRNA"|datatype=="protein"|datatype=="Methylation"){

tmp1=as.data.frame(tmp1)
tmp1[,variable]=as.numeric(tmp1[,variable])
tmp1[,"group"]=factor(tmp1[,"group"])
write.table(tmp1,file = paste("out/boxplot",variable,tm,".txt",sep=""),sep = "\t",quote = F,row.names = F,col.names =F)
box=zh.sig.boxplot(mdat = tmp1,x = "group",y = variable,title = variable,
                   xlab = "tumor vs normal", ylab = "TPM",
                   sig = T)

fname=paste("out/",ca,"_boxplot_",variable,"-",tm,sep = "")
png(file=paste(fname,".png",sep = ""),width = 816*3,height = 641*3,res=300)
print(box)
dev.off()

pdf(file=paste(fname,".pdf",sep = ""),width = 10,height = 8,onefile = FALSE)
print(box)
dev.off()

tiff(file=paste(fname,".tiff",sep = ""),width = 816*3,height = 641*3,res = 300)
print(box)
dev.off()
}


