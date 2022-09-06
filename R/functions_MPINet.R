##########################################################################################################################################
###From MPInet package, identifypathway.R
##########################################################################################################################################
identifypathway<-function(componentList,PSS,pathType="KEGG",method="MPINet",
                          weightnum=6,backgroundcid=Getenvir("getBackground"),annlim=1,bglim=6,order="pvalue",decreasing=FALSE){


  Veg<-Getenvir("getnodeseq")
  pathType<-unique(pathType)
  if((length(pathType)==1)&(length(pathType[which(pathType=="KEGG")])>0))
  {
    pathList<-Getenvir("getpathList('KEGG')")

    if(method=="MPINet")
    {
      componentList<-unique(intersect(componentList,backgroundcid))
      message(c("your input componentList have ",length(componentList)," components in background"))
      compoundListinnet<-intersect(Veg,componentList)
      message(c("your input componentList have ",length(compoundListinnet)," components in network"))

      annList<-list()
      N<-length(unique(backgroundcid)) ###number of background compounds

      if(length(pathList[,1])>0)
      {
        for(i in 1:length(pathList[,1]))
        {
          ann<-list(pathwayId=character(),pathwayName="not known",annComponentList=character(),annComponentNumber=0,
                    annBgComponentList=character(),annBgNumber=0,componentNumber=0,bgNumber=0,pvalue=1,fdr=1,
                    InWeight=0,weight=1,anncompinNetworkNum=0,anncompinNetworkList=character(),riskcompinNetworkNum=0,riskcompinNetworkList=character())

          ann$pathwayId<-pathList[i,1]  ###pathwayId
          ann$pathwayName<-pathList[i,2]###pathwayName
          ann$annBgComponentList<-pathList[i,3]###annBgComponentList
          ann$bgNumber<-N                      ### bgNumber
          ann$componentNumber<-length(componentList) ###componentNumber

          tpathmetabolite<-strsplit(pathList[i,"annBgComponentList"],';')[[1]]####all metabolites in pathway i
          m1<-length(unique(unlist(tpathmetabolite)))
          inter<-unique(intersect(as.character(componentList),unlist(tpathmetabolite)))
          x<-length(unique(inter))

          ann$annBgNumber<-m1                   ####annBgNumber
          ann$annComponentNumber<-x             ####annComponentNumber
          ann$annComponentList<-inter           ####annComponentList

          #######################################calcualte the anncompinNetworkNum,anncompinNetworkList,riskcompinNetworkNum,riskcompinNetworkList of pathway i ##########################
          inter_2<-intersect(unique(inter),Veg)######## differential metabolites annotate in pathway i and network
          riskinnet<-intersect(componentList,Veg)#######differential metabolites in network

          ann$anncompinNetworkNum<-length(inter_2) ###anncompinNetworkNum
          ann$anncompinNetworkList<-inter_2        ###anncompinNetworkList
          ann$riskcompinNetworkNum<-length(riskinnet) ###riskcompinNetworkNum
          ann$riskcompinNetworkList<-riskinnet        ###riskcompinNetworkList

          #######################calculate pathway Inweight################

          strinner<-0

          metaboliteinpath<-as.character(unlist(tpathmetabolite))
          metapathinnet<-intersect(Veg,unlist(tpathmetabolite))

          if(length(metapathinnet)>0)
          {
            pathw<-sum(PSS[metapathinnet,"pss"])
            strinner<-pathw/length(metaboliteinpath)


          }







          ann$InWeight<-strinner ###InWeight




          annList[[i]]<-ann
        }#####for i in 1:length(pathList[,1])


        #######################calculate pathway weight################
        annList<-annList[sapply(annList,function(x) x$annComponentNumber>(annlim-1))]

        if(length(annList)>0)
        {
          annList<-annList[sapply(annList,function(x) x$annBgNumber>(bglim-1))]

          if(length(annList)>0)
          {
            InWeightList<-c()
            for(j in (1:length(annList)))
            {
              InWeightList<-c(InWeightList,annList[[j]]$InWeight)

            }
            meanweight<-mean(InWeightList)

            for(j in (1:length(annList)))
            {

              w<-annList[[j]]$InWeight/meanweight

              annList[[j]]$weight<-w^weightnum

            }

            ############################### calculate Pvalue###########################


            for(i in (1:length(annList)))
            {

              tv<-as.numeric(annList[[i]]["weight"])

              if(tv>0)
              {

                pvalue<- BiasedUrn::pWNCHypergeo((as.numeric(annList[[i]]["annComponentNumber"])-1),as.numeric(annList[[i]]["annBgNumber"]),as.numeric(annList[[i]]["bgNumber"])-as.numeric(annList[[i]]["annBgNumber"]),as.numeric(annList[[i]]["componentNumber"]),tv,precision=1E-100,lower.tail=FALSE)
                annList[[i]]["pvalue"]<-pvalue

              }

            }####### for(i in (1:length(annList)))

          }###########if(length(annList)>0)




        }###############if(length(annList)>0)

      }####if(length(pathList[,1])>0)




    }#####if method=="MPINet"

    if(method=="Hyper")
    {
      componentList<-unique(intersect(componentList,backgroundcid))
      message(c("your input componentList have ",length(componentList)," components in background"))
      annList<-list()
      N<-length(unique(backgroundcid)) ### number of background compounds

      if(length(pathList[,1])>0)
      {
        for(i in 1:length(pathList[,1]))
        {
          ann<-list(pathwayId=character(),pathwayName="not known",annComponentList=character(),annComponentNumber=0,
                    annBgComponentList=character(),annBgNumber=0,componentNumber=0,bgNumber=0,pvalue=1,fdr=1)

          ann$pathwayId<-pathList[i,1]  ###pathwayId
          ann$pathwayName<-pathList[i,2]###pathwayName
          ann$annBgComponentList<-pathList[i,3]###annBgComponentList
          ann$bgNumber<-N                      ### bgNumber
          ann$componentNumber<-length(componentList) ###componentNumber

          tpathmetabolite<-strsplit(pathList[i,"annBgComponentList"],';')[[1]]####all metabolites in pathway i
          m1<-length(unique(unlist(tpathmetabolite)))
          inter<-unique(intersect(as.character(componentList),unlist(tpathmetabolite)))
          x<-length(unique(inter))

          ann$annBgNumber<-m1                   ####annBgNumber
          ann$annComponentNumber<-x             ####annComponentNumber
          ann$annComponentList<-inter           ####annComponentList


          annList[[i]]<-ann
        }#####for i in 1:length(pathList[,1])



        ###############################calculate Pvalue###########################
        annList<-annList[sapply(annList,function(x) x$annComponentNumber>(annlim-1))]
        if(length(annList)>0)
        {

          annList<-annList[sapply(annList,function(x) x$annBgNumber>(bglim-1))]
          if(length(annList)>0)
          {
            for(i in (1:length(annList)))
            {



              pvalue<- BiasedUrn::pWNCHypergeo((as.numeric(annList[[i]]["annComponentNumber"])-1),as.numeric(annList[[i]]["annBgNumber"]),as.numeric(annList[[i]]["bgNumber"])-as.numeric(annList[[i]]["annBgNumber"]),as.numeric(annList[[i]]["componentNumber"]),1,precision=1E-100,lower.tail=FALSE)
              annList[[i]]["pvalue"]<-pvalue



            }####### for(i in (1:length(annList)))

          }########if(length(annList)>0)

        }########if(length(annList)>0)

      }####if(length(pathList[,1])>0)


    }######if method=="Hyper"





  }############# if (pathtype=="KEGG")










  else
  {

    if((length(pathType)==1)&(length(pathType[which(pathType=="consensusPath")])>0))
    {

      pathList<-Getenvir("consensusPath")
    }

    else
    {
      pathListall<-Getenvir("consensusPath")
      pathList<-pathListall[which(as.character(pathListall[,2])%in%pathType),]
    }
    if(method=="MPINet")
    {
      componentList<-unique(intersect(componentList,backgroundcid))
      message(c("your input componentList have ",length(componentList)," components in background"))
      compoundListinnet<-intersect(Veg,componentList)
      message(c("your input componentList have ",length(compoundListinnet)," components in network"))

      annList<-list()
      N<-length(unique(backgroundcid)) ###number of background compounds

      if(length(pathList[,1])>0)
      {
        for(i in 1:length(pathList[,1]))
        {
          ann<-list(pathwayName="not known",pathsource=character(),annComponentList=character(),annComponentNumber=0,
                    annBgComponentList=character(),annBgNumber=0,componentNumber=0,bgNumber=0,pvalue=1,fdr=1,
                    InWeight=0,weight=1,anncompinNetworkNum=0,anncompinNetworkList=character(),riskcompinNetworkNum=0,riskcompinNetworkList=character())

          ann$pathsource<-pathList[i,2]  ###pathwayId
          ann$pathwayName<-pathList[i,1]###pathwayName
          ann$annBgComponentList<-pathList[i,3]###annBgComponentList
          ann$bgNumber<-N                      ### bgNumber
          ann$componentNumber<-length(componentList) ###componentNumber

          tpathmetabolite<-strsplit(pathList[i,"metabolites"],';')[[1]]####all metabolites in pathway i
          m1<-length(unique(unlist(tpathmetabolite)))
          inter<-unique(intersect(as.character(componentList),unlist(tpathmetabolite)))
          x<-length(unique(inter))

          ann$annBgNumber<-m1                   ####annBgNumber
          ann$annComponentNumber<-x             ####annComponentNumber
          ann$annComponentList<-inter           ####annComponentList

          #######################################calculate anncompinNetworkNum,anncompinNetworkList,riskcompinNetworkNum,riskcompinNetworkList of pathway i##########################
          inter_2<-intersect(unique(inter),Veg)########differential metabolites annotate in pathway i and network
          riskinnet<-intersect(componentList,Veg)#######differential metabolites in network

          ann$anncompinNetworkNum<-length(inter_2) ###anncompinNetworkNum
          ann$anncompinNetworkList<-inter_2        ###anncompinNetworkList
          ann$riskcompinNetworkNum<-length(riskinnet) ###riskcompinNetworkNum
          ann$riskcompinNetworkList<-riskinnet        ###riskcompinNetworkList

          #######################calculate pathway InWeigtht################

          strinner<-0

          metaboliteinpath<-as.character(unlist(tpathmetabolite))
          metapathinnet<-intersect(Veg,unlist(tpathmetabolite))

          if(length(metapathinnet)>0)
          {
            pathw<-sum(PSS[metapathinnet,"pss"])
            strinner<-pathw/length(metaboliteinpath)


          }







          ann$InWeight<-strinner ###InWeight




          annList[[i]]<-ann
        }#####for i in 1:length(pathList[,1])

        annList<-annList[sapply(annList,function(x) x$annComponentNumber>(annlim-1))]

        if(length(annList)>0)
        {
          annList<-annList[sapply(annList,function(x) x$annBgNumber>(bglim-1))]

          if(length(annList)>0)
          {
            InWeightList<-c()
            for(j in (1:length(annList)))
            {
              InWeightList<-c(InWeightList,annList[[j]]$InWeight)

            }
            meanweight<-mean(InWeightList)

            for(j in (1:length(annList)))
            {

              w<-annList[[j]]$InWeight/meanweight

              annList[[j]]$weight<-w^weightnum







            }

            ###############################calculate Pvalue###########################


            for(i in (1:length(annList)))
            {

              tv<-as.numeric(annList[[i]]["weight"])

              if(tv>0)
              {

                pvalue<- BiasedUrn::pWNCHypergeo((as.numeric(annList[[i]]["annComponentNumber"])-1),as.numeric(annList[[i]]["annBgNumber"]),as.numeric(annList[[i]]["bgNumber"])-as.numeric(annList[[i]]["annBgNumber"]),as.numeric(annList[[i]]["componentNumber"]),tv,precision=1E-100,lower.tail=FALSE)
                annList[[i]]["pvalue"]<-pvalue

              }

            }####### for(i in (1:length(annList)))

          }########if(length(annList)>0)
        }########if(length(annList)>0)

      }####if(length(pathList[,1])>0)
    }#####if method=="MPINet"

    if(method=="Hyper")
    {
      componentList<-unique(intersect(componentList,backgroundcid))
      message(c("your input componentList have ",length(componentList)," components in background"))
      annList<-list()
      N<-length(unique(backgroundcid)) ###number of background compounds

      if(length(pathList[,1])>0)
      {
        for(i in 1:length(pathList[,1]))
        {
          ann<-list(pathwayName="not known",pathsource=character(),annComponentList=character(),annComponentNumber=0,
                    annBgComponentList=character(),annBgNumber=0,componentNumber=0,bgNumber=0,pvalue=1,fdr=1)

          ann$pathsource<-pathList[i,2]  ###pathwayId
          ann$pathwayName<-pathList[i,1]###pathwayName
          ann$annBgComponentList<-pathList[i,3]###annBgComponentList
          ann$bgNumber<-N                      ### bgNumber
          ann$componentNumber<-length(componentList) ###componentNumber

          tpathmetabolite<-strsplit(pathList[i,"metabolites"],';')[[1]]#### all metabolites in pathway i
          m1<-length(unique(unlist(tpathmetabolite)))
          inter<-unique(intersect(as.character(componentList),unlist(tpathmetabolite)))
          x<-length(unique(inter))

          ann$annBgNumber<-m1                   ####annBgNumber
          ann$annComponentNumber<-x             ####annComponentNumber
          ann$annComponentList<-inter           ####annComponentList


          annList[[i]]<-ann
        }#####for i in 1:length(pathList[,1])



        ###############################calculate Pvalue###########################
        annList<-annList[sapply(annList,function(x) x$annComponentNumber>(annlim-1))]
        if(length(annList)>0)
        {
          annList<-annList[sapply(annList,function(x) x$annBgNumber>(bglim-1))]
          if(length(annList)>0)
          {
            for(i in (1:length(annList)))
            {



              pvalue<- BiasedUrn::pWNCHypergeo((as.numeric(annList[[i]]["annComponentNumber"])-1),as.numeric(annList[[i]]["annBgNumber"]),as.numeric(annList[[i]]["bgNumber"])-as.numeric(annList[[i]]["annBgNumber"]),as.numeric(annList[[i]]["componentNumber"]),1,precision=1E-100,lower.tail=FALSE)
              annList[[i]]["pvalue"]<-pvalue



            }####### for(i in (1:length(annList)))

          }###########if(length(annList)>0)
        }####if(length(annList)>0)
      }####if(length(pathList[,1])>0)


    }######if method=="Hyper"


  }###########else


  ###############################calculate fdr#####################################
  if(length(annList)>0)
  {
    p_value<-sapply(annList,function(x) return(x$pvalue))
    #fdrtool.List<-fdrtool(p_value,statistic="pvalue",plot=FALSE,verbose=FALSE)

    #print(fdrtool.List$qval)
    #for(i in seq(annList)){
    #   annList[[i]]$qvalue<-fdrtool.List$qval[i]
    #	annList[[i]]$lfdr<-fdrtool.List$lfdr[i]
    #}
    fdr.List<-fdr.est(p_value)
    for(i in seq(annList))
    {
      annList[[i]]$fdr<-fdr.List[i]
    }
    #names(annList)<-sapply(graphList,function(x) x$number)

    annList<-annList[order(sapply(annList,function(x) x[[order]]),decreasing=decreasing)]
  }

  return(annList)


}######function identifypathway


#####################################################################
fdr.est<-function(p)
{
  m <- length(ind <- which(!is.na(p)))
  fdr <- rep(NA, length(p))
  stat <- cbind(1:length(p), p, fdr)
  stat[ind, 3] <- unlist(lapply(stat[ind, 2], function(x) {
    c <- length(which(stat[ind, 2] <= x))
    m * x/c
  }))
  stat[ind, ] <- stat[ind, ][order(stat[, 2], decreasing = TRUE),
  ]
  stat[ind, 3] <- cummin(stat[ind, 3])
  fdr <- stat[order(stat[, 1]), 3]
  fdr[which(fdr > 1)] <- 1
  return(fdr)
}

##########################################################################################################################################
###From MPInet package, getPSS.R
##########################################################################################################################################

getPSS<-function(riskmeta,plot=TRUE,binsize=400){


  nodeseq<-Getenvir("getnodeseq")
  strnew<-Getenvir("getStr")

  riskmetainnet<-intersect(nodeseq,riskmeta)
  riskmetainnet<-unique(riskmetainnet)



  infor<-list()
  for(i in 1:length(strnew[,1]))
  {

    infor_i<-list(nodeId="000",meanvalue=0,diffnode=0)

    infor_i$nodeId<-strnew[i,1]
    infor_i$meanvalue<-strnew[i,2]

    if(strnew[i,1]%in%riskmetainnet)
    {
      infor_i$diffnode<-1
    }


    infor[[i]]<-infor_i

  }
  nodeId<-sapply(infor,function(x) x$nodeId)
  meanvalue<-sapply(infor,function(x) x$meanvalue)
  diffnode<-sapply(infor,function(x) x$diffnode)

  meanstr<-data.frame(nodeId=nodeId,meanvalue=meanvalue,diffnode=diffnode)
  meanstr<-meanstr[order(meanstr[,2],decreasing = FALSE),]

  pss<-performpcls(meanstr[,"meanvalue"],meanstr[,"diffnode"])
  CGNB<-1-pss
  re_out<-data.frame(riskmeta=meanstr[,"diffnode"],meanstrvalue=meanstr[,"meanvalue"],pss=pss,CGNB=CGNB,stringsAsFactors=FALSE)
  rownames(re_out)=meanstr[,"nodeId"]






  if(plot)
  {
    splitter=ceiling(1:length(re_out$riskmeta)/binsize)
    de=sapply(split(re_out$riskmeta,splitter),mean)
    binlen=sapply(split(as.numeric(re_out$meanstrvalue),splitter),mean)
    plot(binlen,de,xlab=paste("Bias Data in ",binsize," metabolite bins",sep=""),ylab="Proportion risk metabolite")
    lines(re_out$meanstrvalue,re_out$pss,col=3,lwd=2)

  }#####plot


  return(re_out)

}###########function()

##########################################################################################################################################
###From MPInet package, getEnvironmentData.R
##########################################################################################################################################

initializeMPINet<-function(){
  utils::data("MPINetData",package="famat")##package="MPINet" => package="FAMAT"
}

Getenvir<-function(envirdata){

  if(!exists("MPINetData")) initializeMPINet()
  return(get(envirdata,envir=MPINetData))

}

##########################################################################################################################################
###From MPInet package, performpcls.R
##########################################################################################################################################

performpcls<-function(x, y, nKnots = 6){

  f.ug <- gam(y ~ s(x, k = nKnots, bs = "cr"),family=binomial())
  dat <- data.frame(x = x, y = y)
  sm <- smoothCon(s(x, k = nKnots, bs = "cr"), dat, knots = NULL)[[1]]
  if (length(sm$xp) < nKnots)
    warning("Few than 6 nKnots were specified!\n Possible inaccurate fitting!\n")
  F <- mono.con(sm$xp,TRUE)

  G<-list(X=sm$X,C=matrix(0,0,0),sp=f.ug$sp,p=sm$xp,y=y,w=y*0+1)
  G$Ain<-F$A;G$bin<-F$b;G$S<-sm$S;G$off<-0
  p <- pcls(G)
  fv <- Predict.matrix(sm, data.frame(x = x)) %*% p
  fv <- as.vector(fv)
  #If pcls still fails for some reason and returns negative values (or values that are so low that they'll have an effective relative weight of Inf
  #then we need to add a little bit to every weight to make things non-negative, the choice of minimum weight is somewhat arbitrary and serves only to
  #ensure that we don't have positive weights that will be infinitely prefered as a ratio.
  lower_bound=10^-3
  if(min(fv)<lower_bound)
    fv=fv-min(fv)+lower_bound
  fv2<-as.matrix(fv)
  return(fv2)

}######function
