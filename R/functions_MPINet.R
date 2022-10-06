##########################################################################################################################################
###From MPInet package, identifypathway.R
##########################################################################################################################################
identifypathway<-function(componentList,pathType="KEGG",method="MPINet",
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

          annList[[i]]<-ann
        }#####for i in 1:length(pathList[,1])

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

          annList[[i]]<-ann
        }#####for i in 1:length(pathList[,1])

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


      }####if(length(pathList[,1])>0)


    }######if method=="Hyper"


  }###########else

  return(annList)


}######function identifypathway



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


