#R Function to make a .tag file that includes the landamrks of two specimens for registration purposes
#Will read in two .tag files and combine them into a two volume .tag file
#XFM.dir is the directory with the directory with the .tag files to combine
#BaseSpec is the name of the specimen that all others are registered to
#SpecName is the name of the specimen to be registered
#We assume 3D landmarks
combo.tag=function(XFM.dir,BaseSpec,SpecName){
  Base.LMs=read.table(paste(MNC.dir,paste(BaseSpec,".tag",sep=""),sep="/"), sep=c(" ",";"), skip=4,header=FALSE,stringsAsFactors=FALSE)[,2:4]
  Spec.LMs=read.table(paste(MNC.dir,paste(SpecName,".tag",sep=""),sep="/"), sep=c(" ",";"), skip=4,header=FALSE,stringsAsFactors=FALSE)[,2:4]
  n.Base.LM=nrow(Base.LMs)
  #Create Output File and add Header
  FileName=paste("Tag_",BaseSpec,"_to_",SpecName,".tag",sep="")
  cat("MNI Tag Point File\n",file=FileName,append=FALSE)
  cat("Volumes = 2;\n",file=FileName,append=TRUE)
  cat(paste("%VIO_Volume: ",BaseSpec,".mnc\n",sep=""),file=FileName,append=TRUE)
  cat(paste("%VIO_Volume: ",SpecName,".mnc\n",sep=""),file=FileName,append=TRUE)
  cat("\n",file=FileName,append=TRUE)
  cat("Points = \n",file=FileName,append=TRUE)
  
  #Because we must include numbers and characters in each line, we will cat() each line instead of making a matrix first
  for (i in 1:n.Base.LM) {
    cat(paste(Base.LMs[i,1],Base.LMs[i,2],Base.LMs[i,3],Spec.LMs[i,1],Spec.LMs[i,2],Spec.LMs[i,3],"\n",sep=" "),file=FileName,append=TRUE)
  }
}



#Read in the landmark file from which tags will be created
Source.dir="/media/bhlab/MVG/Init2/Tag/JD_Temp"
MNC.dir="/media/bhlab/MVG/Init2/Tag/JD_Temp"
BaseSpec="cl_awse115_pax68" # Target
Spec.List=read.table(paste(Source.dir,"spec_list_JD_rev2.txt",sep="/"),stringsAsFactors=FALSE)  
setwd(MNC.dir)
for(j in 1:nrow(Spec.List)){
  combo.tag(Source.dir,BaseSpec,Spec.List[j,1])
}





