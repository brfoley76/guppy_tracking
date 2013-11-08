
# called from complexMatch()
# This function takes a single blob which has unambiguously split into several sub-blobs, and splits the flies up 
# (possibly creating new flies by mitosis)

splitComplex=function(muddledFlies, muddledBlobs,nFlies, muddledFlyIndex, nBlobs){
  yesNo=c(0,1)
  #I create a matrix with all possible combinations of 1 and 0 for flies within blobs
  mCombinatorium=lapply(numeric(nFlies*nBlobs), function(x) yesNo)
  mCombinatorium=as.matrix(expand.grid(mCombinatorium))
  #these indices let me keep track of blobs and flies in the assortment matrix. It's kind of an ugly way to create them, though!
  blobIndex=matrix(rep(muddledBlobs$blobIndex, nFlies), nrow=nBlobs, ncol=nFlies)
  blobIndex=as.vector(t(blobIndex))
  flyIndex=matrix(rep(muddledFlyIndex, nBlobs), nrow=nFlies, ncol=nBlobs)
  flyIndex=as.vector(flyIndex)

  mCombinatorium=cbind(blobIndex, flyIndex, t(mCombinatorium))
  
  mSplitorium=unlist(lapply(1:nrow(mCombinatorium),  function (k) {muddledBlobs[muddledBlobs[,"blobIndex"]==mCombinatorium[k,"blobIndex"],"blobArea",drop=F]}))
  mSplitorium= mSplitorium*mCombinatorium
  mSplitorium[,1:2]=mCombinatorium[,1:2]
  
  mSplitorium=do.call(rbind, lapply(1:nFlies,  function (k) {t(t(mSplitorium[mSplitorium[,"flyIndex"]==flyIndex[k],,drop=F])/colSums(mSplitorium[mSplitorium[,"flyIndex"]==flyIndex[k],,drop=F]))}))
  mSplitorium[,1:2]=mCombinatorium[,1:2]
  mSplitorium[is.na(mSplitorium)]=0
 
  noDuples=do.call(rbind, lapply(1:nFlies,  function (k) {colSums(mCombinatorium[mCombinatorium[,"flyIndex"]==flyIndex[k],,drop=F])}))
  droppedFlyIndex=noDuples==0
  #noDuples=noDuples<2


  #retainers=colSums(noDuples)==nFlies
  #retainers[1:2]=TRUE
  #droppedFlyIndex=droppedFlyIndex[,retainers, drop=F]
  #mCombinatorium=mCombinatorium[,retainers, drop=F]

  #here I calculate the predicted area of each blob, given each possible combination of flies. 
  #I pick the combination that minimizes the difference from expected

  flyAreaVec=muddledFlies$flyArea
  targetAreaVec=muddledBlobs$blobArea

  flyLeavePenalty=colSums(droppedFlyIndex*flyAreaVec)
  #flyLeavePenalty=t(as.data.frame(flyLeavePenalty))

  flyAreaCombos=cbind(mSplitorium[,1:2, drop=F], mSplitorium[,3:ncol(mSplitorium), drop=F]*flyAreaVec)
  #flyAreaCombos=cbind(mCombinatorium[,1:2, drop=F], mCombinatorium[,3:ncol(mCombinatorium), drop=F]*flyAreaVec)
  targetMatch=aggregate(flyAreaCombos, by=list(blobIndex), FUN="sum")
  #flyLeavePenalty=data.frame(targetMatch[,1,drop=F], flyLeavePenalty)

  targetMatch=abs(targetMatch-targetAreaVec)
  diffVec=colSums(targetMatch[c(-1, -2, -3)])+flyLeavePenalty[c(-1, -2)]
  minDiff=min(diffVec)

  diffVec=diffVec==minDiff
  diffVec=c(TRUE, TRUE, diffVec)
  diffMat=mSplitorium[,diffVec, drop=F]

  return(diffMat)
}


# This is an exhaustive fly-blob matching algorithm, which attempts to assign flies to blobs by minimizing the difference between combined-fly area
# and blob area across all blobs. The process is fine for up to 15 fly and blob combinations, but the assignment matrix is 2^(blob*fly) in size, 
# so the process gets exponentially slower as the number increases.
tractablyComplex=function(muddledFlies, muddledBlobs,nFlies, muddledFlyIndex, nBlobs){
  yesNo=c(0,1)
  #I create a matrix with all possible combinations of 1 and 0 for flies within blobs
  mCombinatorium=lapply(numeric(nFlies*nBlobs), function(x) yesNo)
  mCombinatorium=as.matrix(expand.grid(mCombinatorium))
  #these indices let me keep track of blobs and flies in the assortment matrix. It's kind of an ugly way to create them, though!
  blobIndex=matrix(rep(muddledBlobs$blobIndex, nFlies), nrow=nBlobs, ncol=nFlies)
  blobIndex=as.vector(t(blobIndex))
  flyIndex=matrix(rep(muddledFlyIndex, nBlobs), nrow=nFlies, ncol=nBlobs)
  flyIndex=as.vector(flyIndex)

  mCombinatorium=cbind(blobIndex, flyIndex, t(mCombinatorium))
  noDuples=do.call(rbind, lapply(1:nFlies,  function (k) {colSums(mCombinatorium[mCombinatorium[,"flyIndex"]==flyIndex[k],,drop=F])}))
  droppedFlyIndex=noDuples==0
  noDuples=noDuples>1
  retainers=colSums(noDuples)<1
  retainers[1:2]=TRUE
  
  mCombinatorium=mCombinatorium[,retainers, drop=F]
  droppedFlyIndex=droppedFlyIndex[,retainers, drop=F]
  
  #here I calculate the predicted area of each blob, given each possible combination of flies. 
  #I pick the combination that minimizes the difference from expected

  flyAreaVec=muddledFlies$flyArea
  targetAreaVec=muddledBlobs$blobArea

  flyLeavePenalty=colSums(droppedFlyIndex*flyAreaVec)
  flyLeavePenalty=t(as.data.frame(flyLeavePenalty))
  flyAreaCombos=cbind(mCombinatorium[,1:2, drop=F], mCombinatorium[,3:ncol(mCombinatorium), drop=F]*flyAreaVec)
  targetMatch=aggregate(flyAreaCombos, by=list(blobIndex), FUN="sum")
  #flyLeavePenalty=data.frame(targetMatch[,1,drop=F], flyLeavePenalty)

  targetMatch=abs(targetMatch-targetAreaVec)
  diffVec=colSums(targetMatch[c(-1, -2, -3)])+flyLeavePenalty[c(-1, -2)]
  minDiff=min(diffVec)

  diffVec=diffVec==minDiff
  diffVec=c(TRUE, TRUE, diffVec)
  diffMat=mCombinatorium[,diffVec, drop=F]
  diffMat=diffMat[,1:3,drop=F]

  return(diffMat)
}

# I try to clean up a really horrible superblob before implementing tractablyComplex()
perplexComplex=function(muddledFlies, muddledBlobs,nFlies, muddledFlyIndex, nBlobs){
  #print("GODDAMIT!")
  #The following is *terrible*, but I need to get it into the same format as in the other functions.
  blobIndex=matrix(rep(muddledBlobs$blobIndex, nFlies), nrow=nBlobs, ncol=nFlies)
  blobIndex=as.vector(t(blobIndex))
  flyIndex=matrix(rep(muddledFlyIndex, nBlobs), nrow=nFlies, ncol=nBlobs)
  flyIndex=as.vector(flyIndex)
  dummyMat=cbind(blobIndex, flyIndex, 0)

  diffMat=c()
  
  
  if(sum(muddledBlobs$nFlies==1)>0){
    diffMat=dealWithSingletons(muddledFlies, muddledBlobs,nFlies, muddledFlyIndex, nBlobs)
    colnames(diffMat)=c("blobIndex", "flyIndex","")
    muddledFlies=muddledFlies[!muddledFlies$muddledFlyIndex %in% diffMat[,2],,drop=F]
    muddledBlobs=muddledBlobs[!muddledBlobs$blobIndex %in% diffMat[,1],,drop=F]
    nFlies=nrow(muddledFlies)
    muddledFlyIndex=muddledFlies$muddledFlyIndex
    nBlobs=nrow(muddledBlobs)
  }
  if(nBlobs*nFlies<16 & nBlobs*nFlies>0){
    diffMat=rbind(diffMat, tractablyComplex(muddledFlies, muddledBlobs,nFlies, muddledFlyIndex, nBlobs))
  }else{
    #print("REALLY goddamit")
    #print(muddledFlies)
    diffMat=rbind(diffMat, tractablyComplex(muddledFlies, muddledBlobs,nFlies, muddledFlyIndex, nBlobs))
  }
  diffMat=diffMat[order(diffMat[,"blobIndex"]),]
  dummyMat=rbind(dummyMat, diffMat)
  dummyAgg=aggregate(dummyMat, by=list(dummyMat[,"flyIndex"],dummyMat[,"blobIndex"]), FUN="sum")

  diffMat=cbind(blobIndex, flyIndex, dummyAgg[,5])
  return(diffMat)
}

#called from perplexComplex(). An attempt to clear up the possible singletons in a complicated superblob
# so that the exhaustive matching doesn't take forever.
dealWithSingletons=function(muddledFlies, muddledBlobs,nFlies, muddledFlyIndex, nBlobs){
  singleTons=muddledBlobs[muddledBlobs$nFlies==1,,drop=F]
  #changing the following at the moment because otherwise it'll screw up my call to matchPair#
  muddledFlies$nFlies=1
  ##I need to do something here to triple check I don't split mating pairs##
  mSingleAreas=matrix(cbind(singleTons$blobIndex, singleTons$blobArea), nrow=nrow(singleTons), ncol=2)
  mFlyAreas=matrix(cbind(muddledFlies$muddledFlyIndex, muddledFlies$flyArea), nrow=nrow(muddledFlies), ncol=2) 
  mDifference=outer(mSingleAreas[,2], mFlyAreas[,2], "-")
  
  unSquare=ncol(mDifference)-nrow(mDifference)
  if(unSquare>0){
    enSquare=matrix(rep(999999, unSquare*ncol(mDifference)), nrow=unSquare, ncol=ncol(mDifference))
    areaDistMatrix=rbind(mDifference, enSquare)
  }else if(unSquare<0){
    enSquare=matrix(rep(999999, (abs(unSquare)*nrow(mDifference))), nrow=nrow(mDifference), ncol=abs(unSquare)) 
    areaDistMatrix=cbind(areaDistMatrix, enSquare)
  }
  
  areaHungarian=lp.assign(areaDistMatrix^2)
  singleSolution=(areaHungarian$solution[1:nrow(mDifference), 1:ncol(mDifference), drop=F])>0
  
  diffMat=do.call(rbind, lapply(1:nrow(singleTons), function(k) {cbind(singleTons[k,]$blobIndex, muddledFlyIndex[singleSolution[k,]], 1)}))
  return(diffMat)
}


# this function takes the complicated matches between a set of 
# blobs, and updates blob and fly data accordingly
# called from complexMatch() and findHomes()
mixoMatchosis=function(dfNextBlob, dfPrevFlies){
   areaDeviance=sum(dfPrevFlies$flyArea)-dfNextBlob$blobArea
   propDeviance=abs(areaDeviance)/dfNextBlob$blobArea
   nFlies=nrow(dfPrevFlies)
   dfNextBlob=dfNextBlob[rep(1, nFlies),]
   dfNextBlob$nFlies=nFlies
   dfNextBlob$flyFrame=dfPrevFlies$flyFrame+1
   dfNextBlob[,c("flyID", "deltaX", "deltaY", "flySpeed", "flyArea", "flyAspect", "framesDeviance")]=dfPrevFlies[,c("flyID", "deltaX", "deltaY", "flySpeed", "flyArea", "flyAspect", "framesDeviance")]
   #dfNextBlob$framesDeviance=dfNextBlob$framesDeviance*0.8
   dfNextBlob$framesDeviance=0  #If there is a complex match, it seems impossible to reliably keep track of framesDeviance
   dfNextBlob$areaDeviance=areaDeviance
   if(propDeviance>0.1){
    dfNextBlob$framesDeviance=dfNextBlob$framesDeviance+1
   }

   return(dfNextBlob)
}
