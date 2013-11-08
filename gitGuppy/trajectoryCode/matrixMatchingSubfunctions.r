#########################################
#### The following two functions are called in the neighbourhood ID part of matrixMatching()

# called from matrixMatching(). Identifies whether mismatches between the (permissive) pathMatrix and
# (restrictive) recipMatrix are justified. In cases where the most restrictive conditions aren't warranted, we consider
# more permissive matches.
identifyJoins=function(pathMatrix, recipMatrix, compAreaMatrix, unRequited){
  for(j in 1:length(unRequited)){
    z=unRequited[j]
    addiMatrix=recipMatrix
    addiMatrix[z,]=pathMatrix[z,]
    focusMatrix=t(t(addiMatrix)*addiMatrix[z,])
    
 
    withProblem=t(t(compAreaMatrix)+compAreaMatrix[z,]) #add the area of the problematic fly to the rest of the areas 
  #  focusMatrix=t(t(pathMatrix)*pathMatrix[z,]) #describe which cells we want to consider
    withProblem=withProblem*focusMatrix
    withoutProblem=compAreaMatrix*focusMatrix
    isItCloser=abs(withProblem-1)
    if(length(isItCloser[,1])==2){
      keeper=(isItCloser[-z,, drop=F]<abs(withoutProblem[-z,, drop=F]-1))*1
    }else {
      #keeper=colSums(isItCloser[-z,]<abs(withoutProblem[-z,]-1)*1) #does this work better?
      keeper=(abs(colSums(withProblem[-z,, drop=F])-1)<abs((colSums(withoutProblem[-z,, drop=F]))-1))*1
    }
    recipMatrix[z,]=keeper
    #recipMatrix[isItCloser<abs(withoutProblem-1)]=1
   # recipMatrix[z, c(colSums((isItCloser<abs(withoutProblem-1))*1)>0)]=1 #weird, I know, but just in case the z entry in question is greater than 0.5
    
  }
  return(recipMatrix)
}

# called from matrixMatching(). Identifies whether mismatches between the (permissive) pathMatrix and
# (restrictive) recipMatrix are justified. In cases where the most restrictive conditions aren't warranted, we consider
# more permissive matches.

identifyLeaves=function(pathMatrix, recipMatrix, compAreaMatrix, unRequited){
  compAreaMatrix=1/compAreaMatrix
  for(j in 1:length(unRequited)){
    z=unRequited[j]
    
    addiMatrix=recipMatrix
    addiMatrix[,z]=pathMatrix[,z]
    focusMatrix=addiMatrix*addiMatrix[,z]
    
    withProblem=compAreaMatrix+compAreaMatrix[,z] #add the problematic column to all other columns
   # focusMatrix=pathMatrix*pathMatrix[,z] 
    withProblem=withProblem*focusMatrix
    withoutProblem=compAreaMatrix*focusMatrix
    isItCloser=abs(withProblem-1)
    if(length(isItCloser[1,])==2){
      keeper=(isItCloser[,-z, drop=F]<abs(withoutProblem[,-z, drop=F]-1))*1
    }else {
      #keeper=rowSums(isItCloser[,-z]<abs(withoutProblem[,-z]-1)*1)
      keeper=(abs(rowSums(withProblem[,-z, drop=F])-1)<abs((rowSums(withoutProblem[,-z, drop=F]))-1))*1
    }
    recipMatrix[,z]=keeper
 
   # recipMatrix[isItCloser<abs(withoutProblem-1)]=1
   # recipMatrix[c(rowSums((isItCloser<abs(withoutProblem-1))*1)>0),z]=1 #weird, I know, but just in case the z entry in question is greater than 0.5

  }
  return(recipMatrix)
}


###################################################################
################ These next functions represent matches between combinations 
####### of blobs of varying complexity.


# called from matrixMatching()
# calls computeRunningStats()
# This gets called when there is a single, unambiguous, best match between blobs. A single next-blob inherits the fly (or all the flies)
# from a previous blob
matchPair=function(dfNextBlob, dfPrevBlob){
  #the following two diagnostics may help us identify when we're assigning the wrong number of flies to a blob
  areaDeviance=sum(dfPrevBlob$flyArea)-dfNextBlob$blobArea
  propDeviance=abs(areaDeviance)/dfNextBlob$blobArea
  twoOne=dfNextBlob$nFlies[1]==1 
  
  if(dfPrevBlob$nFlies[1]>1){
    dfNextBlob=dfNextBlob[rep(1, dfPrevBlob$nFlies[1]),]
    dfNextBlob$nFlies=dfPrevBlob$nFlies[1]
    dfNextBlob$framesDeviance=dfPrevBlob$framesDeviance
    if(twoOne){
      dfNextBlob$framesDeviance=dfNextBlob$framesDeviance+1
    }
  }else{
    dfNextBlob$nFlies=1
    dfNextBlob$framesDeviance=dfPrevBlob$framesDeviance
  }
  
  dfNextBlob$flyID=dfPrevBlob$flyID
  dfNextBlob$areaDeviance=areaDeviance
  dfNextBlob=computeRunningStats(dfNextBlob, dfPrevBlob)
  
  if((propDeviance>0.1)){
    dfNextBlob$framesDeviance=dfNextBlob$framesDeviance+1
  }
  
  return(dfNextBlob)
}

# called from within matchPair(), and computes update stats in cases of simple blob-to-blob matches.
# corresponds to mixoMatchosis(), used for more complex matches
#This takes a pair of blobs, which match with no problems, 
# and update the next frame flies. If the blobs have more than 
# a single fly, running statistics like area and stuff don't get 
# recalculated.
computeRunningStats=function(dfNextBlob, dfPrevBlob){
  # we won't update the fly stats if we think the fly is coblobulated
  areaDeviance=sum(dfPrevBlob$flyArea)-dfNextBlob$blobArea
  
  if((dfNextBlob$nFlies[1]==1) & !is.na(dfPrevBlob$deltaX[1])){
    if((abs(areaDeviance[1]) < 100) | dfNextBlob$framesDeviance>5){
       #I want to be a bit careful about updating flies too quickly, if they change area rapidly. But not too careful
      dfNextBlob$deltaX=dfPrevBlob$deltaX*0.7 + (dfNextBlob$blobX-dfPrevBlob$blobX)*0.3
      dfNextBlob$deltaY=dfPrevBlob$deltaY*0.7 + (dfNextBlob$blobY-dfPrevBlob$blobY)*0.3
      dfNextBlob$flySpeed=dfPrevBlob$flySpeed*0.8 + 0.2*(dfNextBlob$deltaX*dfNextBlob$deltaX + dfNextBlob$deltaY*dfNextBlob$deltaY)^0.5
      dfNextBlob$flyArea=dfPrevBlob$flyArea*0.8 + dfNextBlob$blobArea*0.2
      dfNextBlob$flyAspect=dfPrevBlob$flyAspect*0.9 + (dfNextBlob$blobB/dfPrevBlob$blobA)*0.1
      dfNextBlob$areaDeviance=areaDeviance
      dfNextBlob$framesDeviance=dfNextBlob$framesDeviance*0.9
    }else{
      dfNextBlob[,c("deltaX", "deltaY", "flySpeed", "flyArea", "flyAspect")]=dfPrevBlob[,c("deltaX", "deltaY", "flySpeed", "flyArea", "flyAspect")]
      dfNextBlob$areaDeviance=areaDeviance    
    }
  }else if(is.na(dfPrevBlob$deltaX[1])){
    #this should only happen on frame 2 of any trajectory
    dfNextBlob$deltaX=dfNextBlob$blobX-dfPrevBlob$blobX
    dfNextBlob$deltaY=dfNextBlob$blobY-dfPrevBlob$blobY
    dfNextBlob$flySpeed=(dfNextBlob$deltaX*dfNextBlob$deltaX + dfNextBlob$deltaY*dfNextBlob$deltaY)^0.5
    dfNextBlob$flyArea= dfPrevBlob$flyArea*0.8 + dfNextBlob$blobArea*0.2
    dfNextBlob$flyAspect=dfPrevBlob$flyAspect*0.8 + (dfNextBlob$blobB/dfPrevBlob$blobA)*0.2
    dfNextBlob$areaDeviance=areaDeviance
    dfNextBlob$framesDeviance=dfNextBlob$framesDeviance*0.9
  }else{
    if((max(dfNextBlob$framesDeviance)<8) | (max(areaDeviance)<0)){
      #this should only happen if there are multiple flies, all of which have entries for the values we're carrying through
      dfNextBlob[,c("deltaX", "deltaY", "flySpeed", "flyArea", "flyAspect")]=dfPrevBlob[,c("deltaX", "deltaY", "flySpeed", "flyArea", "flyAspect")]
      dfNextBlob$areaDeviance=areaDeviance
     }else{
       # I don't really arbitrarily want to send some fly off into the ether, but, if it gets to this point, every opportunity has been given to sort out 
       # where the fly could have gone. Something is awry, and keeping it longer will only make things worse
       dfNextBlob[,c("deltaX", "deltaY", "flySpeed", "flyArea", "flyAspect")]=dfPrevBlob[,c("deltaX", "deltaY", "flySpeed", "flyArea", "flyAspect")]
       dfNextBlob$areaDeviance=areaDeviance      
     }
  }
  dfNextBlob$flyFrame=dfPrevBlob$flyFrame+1      
  return(dfNextBlob)
}


# called from matrixMatching()
# an attempt to pull out simple, good matches from a messy bunch of nearby blobs 
# using both distance and area data and a hungarian algorithm approach before going to the full all-combinations
# assignment algorithm.
identifyBest=function(groupNextIndex, groupPrevIndex, pathMatrix, compAreaMatrix, mDistMatrix){
  euclidDistMatrix=mDistMatrix[groupNextIndex, groupPrevIndex, drop=F]
  areaDistMatrix=compAreaMatrix[groupNextIndex, groupPrevIndex, drop=F]
  areaDistMatrix=(log(areaDistMatrix)*10)^2 #this is probably not strictly necessary, but if I ever want to combine the area and euclidean 
					#distances into a single measure, we'll need to multiply by a constant to get them into something like
					#an equivalent scale
					###This is ugly
  areaDistMatrix=areaDistMatrix^2
 ##Pretty lenient threshold. ~15% change in area.

  
  unSquare=length(groupNextIndex)-length(groupPrevIndex)
  if(unSquare>0){
    enSquare=matrix(rep(999999, unSquare*length(groupNextIndex)), nrow=length(groupNextIndex), ncol=unSquare)
    areaDistMatrix=cbind(areaDistMatrix, enSquare)
    euclidDistMatrix=cbind(euclidDistMatrix, enSquare)
  }else if(unSquare<0){
    enSquare=matrix(rep(999999, ((0-unSquare)*length(groupPrevIndex))), nrow=(0-unSquare), ncol=length(groupPrevIndex)) 
    areaDistMatrix=rbind(areaDistMatrix, enSquare)
    euclidDistMatrix=rbind(euclidDistMatrix, enSquare)
  }
  
  threshCutoffMatrix=areaDistMatrix<2 
  
  areaHungarian=lp.assign(areaDistMatrix)
  euclidTheHungarian=lp.assign(euclidDistMatrix)
  
  hungarianAgreement=(areaHungarian$solution>0)&(euclidTheHungarian$solution>0)
  hungarianAgreement=hungarianAgreement&threshCutoffMatrix
  if(unSquare !=0){
    hungarianAgreement=hungarianAgreement[1:length(groupNextIndex), 1:length(groupPrevIndex)]
  }
  prevAgreeable=rowSums(t(t(hungarianAgreement)*groupPrevIndex))
  bestUnambiguousMatches=cbind(groupNextIndex, prevAgreeable)
  
  bestUnambiguousMatches=bestUnambiguousMatches[c(prevAgreeable>0),, drop=F]
 
  return(bestUnambiguousMatches)
}

# called from matrixMatching()
# calls out to subfunctions splitComplex(), tractablyComplex() and perplexComplex() 
# (complexMatchingSubfunctions.r)
# in order of difficulty.
complexMatch=function(dfNextFrame, dfPrevFrame, groupNextIndex, groupPrevIndex){
   
    muddledFlies=dfPrevFrame[dfPrevFrame$blobIndex %in% groupPrevIndex,]
    muddledBlobs=dfNextFrame[dfNextFrame$blobIndex %in% groupNextIndex,]
    nFlies=nrow(muddledFlies)
    muddledFlyIndex=c(1:nFlies)
    muddledFlies=cbind(muddledFlies, muddledFlyIndex)
    nBlobs=nrow(muddledBlobs)
    
    
    if((nFlies==1) && (nBlobs >1)){
      diffMat=splitComplex(muddledFlies, muddledBlobs,nFlies, muddledFlyIndex, nBlobs)
    }else if(nFlies*nBlobs<16){
      diffMat=tractablyComplex(muddledFlies, muddledBlobs,nFlies, muddledFlyIndex, nBlobs)
    }else{
      diffMat=perplexComplex(muddledFlies, muddledBlobs,nFlies, muddledFlyIndex, nBlobs)
   
    }
        
      #so sorry for the egregiously nested syntax :(
      flyList=unlist(lapply(1:nBlobs, function(k) {muddledFlies$flyID[diffMat[(diffMat[,"blobIndex"]==muddledBlobs$blobIndex[k]),3]>0]}))
      
      blobuList=diffMat[diffMat[,3]>0, "blobIndex", drop=F]
      propArea=diffMat[diffMat[,3]>0, 3, drop=F]
      flyApportion=data.frame(blobuList, flyList, propArea)
      
      #bestMuddledMatches=mixoMatchosis(dfNextBlobs, dfPrevFlies)
      euNique=unique(flyApportion$blobIndex)
      #dear whomever has to parse the following line: I am truly sorry from the bottom of my heart.
      bestMuddledMatches=do.call(rbind, lapply(1:length(euNique),  function (k) 
	    {mixoMatchosis(dfNextFrame[(dfNextFrame$blobIndex==euNique[k]),], dfPrevFrame[(dfPrevFrame$flyID %in% flyApportion$flyList[flyApportion$blobIndex==euNique[k]]),])}))
      bestMuddledMatches$flyArea=bestMuddledMatches$flyArea*flyApportion[,3]
      bestMuddledMatches$areaDeviance=unlist( lapply(1:length(euNique), function(k) 
	    {sum(bestMuddledMatches[bestMuddledMatches$blobIndex==euNique[k], "flyArea",drop=F])-bestMuddledMatches[bestMuddledMatches$blobIndex==euNique[k],,drop=F]$blobArea}))
      vSplitization=duplicated(bestMuddledMatches$flyID)
      if(sum(vSplitization>0)){
	bestMuddledMatches$flyID[vSplitization]="new"
	bestMuddledMatches=createNewFly(bestMuddledMatches)
	}
    bestMuddledMatches$framesDeviance=bestMuddledMatches$framesDeviance/2
    return(bestMuddledMatches)
}

