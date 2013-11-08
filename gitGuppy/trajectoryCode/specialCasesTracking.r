# wendyLostBoys() and these subfunctions take a last stab at identifying and fixing problems 
# before moving on to the next frame
# called from matrixMatching()
#lostBoys is a dataframe of flies from dfPrevFrame that did not get assigned to any blob
# vDeviantBoys is a vector of the flies in dfNextFrame that show evidence of having been in an inappropriate blob across several frames
wendyLostBoys=function(matchedLines, vDeviantBoys, lostBoys){
  
  if(nrow(lostBoys)>0){
     matchedLines=findHomes(matchedLines, lostBoys)
     vDeviantBoys=matchedLines$framesDeviance>1
  }
  
  if(sum(vDeviantBoys)>0){
    matchedLines=dealWithSwaps(matchedLines, vDeviantBoys)
    }
  vDeviantBoys=matchedLines$framesDeviance>8
  if(sum(vDeviantBoys)>0){
    matchedLines=clearlySomeoneNeedsToLeave(matchedLines, vDeviantBoys)
  }

  return(matchedLines)
}

# called from wendyLostBoys()
# check to see if there are any unassigned flies from dfPrevFrame (lostBoys) that can be assigned to a blob and 
# improve fit
findHomes=function(matchedLines, lostBoys){
  blobIt=matchedLines[!duplicated(matchedLines$blobIndex),]
  blobItIndex=blobIt$blobIndex
  
  strangeFit=matrix(rep(blobIt$areaDeviance, nrow(lostBoys)), nrow=nrow(blobIt), ncol=nrow(lostBoys), byrow=F)
  changeFit=abs(t(t(cbind(strangeFit)+lostBoys$flyArea)))
  
  
  if(ncol(changeFit)==1){
    fosterHome=changeFit==min(changeFit)
  }else{
    unSquare=nrow(changeFit)-ncol(changeFit)

    if(unSquare>0){
      enSquare=matrix(rep(999999, unSquare*nrow(changeFit)), nrow=nrow(changeFit), ncol=unSquare)    
      findFit=cbind(changeFit, enSquare)
    }else if(unSquare<0){
      enSquare=matrix(rep(999999, ((0-unSquare)*ncol(changeFit))), nrow=(0-unSquare), ncol=ncol(changeFit)) 
      findFit=rbind(changeFit, enSquare)
    }else{
      findFit=changeFit
    }

    fosterHome=lp.assign(findFit)$solution[1:nrow(changeFit), 1:ncol(changeFit), drop=F]

  }
  fosterHome=fosterHome>0
  fosterCutoff=changeFit/lostBoys$flyArea
  fosterCutoff=fosterCutoff<0.33
  if(sum(fosterCutoff[fosterHome])>0){
    fosterHome=t(t(fosterHome)*blobItIndex)
    fosterAddress=colSums(fosterHome*fosterCutoff)
    lostBoys$blobIndex=fosterAddress
    lostBoys=lostBoys[lostBoys$blobIndex!=0,]
    if(nrow(lostBoys)>0){
      euNique=lostBoys$blobIndex
      dfPrevFrame=rbind(matchedLines, lostBoys) 
      bestMuddledMatches=do.call(rbind, lapply(1:length(euNique),  function (k) 
	      {mixoMatchosis(blobIt[(blobIt$blobIndex==euNique[k]),], dfPrevFrame[dfPrevFrame$blobIndex==euNique[k],])}))	    
      matchedLines=matchedLines[(!matchedLines$blobIndex %in% euNique),]
      matchedLines=rbind(matchedLines, bestMuddledMatches)
      matchedLines=matchedLines[order(matchedLines$blobIndex),]
    }
  }
  return(matchedLines)
}

# called from wendyLostBoys()
# when combined fly area has significantly exceeded blob area for many frames, we check to see if eliminating one of those flies 
# markedly improves fit.
clearlySomeoneNeedsToLeave=function(matchedLines, vDeviantBoys){
  options(warn=2)
  #I'll only evaluate one leaver per loop. Arbitrary, but it'll be quicker, and I don't want tons of flies pissing off every frame
  deviantBlobList=unique(matchedLines[vDeviantBoys,]$blobIndex)
  deviantBlobIndex=matchedLines$blobIndex %in% deviantBlobList  #as far as I know, these two lists should always be identical, because I'm trying to
								    #always blank-out framesDeviance when blobs recombine in funny ways
  changeFit=matchedLines$areaDeviance-matchedLines$flyArea
  fitImprovement=abs(changeFit)
  bestImprovement=min(fitImprovement[deviantBlobIndex])
  weakestLink=fitImprovement==bestImprovement
  if(sum(weakestLink)>1){ #this really is only going to ever happen when I impose a floor on blobSize
    weakestIndex=c(1:length(weakestLink))
    firstWeak=weakestIndex[weakestLink][1]
    weakestLink=weakestIndex==firstWeak
  }
  changeImprovement=changeFit[weakestLink][1]
  weakestLinkArea=matchedLines[weakestLink,]$flyArea[1]
  weakestBlobArea=matchedLines[weakestLink,]$blobArea[1]
  weakestBlobIndex=matchedLines$blobIndex==matchedLines[weakestLink,]$blobIndex
  if((((-changeImprovement)/weakestLinkArea)<0.1) & ((weakestLinkArea/weakestBlobArea)>0.2) ){ #there is surely blob measurement error, and we don't want to spuriously throw away small flies in big groups
      
      matchedLines[weakestBlobIndex,"nFlies"]=matchedLines[weakestBlobIndex,"nFlies"]-1
      matchedLines[weakestBlobIndex,"areaDeviance"]=bestImprovement
      matchedLines[weakestBlobIndex,"framesDeviance"]=0
      matchedLines[weakestLink,c("blobIndex","nFlies","flyFrame","blobX","blobY", "blobArea")]=NA
      
    }else{
       matchedLines[weakestBlobIndex,"framesDeviance"]=matchedLines[weakestBlobIndex,"framesDeviance"]/2
    }
  
  return(matchedLines)
}


# called from wendyLostBoys()
# are there any pairs of blobs with poor fly-to-area fits that can be improved by fly-donation?
dealWithSwaps=function(matchedLines, vDeviantBoys){
  #matchedLines[vDeviantBoys,]$framesDeviance= matchedLines[vDeviantBoys,]$framesDeviance+2
  recippityBlobs=matchedLines[c(matchedLines[,"areaDeviance"]<0 & vDeviantBoys),, drop=F]
  recippityBlobs=recippityBlobs[!duplicated(recippityBlobs$blobIndex),, drop=F] 
  donatyFlies=matchedLines[c(matchedLines[,"areaDeviance"]>0 & vDeviantBoys),, drop=F]
  donatyBlobs=donatyFlies[!duplicated(donatyFlies$blobIndex),, drop=F] 
  if(nrow(donatyBlobs)>0 & nrow(recippityBlobs)>0){
    reductionClinic=1
    flyMover=c()

    while(reductionClinic==1){
      weightMetric=sum(abs(donatyBlobs$areaDeviance))+sum(abs(recippityBlobs$areaDeviance))
      recippityMatrix=matrix(rep(recippityBlobs$areaDeviance, nrow(donatyFlies)), nrow=nrow(donatyFlies), ncol=nrow(recippityBlobs), byrow=T)
      donatyMatrix=matrix(rep(donatyBlobs$areaDeviance, nrow(donatyFlies)), nrow=nrow(donatyFlies), ncol=nrow(donatyBlobs), byrow=T)
      donatyBelonging=matrix(rep(donatyBlobs$blobIndex, nrow(donatyFlies)), nrow=nrow(donatyFlies), ncol=nrow(donatyBlobs), byrow=T)
      donatyBelonging=donatyBelonging==donatyFlies$blobIndex
           
      recippityImprovement=recippityMatrix+donatyFlies$flyArea
      donatyImprovement=donatyMatrix-donatyFlies$flyArea
      recippitySumprovement=do.call(cbind, lapply(1:ncol(recippityMatrix), function(k) {(abs(recippityImprovement[,k, drop=F]))+(rowSums(abs(recippityMatrix[,-k, drop=F])))}))
      donatyMatrix[donatyBelonging]=donatyImprovement[donatyBelonging]
      donatyMatrix=abs(donatyMatrix)
      
      donatyVec=rowSums(donatyMatrix)
      
      weightWatcher=recippitySumprovement+donatyVec
      squareWatcher=recippitySumprovement^2+donatyVec^2 #the correct pick should minimise squared error. However, 
							  # for the purposes of deciding whether the fit is good enough to worry about, I guess the absolute improvement in fit is the right
							  #metric????
      biggestLoser=min(squareWatcher)
      weightLoser=weightWatcher[squareWatcher==biggestLoser][1]      
      biggDiff=donatyFlies[rowSums(squareWatcher==biggestLoser)>0,,drop=F]$flyArea[1] #if there were a perfect swap, biggDif would be 1/2 weightLoser
 
      if(biggDiff*1.5 < (weightMetric-weightLoser)){
	donatyBlob=donatyFlies[rowSums(squareWatcher==biggestLoser)>0,,drop=F]$blobIndex[1]
	recippityBlob=recippityBlobs[colSums(squareWatcher==biggestLoser)>0,,drop=F]$blobIndex[1]
	leanFly=c(donatyFlies[rowSums(squareWatcher==biggestLoser)>0,,drop=F]$flyID[1], recippityBlob)
	flyMover=rbind(flyMover, leanFly)
	recippityBlobs[recippityBlobs$blobIndex==as.numeric(leanFly[2]),]$areaDeviance=recippityBlobs[recippityBlobs$blobIndex==as.numeric(leanFly[2]),,drop=F]$areaDeviance+donatyFlies[donatyFlies$flyID==leanFly[1],]$flyArea
	
	donatyFlies[donatyFlies$blobIndex==donatyBlob,]$areaDeviance=donatyFlies[donatyFlies$blobIndex==donatyBlob,,drop=F]$areaDeviance-donatyFlies[donatyFlies$flyID==leanFly[1],]$flyArea
	donatyFlies=donatyFlies[!donatyFlies$flyID %in% leanFly[1],,drop=F]
	donatyBlobs=donatyFlies[!duplicated(donatyFlies$blobIndex),, drop=F] 
	
	matchedLines=rearrangeFlies(matchedLines, leanFly, donatyBlob)
	
	if(nrow(donatyBlobs)<1){
	  reductionClinic=0
	}
      }else{
	reductionClinic=0
      }
    } 
  
  }

  return(matchedLines)
}
  
#called from dealWithSwaps()  
rearrangeFlies=function(matchedLines, leanFly, donatyBlob){

      matchedLines[matchedLines$blobIndex==as.numeric(leanFly[2]),]$areaDeviance=matchedLines[matchedLines$blobIndex==as.numeric(leanFly[2]),,drop=F]$areaDeviance+matchedLines[matchedLines$flyID==leanFly[1],]$flyArea
      matchedLines[matchedLines$blobIndex==as.numeric(leanFly[2]),]$nFlies=matchedLines[matchedLines$blobIndex==as.numeric(leanFly[2]),]$nFlies+1
      passInfo=matchedLines[matchedLines$blobIndex==as.numeric(leanFly[2]),]
      matchedLines[matchedLines$blobIndex==as.numeric(leanFly[2]),]$framesDeviance=0
      matchedLines[matchedLines$blobIndex==donatyBlob,]$areaDeviance=matchedLines[matchedLines$blobIndex==donatyBlob,,drop=F]$areaDeviance-matchedLines[matchedLines$flyID==leanFly[1],]$flyArea
      matchedLines[matchedLines$blobIndex==donatyBlob,]$nFlies=matchedLines[matchedLines$blobIndex==donatyBlob,]$nFlies-1
      matchedLines[matchedLines$blobIndex==donatyBlob,]$framesDeviance=0
      
      matchedLines[matchedLines$flyID==leanFly[1],c("blobIndex", "nFlies", "blobColour", "blobX", "blobY", "blobArea", "blobAngle", "blobA", "blobB","areaDeviance")]=
	passInfo[1,c("blobIndex", "nFlies", "blobColour", "blobX", "blobY", "blobArea", "blobAngle", "blobA", "blobB", "areaDeviance")]
	
      return(matchedLines)
}


