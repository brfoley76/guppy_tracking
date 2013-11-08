# These subfunctions are called from both 
# preprocessVideoData() and mainLoopBlobProcessing()


#organises the dataset such that every blob has a line, indexed by frame and blob index
rearrangeData=function(myData, maxTrax, myRangeMin, dataLength, centroidX, centroidY, patchRadius){
  mySubNames=c("frame","blobIndex","blobColour","blobX","blobY","blobArea","blobAngle","blobA", "blobB")
  #I need to restructure the data so each blob is its own row. It will be more efficient to do this all at once, than for each line one at a time I suspect
  subData=myData
  myData=c()
  for(i in 1:maxTrax){
    myRangeMin=3+(i-1)*dataLength
    myRangeMax=myRangeMin+(dataLength-1)
    mySub=subData[,c(1:2, myRangeMin:myRangeMax)]
    names(mySub)=mySubNames
    mySub[,2]=i
    myData=rbind(myData, mySub)
  }
  myData=myData[order(myData[,1]),]
  myData=cliquey(myData,centroidX, centroidY, patchRadius)
  return(myData)
}


# called from rearrangeData(). removes flies that are off patch. 
cliquey=function(myData, centroidX, centroidY, patchRadius){
  distanceFromCentre=(((myData$blobX-centroidX)^2)+(myData$blobY-centroidY)^2)^0.5
  inOut=distanceFromCentre/(patchRadius)
  exclusionAry=inOut>1
  exclusionAry[is.na(exclusionAry)]=FALSE
  myData[exclusionAry,c("blobIndex", "blobX","blobY", "blobArea")]=NA
  
  return(myData)
}

#calls subfunctions remove_fuckups(), identify_singletons() for each fly individually, and then identify_duplicates() on t.
preProcessFrame=function(dfNextFrame, blankBlobJect, maxTrax, dataLength, min_area, max_area, min_ratio, max_ratio){
  blobDataFrame=c()
  for(i in 1:maxTrax){
    vFly=dfNextFrame[i,]
    if(!is.na(vFly[3])){
      vFly=remove_fuckups(vFly, dataLength)
      if(!is.na(vFly[3])){
	vFly=identify_singletons(vFly, min_area, max_area, min_ratio, max_ratio)  
	newLine=blankBlobJect
	newLine[names(vFly)]=vFly
	blobDataFrame=rbind(blobDataFrame, newLine)
      }
    }
  }
  return(blobDataFrame)
}


 

# called from preProcessFrame()
# there are occasional clear screw-ups in the video processing. This eliminates those.
remove_fuckups=function(vFly, dataLength){
    if(sum(is.na(vFly))==0){
      #get rid of cases where area and angle are ~0 aka fuckups
      if(sum(as.numeric(vFly[c("blobX", "blobY", "blobArea")]))<0.1){
	vFly=rep(NA, (dataLength+2))
	}else{
	  if(as.numeric(vFly["blobArea"])<80){
	    vFly=rep(NA, (dataLength+2))
	  }
      }
    }
  return(vFly)
}

# called from preProcessFrame()
#identify all pretty-unambiguous single-flies by area and aspect ratio. 
#We so far only use this information to flag when a presumed singleton fly seems to have too many flies
#in the function matchPair()
identify_singletons=function(vFly, min_area, max_area, min_ratio, max_ratio){

    flyRatio=vFly$blobB/vFly$blobA
    flyTrueEllipse=vFly$blobArea/(vFly$blobA * vFly$blobB * 3.14)
    #only if all 4 conditions are met
    if(flyTrueEllipse>0.95 && ((vFly$blobArea > min_area) && (vFly$blobArea < max_area)) && ((flyRatio > min_ratio) && (flyRatio < max_ratio))){
      nFlies=1
    } else {
      nFlies=0
    }
  vFly=data.frame(vFly, nFlies)
  return(vFly)
}

# called from preProcessFrame()
#check to see if any of the (non-singleton?) blobs are identical
identify_duplicates=function(dfNextFrame){
  sorTable=cbind(dfNextFrame$blobX, dfNextFrame$blobY, dfNextFrame$blobArea, dfNextFrame$blobIndex)
  sorTableOrder=order(dfNextFrame$blobX, dfNextFrame$blobY, dfNextFrame$blobArea, na.last = NA)
  sorTable=sorTable[sorTableOrder,]
  subTable=rbind(c(0,0,0,0), sorTable[1:(length(sorTable[,1])-1),])
  repVec=c(abs(sorTable[,1]-subTable[,1])+abs(sorTable[,2]-subTable[,2])+abs(sorTable[,3]-subTable[,3]))
  repVec=c(repVec<10) #might need to play around with this threshold. Also, might need to divide the area by something, to avoid giving it undue weight?
  dupVec=sorTable[repVec,4]
  
  if(length(dupVec>0)){
    dfNextFrame=dfNextFrame[!(dfNextFrame$blobIndex %in% dupVec),]
  }
  dfNextFrame$blobIndex=c(1:nrow(dfNextFrame))
  return(dfNextFrame)
}

# called from preprocessVideoData(), mainLoopBlobProcessing(), matrixMatching(), and complexMatch(). *this seems excessive*
# any time existing fly indices can't be assigned to a blob, we'll create a new fly index and assume a fly has joined a patch 
createNewFly=function(dfNextFrame){
  orphanedMisfits=c(dfNextFrame$flyID)=="new"
  orphanedMisfits[is.na(orphanedMisfits)]=FALSE
  n0s=6-nchar(dfNextFrame[1,"frame"])
  n0vec=rep(0, n0s)
  newFlyNameVec=paste("fly", paste(n0vec, collapse=""), dfNextFrame$frame[orphanedMisfits], "_", dfNextFrame$blobIndex[orphanedMisfits], sep="")
  dfNextFrame$flyID[orphanedMisfits]=newFlyNameVec
  dfNextFrame$flyFrame[orphanedMisfits]=1
  dfNextFrame$flyArea[orphanedMisfits]=dfNextFrame$blobArea[orphanedMisfits]
  dfNextFrame$flyAspect[orphanedMisfits]=dfNextFrame$blobB[orphanedMisfits]/dfNextFrame$blobA[orphanedMisfits]
  dfNextFrame$areaDeviance[orphanedMisfits]=0
  dfNextFrame$nFlies[orphanedMisfits]=1
  return(dfNextFrame)
}
