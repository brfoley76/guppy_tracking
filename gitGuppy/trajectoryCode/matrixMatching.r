# This function first calculates what blobs are close enough within and between frames to potentially
# consider "neighbors" based on blob size and distance between blobs
# it then attempts increasingly complex methods to assign flies from the previous frame
# to the neighboring blobs in the next frame.
# It is called from function mainLoopBlobProcessing()
# and calls out to subfunctions identifyJoins() and identifyLeaves()
# matchPair(), identifyBest(), and complexMatch() in matrixMatchingSubfunctions.r
# createNewFly() in preprocessFunctions.r -> maybe we can change this? probably it would be best to do all the new-fly creation from within mainLoopBlobProcessing().
# and all the crazy ad hoc shit for shuffling flies and stuff that happens in wendyLostBoys()
# is then called by the functions in specialCasesTracking.r


matrixMatching=function(dfNextFrame, dfPrevFrame){
   
  #make couple matrices. First, a pariwise distance matrix between newvec and oldvec.
  # then - an overlap matrix, based on distances between blobs, and their long-radii
  # then - an old nFly vector 
  # then - identify unambiguous singleton matches (NA out those columns & rows in the matrix)
  # then - identify unambiguous joins (NA out those columns in the matrix)
  # then - identify unabiguous singleton splits in the matrix (NA out those columns, adjust nFLy)
  # then - starting with the lowest index nextBlob, moving to the highest, identify overlapping (maybe serially overlapping) old/new blobs, then allocate flies by area. 
  #indiVecPrev %in% c(1, 2, 8)
  
  #I'll do the matching on a blob-by-blob basis, not blobs-to-flies, for simplicity.
  nPrevFlies=length(dfPrevFrame$blobIndex)
  if(nPrevFlies>1){
    dfPrevBlobIndex=c(TRUE, c(dfPrevFrame$blobIndex[2:nPrevFlies]!=dfPrevFrame$blobIndex[1:(nPrevFlies-1)]))
  }else{
    dfPrevBlobIndex=c(TRUE)
  }

  #print(dfNextFrame)
  #print(dfPrevFrame)
  #this dataframe will hold the next frame matches
  matchedLines=c()
  
  nextXY=matrix(cbind(dfNextFrame$blobIndex, dfNextFrame$blobX, dfNextFrame$blobY), nrow=length(dfNextFrame$blobIndex), ncol=3)
  prevXY=matrix(cbind(dfPrevFrame$blobIndex[dfPrevBlobIndex], dfPrevFrame$blobX[dfPrevBlobIndex], dfPrevFrame$blobY[dfPrevBlobIndex]), nrow=sum(dfPrevBlobIndex), ncol=3)
 
  mDistMatrix=rdist(nextXY[,2:3, drop=F], prevXY[,2:3, drop=F])  #rdist() in package fields
  #the rows are nextLine, the columns are prevLine
  mDistMatrix=mDistMatrix+1 #no division by zero!  
  #the next 2 matrices are used to calculate a maximum radius for determining potential overlap

  nextRadMatrix=matrix(rep(dfNextFrame$blobA, sum(dfPrevBlobIndex)), nrow=length(dfNextFrame$blobA), ncol=sum(dfPrevBlobIndex))
  prevRadMatrix=matrix(rep(dfPrevFrame$blobA[dfPrevBlobIndex], length(dfNextFrame$blobA)), nrow=length(dfNextFrame$blobA), ncol=sum(dfPrevBlobIndex), byrow=T)

  nextAreaMatrix=matrix(rep(dfNextFrame$blobArea, sum(dfPrevBlobIndex)), nrow=length(dfNextFrame$blobA), ncol=sum(dfPrevBlobIndex))
  prevAreaMatrix=matrix(rep(dfPrevFrame$blobArea[dfPrevBlobIndex], length(dfNextFrame$blobA)), nrow=length(dfNextFrame$blobA), ncol=sum(dfPrevBlobIndex), byrow=T)
  

  compAreaMatrix=nextAreaMatrix/prevAreaMatrix

  littleMatrix=prevRadMatrix<=nextRadMatrix
  axisMatrix=prevRadMatrix
  axisMatrix[littleMatrix]=nextRadMatrix[littleMatrix]
  axisMatrix=axisMatrix*2
  
  neighborhoodMatrix=axisMatrix/mDistMatrix
  #try this: create a row and column disambiguity matrix. Anything more than half the max of that row or column is counted as a hit. Use this as the path matrix.
  colMax=apply(neighborhoodMatrix, 2, max)
  colMax[colMax<0.5]=0 # anything more than two major axes away from the nearest neighbor should almost certainly be ignored. Maybe I should only worry about this in the join/leave functions?
  rowMax=apply(neighborhoodMatrix, 1, max)
  rowMax[rowMax<0.5]=0
  neighbColMax=t(t(neighborhoodMatrix) / colMax)
  neighbRowMax=neighborhoodMatrix / rowMax
  neighbColMaxCut=neighbColMax*0
  neighbRowMaxCut=neighbRowMax*0
  neighbColMaxCut[neighbColMax > 0.5]=1
  neighbRowMaxCut[neighbRowMax > 0.5]=1
  
  pathMatrix=neighbColMaxCut+neighbRowMaxCut
  pathMatrix[pathMatrix > 0]=1
  recipMatrix=(neighbColMaxCut+neighbRowMaxCut)/2
  recipMatrix[recipMatrix < 1]=0
  
  removalIndex=dfNextFrame$blobIndex
  

  #prevOrphans=prevIndex
    
  #I don't really want to do the following here, but if I don't, the leavers and joiners fuck up my neighborhood matrices
    unRequited=dfNextFrame$blobIndex[(rowSums(recipMatrix) == 0)]
    if(length(unRequited)>0){
      recipMatrix=identifyJoins(pathMatrix, recipMatrix, compAreaMatrix, unRequited)
      joiners=dfNextFrame$blobIndex[rowSums(recipMatrix)==0]
      if(length(joiners)>0){
	joinLine=dfNextFrame[(dfNextFrame$blobIndex %in% joiners),]
	joinLine$flyID="new"
	matchedLines=rbind(matchedLines, joinLine)
	removalIndex=removalIndex[!(removalIndex %in% joiners)]
      }
    }

    unRequited=dfPrevFrame$blobIndex[dfPrevBlobIndex][(colSums(recipMatrix) == 0)]
    if(length(unRequited)>0){
      recipMatrix=identifyLeaves(pathMatrix, recipMatrix, compAreaMatrix, unRequited)
      #this should be sufficient to ensure the trajectory ends here? We can try and splice trajectories exogenously if there are short gaps
      #leavers=prevIndex[rowSums(recipMatrix)==0]
      #prevIndex=prevIndex[!(prevIndex %in% leavers)]  
    }
 
  pathMatrix=recipMatrix
  
  while(length(removalIndex)>0){
    checkNext=removalIndex[1] #initiate the neighborhood finding
    groupNextIndex=checkNext # the index vector of the current line neighborhood
    groupPrevIndex=c() # the index vector of the previous line in the neighborhood
    k=0
    
    #the following loop is awesome. I iteratively check which old blobs overlap with which new blobs and vice versa, 
    #until all the serially overlapping blobs are accounted for
    while(k==0){
      OldGroupNextIndex=groupNextIndex
      OldGroupPrevIndex=groupPrevIndex
      #explore all possible paths connecting the blob checkNext and other blobs in the matrix
      #use colSums and rowSums
      #and a matrix with 1s and 0s where there are overlaps

      prevJoin=colSums(pathMatrix[c(dfNextFrame$blobIndex %in% groupNextIndex),, drop=F])
      groupPrevIndex=dfPrevFrame$blobIndex[dfPrevBlobIndex][prevJoin>0]
      nextJoin=rowSums(pathMatrix[,c(dfPrevFrame$blobIndex[dfPrevBlobIndex] %in% groupPrevIndex),drop=F])
      groupNextIndex=dfNextFrame$blobIndex[nextJoin>0]
      
      if(setequal(groupNextIndex, OldGroupNextIndex) && setequal(groupPrevIndex, OldGroupPrevIndex)){
	k=1     
      }
    }
    removalIndex=removalIndex[!removalIndex %in% groupNextIndex] 
    #prevOrphans=prevOrphans[!prevOrphans %in% groupPrevIndex]
    
 ##################################################
 ### This is the guts of the blob-> fly matching.
 ######## I go from the simplest to harder cases. 
 ### I start with unambiguous 1:1 matches, 
 #### then on to multiple nearby blobs, still with pretty clear matches (nPrevBlobs == nNextBlobs)
 #### then finally to tough cases.
 ####
    
    if((length(groupNextIndex)==1)&&(length(groupPrevIndex)==1)){
      #easiest case: there is a single match between lines  
	dfNextBlob=dfNextFrame[dfNextFrame$blobIndex==groupNextIndex,]
	dfPrevBlob=dfPrevFrame[c(dfPrevFrame$blobIndex == groupPrevIndex),]
	dfNextBlob=matchPair(dfNextBlob, dfPrevBlob)
	matchedLines=rbind(matchedLines, dfNextBlob)
    }else {

	#this is a grab bag of possible complication
	#Maybe (e.g.) a complex blob got mistakenly flagged as a 1?
	#can hopefully break it down into manageable bites.
	

	bestUnambiguousMatches=identifyBest(groupNextIndex, groupPrevIndex, pathMatrix, compAreaMatrix, mDistMatrix)
	groupNextIndex=groupNextIndex[!(groupNextIndex %in% bestUnambiguousMatches[,1])]
	groupPrevIndex=groupPrevIndex[!(groupPrevIndex %in% bestUnambiguousMatches[,2])]
	if(nrow(bestUnambiguousMatches)>0){
	  #deal with the best 1:1 unambiguous guys first
	  dfNextBlob=do.call(rbind, lapply(1:nrow(bestUnambiguousMatches),  function (k) 
	    {matchPair(dfNextFrame[(dfNextFrame$blobIndex==bestUnambiguousMatches[k,1]),], dfPrevFrame[c(dfPrevFrame$blobIndex==bestUnambiguousMatches[k,2]),])}))
	  matchedLines=rbind(matchedLines, dfNextBlob)  
	}
	if(length(groupNextIndex)>0 & length(groupPrevIndex)>0){
	  bestMuddledMatches=complexMatch(dfNextFrame, dfPrevFrame, groupNextIndex, groupPrevIndex)
	  matchedLines=rbind(matchedLines,bestMuddledMatches) 
	}     
    }
  }
  
   unMatchedLines=dfNextFrame[!dfNextFrame$blobIndex %in% matchedLines$blobIndex,,drop=F]
  if(nrow(unMatchedLines)>0 | (sum(matchedLines$flyID=="new") > 0) ){
    if(nrow(unMatchedLines)>0){
      unMatchedLines$flyID="new"
      matchedLines=rbind(matchedLines, unMatchedLines)
    }
    matchedLines=createNewFly(matchedLines)
    matchedLines=matchedLines[order(matchedLines$blobIndex),]  
  }
  
  
  vDeviantBoys=matchedLines$framesDeviance>1
  lostBoys=dfPrevFrame[!dfPrevFrame$flyID %in% matchedLines$flyID,,drop=F]  #this isn't elegant, but it should recover flies that were dropped (say) in the complex matching
  
  if(sum(vDeviantBoys)| nrow(lostBoys)>0){
    matchedLines=wendyLostBoys(matchedLines, vDeviantBoys, lostBoys)
  }
  
  dfNextFrame=matchedLines[order(matchedLines$blobIndex),]  
  dfNextFrame=dfNextFrame[!is.na(dfNextFrame$blobIndex),]
  return(dfNextFrame)
}



