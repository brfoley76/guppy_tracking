
# This function calls the subfunction rearrangeData()
# to initialise the raw videoGrabber data
# then calls out to subfunctions in "preprocessFunctions.R"
# preProcessFrame(),identify_duplicates(), and createNewFly()
# to initialise the first frame.
# It then calls out to the main loop, mainLoopBlobProcessing()
# in "videoProcessMain.r" 


#library(plyr)
#library(aspace) # change degree angle to radians: as_radians()

preprocessVideoData=function(path, infile, outfile, centroidX, centroidY, patchRadius){
  setwd(path)
  if(!exists("lp.assign", mode="function")){
    library(lpSolve) # for lp.assign(), hungarian algorithm
    #library(fields)
    source("/media/Data/Documents/VideoAnalysis/check_assigment_algorithm/mainLoopBlobProcessing.r")
    source("/media/Data/Documents/VideoAnalysis/check_assigment_algorithm/preprocessFunctions.r")
    source("/media/Data/Documents/VideoAnalysis/check_assigment_algorithm/matrixMatching.r")
    source("/media/Data/Documents/VideoAnalysis/check_assigment_algorithm/matrixMatchingSubfunctions.r")
    source("/media/Data/Documents/VideoAnalysis/check_assigment_algorithm/specialCasesTracking.r")
    source("/media/Data/Documents/VideoAnalysis/check_assigment_algorithm/complexMatchingSubfunctions.r")

  }


  setwd(path)
  options(width=250)

  myFile=infile
  n_col=max(count.fields(file=myFile, sep = ","))
  readvecEnd=c(n_col)
  readvecStart=c(1)

  myData = read.table(myFile, stringsAsFactors=F, sep=",", fill = TRUE, strip.white=T, header = FALSE, col.names = paste("V", seq_len(n_col), sep = ","))

  fillInnaBlanks=myData[,]==""
  myData[fillInnaBlanks]=NA
  fillInnaBlanks=myData[,]=="NA"
  myData[fillInnaBlanks]=NA
  nextLine=myData[1,]


  #important variables and thresholds
  lineLength=length(nextLine)
  dataLength=7
  #for identifySingletons()
  min_area = 30
  max_area = 900
  min_ratio = 0.3
  max_ratio = 0.45
  maxTrax=(lineLength-2)/dataLength
  #AreaBeta=0.001 # multiplier for guessing the number of flies in a blob. Just to initialise the algorithm

  #create generic vectors I'll need for later functions, 
  #so I don't recreate them every time I call the function
  blobJectNames=c("frame", "blobIndex", "nFlies", "flyFrame", "blobColour", "blobX", "blobY","blobArea", "blobAngle", "blobA", "blobB", "flyID", "deltaX", "deltaY", "flySpeed", "flyArea", "flyAspect", "areaDeviance", "framesDeviance")
  blankBlobJect=data.frame(integer(), integer(), integer(), integer(), I(character()), integer(), integer(), numeric(), integer(), numeric(), numeric(), I(character()), numeric(), numeric(), numeric(), numeric(), numeric(), numeric(), numeric())
  names(blankBlobJect)=blobJectNames
  filename=paste(outfile, "out.csv", sep="")
  write.table(blankBlobJect, file=filename, sep=",", row.names=F, quote=F)
  
  blankBlobJect[1,]=NA
  blankBlobJect$framesDeviance=0.0
  #the variables flySpeed, flyArea, and flyAspect will be running averages, to be used to calculate .
  #Possibly deltaX and deltaY will be too, in order to smooth some frame-to-frame noise.

  #puts the data into a format we can easily use
  #preprocessFunctons.r
  myData=rearrangeData(myData, maxTrax, myRangeMin, dataLength, centroidX, centroidY, patchRadius)

  backupData=myData #I'll make a backup of my data, because, to avoid having to search and extract from the whole fucking thing
  #each line, I'll just lop off the first maxTrax rows each time.
  #myData=backupData
  #myData=myData[myData$frame>30,]
  #myData=prunedData
  minArea=300
  setMinArea=0

  dfNextFrame=myData[1:maxTrax,]
  myData=myData[(maxTrax+1):length(myData[,1]),]
  dfNextFrame=na.omit(dfNextFrame)

  if(sum(!is.na(dfNextFrame$blobX))!=0){
    dfNextFrame=preProcessFrame(dfNextFrame, blankBlobJect, maxTrax, dataLength, min_area, max_area, min_ratio, max_ratio)
    nBlobs=sum(!is.na(dfNextFrame$blobIndex))
    if(nBlobs>1){
      dfNextFrame=identify_duplicates(dfNextFrame)
      }
    if(nBlobs>0){
      dfNextFrame$blobIndex=c(1:nrow(dfNextFrame)) 
      dfNextFrame$nFlies=1
      dfNextFrame$flyID="new"
      dfNextFrame=createNewFly(dfNextFrame)
    }else{
      dfNextFrame=blankBlobJect
      dfNextFrame$frame=1
      dfNextFrame$nFlies=0
    }
  }else{
    dfNextFrame=blankBlobJect
    nBlobs=0
    dfNextFrame$frame=1
    dfNextFrame$nFlies=0
  }
  
  dfPrevFrame=dfNextFrame
  outputData=dfPrevFrame

  mainLoopBlobProcessing(myData, outputData, dfPrevFrame, blankBlobJect, nBlobs, maxTrax, dataLength, min_area, max_area, min_ratio, max_ratio, filename, setMinArea, minArea)
   
}

