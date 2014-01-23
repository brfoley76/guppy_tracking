#detect frames with three or more guppies and create a vector of boolean values marking which data points are on
#frames with too many guppies
trio = function (data){
  result = vector(mode = "logic", length = length(data$frame))
  for(i in 1:length(data$frame)){
    if(data$frame[i] == data$frame[i+1] & data$frame[i] == data$frame[i+2] & i <= length(data$frame)-2){
      result[i] = TRUE
      result[i+1] = TRUE
      result[i+2] = TRUE
    }
  }
  return(result)
}

#detect blobs with an aspect ratio that differs dramatically from the usual guppy aspect ratio and create a vector
#of boolean values as in "trio"
aspect = function(data, minAspect, maxAspect){
  result = vector(mode = "logic", length = length(data$frame))
  for(i in 1:length(data$frame)){
    if(data$flyAspect[i] < minAspect | data$flyAspect[i] > maxAspect){
      result[i] = TRUE
    }
  }
  return(result)
}

#detect blobs whose area differs dramatically from the present value of the running average of their area and 
#create a vector of boolean values as in "trio"
areaDiff = function(data, maxDeviance){
  result = vector(mode = "logic", length = length(data$frame))
  for(i in 1:length(data$frame)){
    if(abs(data$areaDeviance[i]) > maxDeviance){
      result[i] = TRUE
    }
  }
  return(result)
}

#detect blobs whose area differs dramatically from the usual range of areas that guppies fall into
#create a vector of boolean values as in "trio"
area = function(data, minArea, maxArea){
  result = vector(mode = "logic", length = length(data$frame))
  for(i in 1:length(data$frame)){
    if(data$blobArea[i] < minArea | data$blobArea[i] > maxArea){
      result[i] = TRUE
   }
  }
  return(result)
}

#detect blobs whose ellipse area to blob area ratio differs dramatically from the usual ratio found in guppies
#create a vector of boolean values as in "trio"
ellipseRatio = function(data, minRatio, maxRatio){
  ellipseArea = pi*data$blobA*data$blobB
  result = vector(mode = "logic", length = length(data$frame))
  for(i in 1:length(data$frame)){
    ratio = ellipseArea[i]/data$blobArea[i]
    if(ratio > maxRatio | ratio < minRatio){
      result[i] = TRUE
    }
  }
  return(result)
}

#detect blobs that exceed a maximum speed
#create a vector of boolean values as in "trio"
speed = function(data, maxSpeed){
  result = vector(mode = "logic", length = length(data$frame))
  for(i in 1:length(data$frame)){
    if(!is.na(data$flySpeed[i])){
      if(data$flySpeed[i] > maxSpeed){
        result[i] = TRUE
      }
    }
  }
  return(result)
}

#detect frames in which only one blob is identified
#create a vector of boolean values as in "trio"
solo = function(data){
  result = vector(mode = "logic", length = length(data$frame))
  for(i in 1:length(data$frame)){
    if(i > 1){ 
      if(data$frame[i] != data$frame[i-1] & data$frame[i] != data$frame[i+1] & i <= length(data$frame)-1){
        result[i] = TRUE
      }
    }
  }
  return(result)
}

#collate the results of several of the tests above into a vector of strings that can then be added to
#a spreadsheet, with a number representing a flag raised by a particular function. (e.g. trios might be flagged
#with a 1, blobs moving too fast with a 2, blobs with unusual aspect ratios with a 3, etc.)
merge = function(data, ...){
  arguments = matrix(c(...), nrow = length(data$frame), ncol = length(list(...)), byrow = FALSE, dimnames = NULL)
  for(i in 1:length(data$frame)){
    for(j in 1:length(list(...))){
      arguments[i,j] = as.character(j)
    }
  }
  result = vector(mode = "character", length = 0)
  for(k in 1:length(data$frame)){
    string = ""
    for(l in 1:length(list(...))){
      string = paste(string, arguments[k,l], sep = " ")
    }
    result = append(result, string)
  }
  return(result)
}
