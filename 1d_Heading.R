get.angle = function(x2, y2, x1, y1){
  angle = atan((y1-y2)/(x1-x2))
  if(angle > 0 & (x1-x2)<0){
    angle = pi - angle
  }
  if(angle < 0){
    angle = 2*pi + angle
    if((x1-x2) < 0){
      angle = 3*pi-angle
    }
  }
  angle = round(angle*(180/pi))
  return(angle)
}

correct = function(angle, deltaX){
  if(!is.na(deltaX)){
    if(angle > 180 & deltaX > 0){
      newAngle = angle - 180
    }else if(angle < 180 & deltaX > 0){
      newAngle = angle + 180
    }else{
      newAngle = angle
    }
  }else{
    newAngle = angle
  }
  return(newAngle)
}

equation = function(x1, x2, y1, y2){
  m = (y1-y2)/(x1-x2)
  b = y1 - m*x1
  return(c(m,b))
}

oneDimension = function(data, centerX, centerY){
  goodRows = which(is.na(data$Flag))
  goodData = data[goodRows,]
  x = goodData$blobX
  y = goodData$blobY
  dx = goodData$deltaX
  dy = goodData$deltaY
  angle = goodData$blobAngle
  goodData$radius = mapply(get.angle, x, y, MoreArgs = list(centerX, centerY))
  guppy2Rows = which(1:length(goodData$frame) %% 2 == 0)
  guppy2Data = goodData[guppy2Rows,]
  guppy1Rows = guppy2Rows - 1
  guppy1Data = goodData[guppy1Rows,]
  guppy1Parsimonious = vector(mode = "numeric", length = 0)
  guppy2Parsimonious = vector(mode = "numeric", length = 0)
  angle1 = guppy1Data$blobAngle
  angle2 = guppy2Data$blobAngle
  deltaX1 = guppy1Data$deltaX
  deltaX2 = guppy2Data$deltaX
  for(i in 1:length(guppy1Rows)){
    corrected1 = correct(angle1[i], deltaX1[i])
    corrected2 = correct(angle2[i], deltaX2[i])
    if(i > 1){
      if(abs(corrected1 - guppy1Parsimonious[i-1]) < abs(angle1[i] - guppy1Parsimonious[i-1])){
        output1 = corrected1
      }else{
        output1 = angle1[i]
      }
      if(abs(corrected2 - guppy2Parsimonious[i-1]) < abs(angle2[i] - guppy2Parsimonious[i-1])){
        output2 = corrected2
      }else{
        output2 = angle2[i]
      }
    }else{
     output1 = angle1[i]
     output2 = angle2[i]
    }
    guppy1Parsimonious = append(guppy1Parsimonious, output1)
    guppy2Parsimonious = append(guppy2Parsimonious, output2)
  }
  guppy1Deviance = guppy1Parsimonious - guppy1Data$radius
  guppy2Deviance = guppy2Parsimonious - guppy2Data$radius
  result = data.frame(guppy1Parsimonious, guppy2Parsimonious, guppy1Deviance, guppy2Deviance, guppy1Deviance-guppy2Deviance)
  colnames(result) = c("Guppy 1 Angle", "Guppy 2 Angle", "Guppy 1 Deviance", "Guppy 2 Deviance", "Difference")
  return(result)
}
