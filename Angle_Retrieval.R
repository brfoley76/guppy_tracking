myWd = "/Users/ozone96/Desktop/"
myData = "Guppy Heading Data 4 copy.csv"
setwd(myWd)
data = read.csv(myData)
correct = function(x1, x2, y1, y2){
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
  return(angle)
}

x11s = data$Guppy.1.u1
y11s = data$Guppy.1.v1
x12s = data$Guppy.1.u2
y12s = data$Guppy.1.v2
x21s = data$Guppy.2.u1
y21s = data$Guppy.2.v1
x22s = data$Guppy.2.u2
y22s = data$Guppy.2.v2
x11s = x11s[!is.na(x11s)]
y11s = y11s[!is.na(y11s)]
x12s = x12s[!is.na(x12s)]
y12s = y12s[!is.na(y12s)]
x21s = x21s[!is.na(x21s)]
y21s = y21s[!is.na(y21s)]
x22s = x22s[!is.na(x22s)]
y22s = y22s[!is.na(y22s)]

Guppy1Headings = vector(mode = "numeric", length = 0)
Guppy2Headings = vector(mode = "numeric", length = 0)
Guppy1Deviation = vector(mode = "numeric", length = 0)
Guppy2Deviation = vector(mode = "numeric", length = 0) 
HeadingDiffs = vector(mode = "numeric", length = 0)
Distances = vector(mode = "numeric", length =0)

for(i in 1:length(x11s)){
  new1Heading = correct(x12s[i], x11s[i], y11s[i], y12s[i])
  new2Heading = correct(x22s[i], x21s[i], y21s[i], y22s[i])
  Guppy1Headings = append(Guppy1Headings, new1Heading)
  Guppy2Headings = append(Guppy2Headings, new2Heading)
}
for(j in 1:length(Guppy1Headings)){
  newDeviationAngle = correct(x21s[j],x11s[j],y21s[j],y11s[j])
  newDiff1Angle = abs(newDeviationAngle-Guppy1Headings[j])
  if(newDiff1Angle > pi){
    newDiff1Angle = newDiff1Angle-pi
  }
  Guppy1Deviation = append(Guppy1Deviation, newDiff1Angle)
  newDeviationAngle2 = correct(x11s[j], x21s[j], y11s[j], y21s[j])
  newDiff2Angle = abs(newDeviationAngle2-Guppy2Headings[j])
  if(newDiff2Angle > pi){
    newDiff2Angle = newDiff2Angle-pi
  }
  Guppy2Deviation = append(Guppy2Deviation, newDiff2Angle)
  newDiffAngle = abs(Guppy1Headings[j]-Guppy2Headings[j])
  if(newDiffAngle > pi){
    newDiffAngle = newDiffAngle-pi
  }
  HeadingDiffs = append(HeadingDiffs, newDiffAngle)
}
for(k in 1:length(x11s)){
  newDist = sqrt((x11s[k]-x21s[k])^2+(y11s[k]-y21s[k])^2)
  Distances = append(Distances, newDist)
}

results = data.frame(Guppy1Headings, Guppy2Headings, Guppy1Deviation, Guppy2Deviation, HeadingDiffs, x11s, x12s, x21s, x22s, y11s, y12s, y21s, y22s, Distances)
write.csv(results, file="Corrected Angle Data.csv")
