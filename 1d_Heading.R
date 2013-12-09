mywd = "/Users/ozone96/Desktop/Internship/Guppy data"
myfile = "Guppy Heading Data 7.csv"
setwd(mywd)
data = read.csv(myfile)

x11s = data$Guppy.1.u1
y11s = data$Guppy.1.v1
x12s = data$Guppy.1.u2
y12s = data$Guppy.1.v2
x21s = data$Guppy.2.u1
y21s = data$Guppy.2.v1
x22s = data$Guppy.2.u2
y22s = data$Guppy.2.v2

Axs = data$A.x.coord
Ays = data$A.y.coord
Bxs = data$B.x.coord
Bys = data$B.y.coord
Cxs = data$C.x.coord
Cys = data$C.y.coord
Dxs = data$D.x.coord
Dys = data$D.y.coord

radii1 = vector(mode = "numeric", length = 0)
radii2 = vector(mode = "numeric", length = 0)
headings1 = vector(mode = "numeric", length = 0)
headings2 = vector(mode = "numeric", length = 0)
deviations1 = vector(mode = "numeric", length = 0)
deviations2 = vector(mode = "numeric", length = 0)
headingDiffs = vector(mode = "numeric", length = 0)
radiusDiffs = vector(mode = "numeric", length = 0)

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

getEquation = function(x1, x2, y1, y2){
  m = (y1-y2)/(x1-x2)
  b = y1 - m*x1
  return(c(m,b))
}

centerxs = vector(mode = "numeric", length = 0)
centerys = vector(mode = "numeric", length = 0)

for(i in 1:length(x11s)){
  diagonal1 = getEquation(Axs[i], Ays[i], Cxs[i], Cys[i])
  diagonal2 = getEquation(Bxs[i], Bys[i], Dxs[i], Dys[i])
  centerxs = append(centerxs, (diagonal2[2]-diagonal1[2])/(diagonal1[1]-diagonal2[1]))
  centerys = append(centerys, diagonal1[1]*centerxs+diagonal1[2])
  radii1 = append(radii1, correct(centerxs[i], centerys[i], x11s[i], y11s[i]))
  radii2 = append(radii2, correct(centerxs[i], centerys[i], x21s[i], y21s[i]))
  headings1 = append(headings1, correct(x12s[i], y12s[i], x11s[i], y11s[i]))
  headings2 = append(headings2, correct(x22s[i], y22s[i], x21s[i], y21s[i]))
  deviations1 = append(deviations1, abs(headings1[i] - radii1[i]))
  deviations2 = append(deviations2, abs(headings2[i] - radii2[i]))
}
for(i in 1:length(headings1)){
  headingDiffs = append(headingDiffs, abs(headings1[i]-headings2[i]))
  radiusDiffAngle = abs(radii1[i]-radii2[i])
  if(radiusDiffAngle > pi){
    radiusDiffAngle = (2*pi)-radiusDiffAngle
  }
  radiusDiffs = append(radiusDiffs, radiusDiffAngle)
}
output = data.frame(radii1, radii2, headings1, headings2, deviations1, deviations2, headingDiffs, radiusDiffs)
write.csv(output, file = "1d_transform.csv");
