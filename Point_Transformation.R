#set myWd to the name of the directory where the data is stored.
#set myRawData to the name of the spreadsheet that contains the untransformed data. 
myWd = "/Users/ozone96/Desktop/Internship/Guppy data/"
myRawData = "Guppy Heading Data 5.csv"

#columnNames = c("Guppy.1.u1", "Guppy.1.v1", "Guppy.1.u2", "Guppy.1.v2", "Guppy.2.u1", "Guppy.2..v1", "Guppy.2..u2", "Guppy.2..v2")

setwd(myWd)
rawData = read.csv(myRawData)

#coordinates denoted (x,y) are in the UNTRANSFORMED coordinate space. 
#coordinates denoted (u,v) are in the TRANSFORMED coordinate space. 
Axs = rawData$A.x.coord
Axs = na.omit(Axs)
Ays = rawData$A.y.coord
Ays = na.omit(Ays)
Bxs = rawData$B.x.coord
Bxs = na.omit(Bxs)
Bys = rawData$B.y.coord
Bys = na.omit(Bys)
Cxs = rawData$C.x.coord
Cxs = na.omit(Cxs)
Cys = rawData$C.y.coord
Cys = na.omit(Cys)
Dxs = rawData$D.x.coord
Dxs = na.omit(Dxs)
Dys = rawData$D.y.coord
Dys = na.omit(Dys)
x11s = rawData$Guppy.1.X1
x11s = na.omit(x11s)
y11s = rawData$Guppy.1.Y1
y11s = na.omit(y11s)
x12s = rawData$Guppy.1.X2
x12s = na.omit(x12s)
y12s = rawData$Guppy.1.Y2
y12s = na.omit(y12s)
x21s = rawData$Guppy.2.X1
x21s = na.omit(x21s)
y21s = rawData$Guppy.2.Y1
y21s = na.omit(y21s)
x22s = rawData$Guppy.2.X2
x22s = na.omit(x22s)
y22s = rawData$Guppy.2.Y2
y22s = na.omit(y22s)

u11s = vector(mode = "numeric", length = 0)
v11s = vector(mode = "numeric", length = 0)
u12s = vector(mode = "numeric", length = 0)
v12s = vector(mode = "numeric", length = 0)
u21s = vector(mode = "numeric", length = 0)
v21s = vector(mode = "numeric", length = 0)
u22s = vector(mode = "numeric", length = 0)
v22s = vector(mode = "numeric", length = 0)

#These coordinates will define the boundaries of the transformed coordinate system
Au = 0
Av = 0
Bu = 1000
Bv = 0
Cu = 1000
Cv = 1000
Du = 0
Dv = 1000

for(i in 1:length(x11s)){
  #S is square, M is arithmetic mean, P is product, going in order of composition, e.g. SM is the square of the means and MS is the mean of the squares
  uM = mean(c(Au, Bu, Cu, Du))
  vM = mean(c(Av, Bv, Cv, Dv))
  xM = mean(c(Axs[i], Bxs[i], Cxs[i], Dxs[i]))
  yM = mean(c(Ays[i], Bys[i], Cys[i], Dys[i]))
  xSM = xM^2
  xMS = mean(c(Axs[i]^2,Bxs[i]^2,Cxs[i]^2,Dxs[i]^2))
  ySM = yM^2
  yMS = mean(c(Ays[i]^2,Bys[i]^2,Cys[i]^2,Dys[i]^2))
  uxMP = mean(c(Axs[i]*Au, Bxs[i]*Bu, Cxs[i]*Cu, Dxs[i]*Du))
  uxPM = xM*uM
  uyMP = mean(c(Ays[i]*Au, Bys[i]*Bu, Cys[i]*Cu, Dys[i]*Du))
  uyPM = yM*uM
  vxMP = mean(c(Axs[i]*Av, Bxs[i]*Bv, Cxs[i]*Cv, Dxs[i]*Dv))
  vxPM = xM*vM
  vyMP = mean(c(Ays[i]*Av, Bys[i]*Bv, Cys[i]*Cv, Dys[i]*Dv))
  vyPM = yM*vM
  D = xMS + yMS - xSM - ySM
  a = (uxMP - uxPM + vyMP - vyPM)/D
  b = (uyPM - vxPM - uyMP + vxMP)/D
  c = uM - a*xM + b*yM
  d = vM - a*yM - b*xM 
  u11 = a*x11s[i] - b*y11s[i] + c
  v11 = a*y11s[i] + b*x11s[i] + d
  u12 = a*x12s[i] - b*y12s[i] + c
  v12 = a*y12s[i] + b*x12s[i] + d
  u21 = a*x21s[i] - b*y21s[i] + c
  v21 = a*y21s[i] + b*x21s[i] + d
  u22 = a*x22s[i] - b*y22s[i] + c
  v22 = a*y22s[i] + b*x22s[i] + d
  
  u11s = append(u11s, u11)
  v11s = append(v11s, v11)
  u12s = append(u12s, u12)
  v12s = append(v12s, v12)
  u21s = append(u21s, u21)
  v21s = append(v21s, v21)
  u22s = append(u22s, u22)
  v22s = append(v22s, v22)
}
transformedData = data.frame(u11s, v11s, u12s, v12s, u21s, v21s, u22s, v22s)
write.csv(transformedData, file = "Transformed Guppy Data.csv") #col.names = columnNames)
