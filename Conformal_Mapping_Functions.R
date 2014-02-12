anamorph.points = function(n, data, headings, areaL, areaW){
  anchors = vector(mode = "character", length = 0)
  fish = vector(mode = "character", length = 0)
  anchorPrompts = c("x-coordinates of the lower left anchor points:", "y-coordinates of the lower 
                    left anchor points", "x-coordinates of the lower right anchor points:",
                    "y-coordinates of the lower right anchor points:", "x-coordinates of the 
                    upper right anchor points:", "y-coordinates of the upper right anchor points:",
                    "x-coordinates of the upper left anchor points:", "y-coordinates of the upper
                    left anchor points")
  if(is.null(headings)){
    for(i in 1:8){
      anchors[i] = readline(paste("Please enter the name of the dataframe column containing the ", anchorPrompts[i], sep = "", collapse = NULL))
      while(anchors[i] %in% colnames(data) == FALSE){
       anchors[i] = readline(paste("The name you entered did not match any data fram columns. ","Please enter the name of the dataframe column containing the ", anchorPrompts[i], sep = "", collapse = NULL))
      }
      head(data[anchors[i]])
    }
    for(j in 1:n){
      fish[(2*j)-1] = readline(paste("Please enter the name of the dataframe column containing the x-coordinates of the first
                                     of point set", as.character(j), ":", sep = "", collapse = NULL))
      while(fish[(2*j)-1] %in% colnames(data) == FALSE){
        fish[(2*j)-1] = readline(paste("The name you entered did not match any data frame columns. ","Please enter the name of the dataframe column containing the x-coordinates of the first
                                     of point set", as.character(j), ":", sep = "", collapse = NULL))
      }
      head(data[fish[(2*j)-1]])
      fish[2*j] = readline(prompt = paste("Please enter the name of the dataframe column containing the y coordinates of point set ", as.character(j), ":", sep = "", collapse = NULL))
      while(fish[2*j] %in% colnames(data) == FALSE){
        fish[2*j] = readline(prompt = paste("The name you entered did not match any data frame columns. ","Please enter the name of the dataframe column containing the y coordinates of point set ", as.character(j), ":", sep = "", collapse = NULL))
      }
      head(data[fish[2*j]])
    }
  }else{
    match = TRUE
    for(i in 1:8+(2*n)){
      if(headings[i] %in% colnames(data) == FALSE){
        match = FALSE
      }
    }
    if(match == TRUE){
      anchors = headings[1:8]
      fish = headings[9:(8+(2*n))]
    }else{
      stop("The vector of headings provided does not match the column headings of the data frame.")
    }
  }
  if(is.null(areaL)){areaL = as.numeric(readline(prompt = "Please enter the desired length of the transformed coordinate area: "))}
  if(is.null(areaW)){areaW = as.numeric(readline(prompt = "Please enter the desired width of the transformed coordinate area: "))}
  lengthGet = nrow(data)
  ordered = matrix(data = NA, nrow = lengthGet, ncol = 8+length(fish), byrow = FALSE, dimnames = NULL)
  out = matrix(data = NA, nrow = lengthGet, ncol = length(fish), byrow = FALSE, dimnames = NULL)
  for(l in 1:8){
    ordered[,l] = data[,anchors[l]]
  }
  for(m in 9:(8+(2*n))){ordered[,m] = data[,fish[m-8]]}
  for(o in 1:lengthGet){
    uM = areaL/2
    vM = areaW/2
    xM = mean(c(ordered[o,1], ordered[o,3], ordered[o,5], ordered[o,7]))
    yM = mean(c(ordered[o,2], ordered[o,4], ordered[o,6], ordered[o,8]))
    xSM = xM^2
    xMS = mean(c(ordered[o,1]^2,ordered[o,3]^2,ordered[o,5]^2,ordered[o,7]^2))
    ySM = yM^2
    yMS = mean(c(ordered[o,2]^2,ordered[o,4]^2,ordered[o,6]^2,ordered[o,8]^2))
    uxMP = mean(c(0, ordered[o,3]*areaL, ordered[o,5]*areaL, 0))
    uxPM = xM*uM
    uyMP = mean(c(0, ordered[o,4]*areaL, ordered[o,6]*areaL, 0))
    uyPM = yM*uM
    vxMP = mean(c(0, 0, ordered[o,5]*areaW, ordered[o,7]*areaW))
    vxPM = xM*vM
    vyMP = mean(c(0, 0, ordered[o,6]*areaW, ordered[o,8]*areaW))
    vyPM = yM*vM
    D = xMS + yMS - xSM - ySM
    a = (uxMP - uxPM + vyMP - vyPM)/D
    b = (uyPM - vxPM - uyMP + vxMP)/D
    c = uM - a*xM + b*yM
    d = vM - a*yM - b*xM 
    for(p in 1:n){
      out[o,(2*p)-1] = a*ordered[o, 7+(2*p)] - b*ordered[o, 8+(2*p)] + c
      out[o, 2*p] = a*ordered[o, 8+(2*p)] + b*ordered[o, 7+(2*p)] + d
    }
  }
  final = data.frame(out)
  colnames(final) = fish
  return(final)
}

anamorph.image = function(image, vertices, areaL, areaW){
  anchors = vector(mode = "numeric", length = 0)
  original = as.array(image)
  anchorPrompts = c("x-coordinate of the lower left anchor points:", "y-coordinate of the lower 
                    left anchor points", "x-coordinate of the lower right anchor points:",
                    "y-coordinate of the lower right anchor points:", "x-coordinate of the 
                    upper right anchor points:", "y-coordinate of the upper right anchor points:",
                    "x-coordinate of the upper left anchor points:", "y-coordinate of the upper
                    left anchor points")
  if(is.null(vertices)){
    for(i in 1:8){
      anchors[i] = as.numeric(readline(paste("Please enter the ", anchorPrompts[i], sep = "", collapse = NULL)))
    }
  }else{
      anchors = vertices[1:8]
  }
  if(is.null(areaL)){areaL = as.numeric(readline(prompt = "Please enter the desired length of the transformed coordinate area: "))}
  if(is.null(areaW)){areaW = as.numeric(readline(prompt = "Please enter the desired width of the transformed coordinate area: "))}
  uM = areaL/2
  vM = areaW/2
  xM = mean(c(anchors[1], anchors[3], anchors[5], anchors[7]))
  yM = mean(c(anchors[2], anchors[4], anchors[6], anchors[8]))
  xSM = xM^2
  xMS = mean(c(anchors[1]^2,anchors[3]^2,anchors[5]^2,anchors[7]^2))
  ySM = yM^2
  yMS = mean(c(anchors[2]^2,anchors[4]^2,anchors[6]^2,anchors[8]^2))
  uxMP = mean(c(0, anchors[3]*areaL, anchors[5]*areaL, 0))
  uxPM = xM*uM
  uyMP = mean(c(0, anchors[4]*areaL, anchors[6]*areaL, 0))
  uyPM = yM*uM
  vxMP = mean(c(0, 0, anchors[5]*areaW, anchors[7]*areaW))
  vxPM = xM*vM
  vyMP = mean(c(0, 0, anchors[6]*areaW, anchors[8]*areaW))
  vyPM = yM*vM
  D = xMS + yMS - xSM - ySM
  a = (uxMP - uxPM + vyMP - vyPM)/D
  b = (uyPM - vxPM - uyMP + vxMP)/D
  c = uM - a*xM + b*yM
  d = vM - a*yM - b*xM 
  transformed = array(data = NA, dim = c(areaW, areaL, 3),dimnames = NULL)
  for(i in round(min(c(vertices[1], vertices[7]))):round(max(c(vertices[3], vertices[5])))){
    for(j in round(min(c(vertices[2],vertices[4]))):round(max(c(vertices[6], vertices[8])))){
      transformed[round(a*j + b*i + d), round(a*i - b*j + c),] = original[j,i,]
    }
  }
  for(k in 1:areaL){
    for(l in 1:areaW){
      if(is.na(transformed[k,l,])){
        transformed[k,l,] = mean(c(transformed[k+1,l,], transformed[k-1,l,], transformed[k,l+1,], transformed[k,l-1,]), na.rm = TRUE)
      }
    }
  }
  return(as.Image(transformed))
}
