#'@export
extras.getXscale <- function(threshold, max, spacing){
  if ((max-(threshold-1)) %% spacing != 0) warning('Spacing does not evenly divide range')

  pos = seq(threshold-1, max, spacing)
  scale = 10*2**pos
  scale[[1]] = paste0('<', scale[[2]])
  #scale = c(paste0('<', threshold), scale)
  list(pos = pos, scale = scale)
}

#'@export
extras.t_col <- function(v, percent)
  {sapply(v,function(color, name = NULL) {
    rgb.val <- col2rgb(color)
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)

    return(t.col)})
  }
