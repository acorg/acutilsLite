#'@export
titerplot.colourFacetLabels <- function(gg, ordered_colours){
  gg <- ggplot_gtable(ggplot_build(gg))
  stripr <- which( grepl('strip-r', gg$layout$name) | grepl('strip-t', gg$layout$name) )
  fills <- ordered_colours
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', gg$grobs[[i]]$grobs[[1]]$childrenOrder))
    gg$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  gg
}



#'@export
titerplot.styleaxes = function(gg, x = 'detect', Yscale, expand_multiplier = 0, gridlines = T){

  if (x=='detect') x = as.character(rlang::quo_get_expr(ggplot_build(gg)$plot$layers[[1]]$mapping$x))

  longtiters = gg$data
  gg = gg+theme_classic() +
    theme(axis.text.x = element_text(size = 11, angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(size = 11),
          panel.border = element_rect(color = "grey50", fill  = NA)) +
    scale_y_continuous(breaks = Yscale[[1]], labels = Yscale[[2]], expand = c(0.02,0.02), limits = range(Yscale[[1]]) + c(-1,1)*expand_multiplier)

if (gridlines){
  gg= gg+
    theme(axis.line=element_line(),
          panel.grid.major.y = element_line(colour="grey80", size=0.25),
          panel.grid.major.x = element_line(colour="grey80", size=0.25),
          strip.text.x = element_text(size = 15))
}
  else{

    gg= gg+
      theme(axis.line=element_line(),
            panel.grid.major.y = element_blank(),
            panel.grid.major.x =element_blank(),
            strip.text.x = element_text(size = 15))

  }

  if (x == 'sr'){
    map <- setNames( longtiters$sr_short, longtiters$sr )
    gg = gg + scale_x_discrete(labels = map[(unique(as.character(longtiters$sr)))]) + xlab('Serum')

  }
  if (x == 'ag'){
    map <- setNames( longtiters$ag_short, longtiters$ag )
    gg = gg + scale_x_discrete(labels = map[(unique(as.character(longtiters$ag)))]) + xlab('Antigen')
  }
  return(gg)

}


#'@export
titerplot.facet <- function(gg, facet = 'detect'){
  if (facet=='detect'){
  x = as.character(rlang::quo_get_expr(ggplot_build(gg)$plot$layers[[1]]$mapping$x))
  facet = c('ag', 'sr')[[which(c('ag', 'sr') != x)]]
  }


  longtiters = gg$data

  map_names <- setNames( unlist(longtiters[,paste0(facet, '_short')]), unlist(longtiters[,facet]) )


  facet_labeller = function(v){
    idx = as.character(unlist(v))
    list(map_names[idx])
  }

  gg + facet_grid(as.formula(paste("~", facet)), labeller = facet_labeller)  # facet by ag
}

#'@export
titerplot.shadeClades <- function(gg, x = 'detect', colors = acutilsLite::h3_clade_colors, transparency = 80){

  if (x=='detect') x = as.character(rlang::quo_get_expr(ggplot_build(gg)$plot$layers[[1]]$mapping$x))

  x_clades = stringr::str_sub(x, end = 2)
  x_clades = paste0(x_clades, '_rootclade')

  longtiters = gg$data
  yrange = ggplot_build(gg)$layout$panel_scales_y[[1]]$limits

  dat = arrange(distinct(longtiters[, c(x, x_clades)]), x_clades)
  v_lines = tibble(x = numeric(), y = numeric())
  r_prev = dat[[1,x_clades]]
  for (i in seq_along(unlist(dat[,x_clades]))){
    r = dat[[i,x_clades]]
    if (isTRUE(r!=r_prev )){
      v_lines = add_row(v_lines, x = i, y = min(yrange) )
      v_lines = add_row(v_lines, x = i, y = max(yrange) )
    }
    r_prev = r
  }

  gg = gg + geom_vline(xintercept = unique(v_lines$x)-.5, size = .15, color = 'grey20')


  gg = gg +
    geom_tile(
    data = dat,
    aes_string(
      x = x,
      y = mean(yrange),
      fill = x_clades
    ) ,
    width = 1,
    height = diff(yrange) - .1
  ) +
    #geom_line(data = v_lines, aes(x = x, y = y)) +
    # background calde colours
    scale_fill_manual(values = extras.t_col(colors, transparency), name = 'Clade')

  gg$layers = c(gg$layers[[length(gg$layers)]], gg$layers[-length(gg$layers)])

  gg
}

#titerplot.cladeText <- function(gg, )

#'@export
titerscatter.faceted <- function(
  longTiters.wt.unmerged,
  xvar = 'sr',
  facet = 'ag',
  yvar = 'plottiter',
  text_subs = F
  ){

  facet_clades = paste0(facet, '_rootclade')
  longTiters.wt.merged = longtiters.merge(longTiters.wt.unmerged, columns = c('ag', 'sr'))
  longTiters.wt.merged$plottiter = logtiter.toPlot(longTiters.wt.merged$logtiter)

  ylims = range(titers.toNumeric(longTiters.wt.merged[,yvar]))

  if (yvar == 'plottiter'){
    ylims[[2]] = ylims[[2]] + 1 # give some space for text
    Yscale = extras.getYscale(ylims,'logtiters')
  }
  else if (str_detect(yvar, 'colbase')){
    ylims[[2]] = ylims[[2]]
    Yscale = extras.getYscale(ylims,'colbase')

  }

  gg.wt <- ggplot( longTiters.wt.merged)+
    geom_line(aes_string(x = xvar, y = yvar), group = 1) + # merged line
    geom_point(data = longTiters.wt.unmerged, aes_string(x= xvar, y = yvar, group = 1, col = 'titertype'), position = position_dodge2(width = .5), size = .7) +  # unmerged points
    scale_colour_manual(values = c('lessthan' = 'red', 'measured' = 'dark green'), name = 'Titer type') + # colour points
    xlab(xvar) + ylab('Titer') +
    geom_text(aes_string(label = facet_clades), size = 4.4, x = 10.3, y = 8.9, vjust = 1, hjust = 1)
  # mutations text

  gg.wt = titerplot.shadeClades(gg.wt, 'sr', longTiters.wt.merged, yrange)

  gg.wt = titerplot.facet(gg.wt, facet, longTiters.wt.merged)


  if (!isFALSE(text_subs))  gg.wt <- gg.wt + geom_text(aes_string(label = text_subs), size = 4.4, x = .6, y = 8.9, vjust = 1, hjust = 0)


  gg.wt %>% titerplot.styleaxes(x = xvar, Yscale, longTiters.wt.unmerged)

}


titerscatter.coloured = function(longtiters.unmerged, xvar = 'sr', colour = 'ag', yvar = 'plottiter'){

  longtiters.merged = longtiters.merge(longtiters.unmerged, columns = c('ag', 'sr'))
  yrange = range(titers.toNumeric(longtiters.merged[,yvar])) # give some space for text

  if (yvar == 'plottiter'){
    yrange[[2]] = yrange[[2]]+1
    Yscale = extras.getYscale(yrange,'logtiters')
  }
  else if (str_detect(yvar, 'colbase')){
    yrange[[2]] = yrange[[2]]+0
    Yscale = extras.getYscale(yrange,'colbase')

  }

  #longTiters.merged$ag = droplevels(longTiters.merged$ag)
  gg <- ggplot( longtiters.merged)+
    geom_line(aes_string(x = xvar, y = yvar, colour = 'ag', group = colour), position = position_dodge(width = .5)) + # merged line
    geom_point(data = longtiters.unmerged, aes_string(x= xvar, y = yvar, group = 1, colour = 'ag'), position = position_dodge2(width = .5), size = .7) +  # unmerged points
    labs(colour = c('ag' = 'Antigen', 'sr' = 'Serum')[colour]) +
    ylab('Titer')

  Yscale = extras.getYscaleFromTiters(longtiters.unmerged$logtiter)

  gg = titerplot.styleaxes(gg, xvar, Yscale, longtiters.merged)
  gg = titerplot.shadeClades(gg, xvar, longtiters.merged, yrange = c(0,10))


  gg
}



#' Make sequence comparison plot
#'
#' sequences_tibble should have:
#' 1. a row per antigen, with its name and id
#' 2. its sequence
#' 3. the id of its parent
#'
#'@export
sequences_plot = function(sequences,
                          mask = rep(1, length(sequences)),
                          mutagenised_pos = NULL,
                          outtype = 'pdf') {
  sequences_split = sequences %>% str_sub(1, min(str_length(sequences))) %>% str_split('')
  sequences_matrix = do.call(rbind, sequences_split)
  colnames(sequences_matrix) = 1:dim(sequences_matrix)[2]
  rownames(sequences_matrix) = names(sequences)

  important_sites = sort(unique(unlist(h3_sites)))
  important_sites_var = important_sites[important_sites %in% which(apply(sequences_matrix, 2, function(col)
    length(unique(col)) > 1))]
  col_order = c(mutagenised_pos,  important_sites_var[!(important_sites_var %in% mutagenised_pos)])

  sequences_matrix_filtered = sequences_matrix[, col_order]
  sequences_matrix_filtered_masked = sequences_matrix_filtered

  for (i in seq_along(mask)) {
    masking_antigen = mask[[i]]
    if (masking_antigen != i) {
      sites_equal =  (sequences_matrix_filtered_masked[i,] == sequences_matrix_filtered_masked[masking_antigen,])
      sequences_matrix_filtered_masked[i, sites_equal] = ''
    }
  }

  if (outtype == 'pdf') {
    sequences_matrix_filtered_masked %>% data.frame() %>% mutate(id = rownames(sequences_matrix_filtered_masked)) %>%
      pivot_longer(-id, names_to = 'site', values_to = 'aa') %>%
      mutate(site = str_sub(site, 2)) %>%
      tibble.factorize(c('site', 'id')) -> sequences_matrix_filtered_masked_long


    gg.seq.grey <- ggplot(sequences_matrix_filtered_masked_long) +
      geom_tile(aes(x = site, y = id,  fill = aa)) +
      scale_fill_manual(values = c(aacolors_al, '.' = 'grey')) +
      geom_text(aes(x = site, y = id,  label = aa),
                size = 7,
                colour = 'white') +
      theme(legend.position = "none",
            axis.text.x = element_text(size = 11))

    if (length(unique(mask)) > 1 &
        length(unique(mask)) != dim(sequences_matrix_filtered_masked)[[1]]) {
      mask_df = data.frame(
        x = rep(1, dim(sequences_matrix_filtered_masked)[[1]]),
        ag = rownames(sequences_matrix_filtered_masked),
        mask = as.factor(mask)
      )
      line_df = tibble(x = rep(c(
        -.1, dim(sequences_matrix_filtered_masked)[[2]] + .5
      ), length(which(diff(
        mask
      ) > 0))), y = which(diff(mask) > 0) + 0.5)

      gg.seq.grey = gg.seq.grey +
        ggnewscale::new_scale_fill() +
        geom_tile(
          data = mask_df,
          aes(x = 0, y = ag, fill = mask),
          width = .2,
          height = 1
        ) +
        geom_line(
          data = line_df,
          mapping = aes(x = x, y = y),
          color = "white",
          size = 2
        )
    }

    gg.seq.grey = gg.seq.grey +
      theme(axis.text.x = element_text(
        size = 11,
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ))

  }
  else if (outtype == 'html') {
    sequences_matrix_filtered_masked[sequences_matrix_filtered_masked == ''] = ' '


    cols = c(aa.colors(), ' ' = 'white')[sequences_matrix_filtered_masked]

    dim(cols) = dim(sequences_matrix_filtered_masked)
    dimnames(cols) = dimnames(sequences_matrix_filtered_masked)
    cols = cbind(rep('white', dim(cols)[[1]]), cols)

    colnames(sequences_matrix_filtered_masked) <-
      paste0(colnames(sequences_matrix_filtered_masked),
             ifelse(
               map_int(colnames(sequences_matrix_filtered_masked), str_length) == 2,
               " ",
               ""
             ))


    simpleTable(sequences_matrix_filtered_masked,
                col_spec = list(background = cols))


  }

}



#'@export
extras.getYscale <- function(range, kind = 'logtiters', upperthreshold = F, threshold_val = 2, spacing = 1){
  threshold = TRUE
  if (kind == 'logtiters'){
    if (diff(range)%%spacing != 0)
      warning("Spacing does not evenly divide range")

    range[1] = max(range[1], threshold_val)

    pos = seq(range[1], range[2], spacing)
    scale = 10 * 2^pos
    if (threshold) {
      pos = c(range[1] - 1, pos)
      scale = c(paste0("<", scale[[1]]), scale)
    }

    if (upperthreshold){
      pos = c(pos, range[2]+1)
      scale = c(scale, paste0("<", scale[[length(scale)]]))
    }

    return(list(pos = pos, scale = scale))
  }
  if (kind == 'colbase'){
    return(list(pos = seq(range[1], range[2]), scale = seq(range[1], range[2])))
  }

}


#'@export
extras.getYscaleFromTiters <- function(logtiters, kind = 'logtiters', default_range = c(3, 8)){
  logtiters = unlist(logtiters)[!is.na(unlist(logtiters))]
  titers.checkThresholds(logtiters)

  if ('morethan' %in% Racmacs::titer_types(logtiters)) {
    MTthreshold = T
  }
  else MTthreshold = F
  yrange = range(titers.toNumeric(logtiters))
  yrange[[1]] = min(yrange[[1]] , default_range)
  yrange[[2]] = max(yrange[[2]] , default_range)

  threshold_val = max(titers.toNumeric(logtiters)[which(str_sub(logtiters, 1,1) == '<')])

  return(extras.getYscale(yrange, kind, MTthreshold, threshold_val))
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

