

parcel_items <- function(df, scales_list){
  out <- list()
  scales_l <- list()
  for(scalename in names(scales_list)){
      #lapply(names(scales_list), function(scalename){
    #scalename <- names(scales_list)[1]
    if(length(scales_list[[scalename]]) > 3){
      tmp_fa <- fa(df[unlist(scales_list[[scalename]])], nfactors = 1)
      tmp_load <- abs(tmp_fa$loadings[,1])
      num_parcels <- floor(length(tmp_load) / 2)
      out_list <- lapply(1:num_parcels, function(x){
        c(names(tmp_load)[order(tmp_load)[x]], names(tmp_load)[order(tmp_load, decreasing = TRUE)[x]])
      })
      if(length(tmp_load) %% 2 > 0){
        out_list <- c(out_list, list(names(tmp_load)[which(!names(tmp_load) %in% unlist(out_list))]))
      }
      names(out_list) <- paste0("P", scalename, 1:length(out_list))
      scales_l[[scalename]] <- names(out_list)
      out <- c(out, out_list)
    } else {
      out <- c(out, scales_list[scalename])
      scales_l[[scalename]] <- scales_list[[scalename]]
    }
  }
  df_out <- data.frame(sapply(out, function(x){
    rowMeans(df[, x, drop = FALSE], na.rm = TRUE)
  }))
  list(parcels = out, scales_list = scales_l, df_parcels = df_out)
}