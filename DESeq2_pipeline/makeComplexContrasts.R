makeComplexContrasts <- function(cmplx,
                                 contrastDT,
                                 useComplexNames = FALSE,
                                 export.as = "df", 
                                 save2file = FALSE,
                                 path = getwd(), 
                                 filename = "complex_contrasts.xlsx"
                                ) {
   
   
   #contrastDT <- data.table(resNames = gsub(columns, "", resultsNames(dds) ))
   setDT(contrastDT)
   dsgn.para <- colnames(contrastDT)
   contrastDT[, (dsgn.para):= lapply(.SD, as.factor), .SDcols = dsgn.para]
   
   for (idx in seq(length(cmplx))) {
      # make all elements in each list element equal (by NA padding)
      cmplxx <- cmplx[[idx]] %>%
         lapply(., `length<-`, max(lengths(.))) %>%
         as.data.table(.)
      
      print(cmplxx)
      
      nullFC <- setdiff(c("leftFC", "rightFC"), names(cmplxx))
      
      if (length(nullFC) == 1 & "rightFC" %in% nullFC) {
         message("rightFC is missing. Even-weighting on the right side of contrast formula.")
         cmplxx[, rightFC := 1 / uniqueN(cmplxx$right)]
         if (prod(sign(cmplxx$leftFC), na.rm = TRUE) < 0) {
            message("leftFC has an odd number of negative factors (delta comparison!).\nFor delta-delta comparison both leftFC and rightFC must be provided.")
            stop()
         }
      }
      
      if (length(nullFC) == 1 & "leftFC" %in% nullFC) {
         message("leftFC is missing. Even weighting on the left side of contrast formula.")
         cmplxx[, leftFC := 1 / uniqueN(cmplxx$left)]
         if (prod(sign(cmplxx$rightFC), na.rm = TRUE) < 0) {
            message("rightFC has an odd number of negative factors (delta comparison!).\nFor delta-delta comparison both leftFC and rightFC must be provided.")
            stop()
         } else if (prod(sign(cmplxx$rightFC), na.rm = TRUE) > 0 & cmplxx$rightFC[1] > 0) {
            message("Provided rightFC are all posetive. It is assumed that the format in mind was -(rightFC.1 + rightFC.2 + ...)")
            message("All rightFC coeffients will therefore be converted to negative! Please check out the final contrast matrix.")
         }
      }
      
      if (length(nullFC) == 0) {
         if (any(c(prod(sign(cmplxx$rightFC), na.rm = TRUE), prod(sign(cmplxx$leftFC), na.rm = TRUE)) < 0)) {
            message("leftFC and rightFC are both provided with negative factors. A delta-delta comparison!")
         }
      }
      
      if (length(nullFC) == 2) {
         cmplxx[, leftFC := 1 / sum(!is.na(cmplxx$left))]
         cmplxx[, rightFC := 1 / sum(!is.na(cmplxx$right))]
      }
      
      
      
      cmplxx[!is.na(left), ll := paste(round(leftFC,2), left, sep = "*")]
      cmplxx[!is.na(right), rr := paste(round(rightFC,2), right, sep = "*")]
      
      ## if even weighting, reomve 1* from string
      if (abs(prod(sign(cmplxx$leftFC), na.rm = TRUE)) == 1)  cmplxx[, ll := gsub("1\\*", "", ll)]
      if (abs(prod(sign(cmplxx$rightFC), na.rm = TRUE)) == 1) cmplxx[, rr := gsub("1\\*", "", rr)]
      
      for (ir in 2:nrow(cmplxx)) {
         if (substring(cmplxx$ll[ir], 1, 1) != "-" & uniqueN(cmplxx$ll) != 1 & !is.na(cmplxx$ll[ir])) cmplxx$ll[ir] <- paste0("+", cmplxx$ll[ir])
         if (substring(cmplxx$rr[ir], 1, 1) != "-" & uniqueN(cmplxx$rr) != 1 & !is.na(cmplxx$rr[ir])) cmplxx$rr[ir] <- paste0("+", cmplxx$rr[ir])
      }
      
      contName <- paste(
         paste(unique(cmplxx$ll), collapse = ""),
         paste(unique(cmplxx$rr), collapse = ""),
         sep = " vs. "
      )
      
      cmplxx[, c("rr", "ll") := c(NULL, NULL)]
      
      if (!useComplexNames) contName <- names(cmplx)[idx] ## whethwe or not use long names for complex contrast
      
      
      cnt.dummy <- merge(cmplxx[!is.na(right), .(right, rightFC)],
                         cmplxx[!is.na(left), .(left, leftFC)],
                         by.x = "right", by.y = "left",
                         all = TRUE) %>%
         .[is.na(rightFC) , rightFC := 0] %>%
         .[is.na(leftFC) , leftFC := 0] %>%
         .[, cnt.dummy := leftFC - rightFC] %>%
         .[, c("leftFC", "rightFC") := c(NULL, NULL)] %>%
         .[,.(cnt.dummy = sum(cnt.dummy)), by = right] %>%
         setnames(., "cnt.dummy", contName)
      
      if (round(sum(cnt.dummy[,2]),6) != 0 ) {
         print(cnt.dummy)
         res <- readline(prompt = "Complex contrast factors do not sum to 0. Do you want to continue? [y/n]:")
         if (res %in% c("n","N","NO")) stop()
      }
      
      contrastDT %<>% merge(., cnt.dummy, by.x = dsgn.para, by.y = "right", all = TRUE) %>%
         .[is.na(get(contName)), (contName) := 0]
   }
   
   message(
      "\n\n!!!WARNING...\nDOUBLE-CHECK THE COMPLEX CONTRASTS BELLOW. ALSO SAVED IN:\n",
      file.path(paste0(path, filename)), "\n\n"
   )
   print(contrastDT)
   
   ## write contrasts into file
   if (save2file) {
      write.xlsx(
         as.data.table(contrastDT, keep.rownames = "Groups"),
         file.path(path, filename),
         overwrite = TRUE
      )
   }
   
   
   if (tolower(export.as) %in% c("data.frame", "df")) {
      rwnm <- as.character(unlist(contrastDT[, 1]))
      contrastDT <- as.data.frame(contrastDT[, -1])
      rownames(contrastDT) <- rwnm
      return(contrastDT)
   } else if (tolower(export.as) %in% c("data.table", "dt")) {
      return(contrastDT)
   }
}
