ExtractMetadata <- function(df) {
  meta <- c(sample.no=df$sample.no[1], sample.name=df$sample.name[1], 
            eco=df$ecozone[1], lat=df$latitude[1], 
            long=df$longitude[1], 
            form=df$life.form[1]
  )
  meta
}
