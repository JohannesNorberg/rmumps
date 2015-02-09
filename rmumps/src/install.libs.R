ff <- "mumpsdrv"
if (file.exists(ff)) {
    dest <- file.path(R_PACKAGE_DIR,paste0('bin',R_ARCH))
    dir.create(dest,recursive=TRUE,showWarnings=FALSE)
    file.copy(ff,dest,overwrite=TRUE)
}
