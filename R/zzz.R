# zzz.R
#
# Package startup and unload functions

.onLoad <- function(libname, pkgname) {

    # # Make list of package parameters and add to global options
    #
    # # filepath of logfile
    # optTrackSig <- list(TrackSig.logfile = logFileName() )
    #
    # # add more options ...
    # optRpt[["nameOfOption"]] <- value
    #
    # optionsToSet <- !(names(optTrackSig) %in% names(options()))
     #
    # if(any(optionsToSet)) {
    #     options(optShi[optionsToSet])
    # }

    invisible()
}


.onAttach <- function(libname, pkgname) {
  # Startup message
  m <- character()
  m[1] <- "\nWelcome to TrackSig.\n"

  packageStartupMessage(paste(m, collapse=""))
}


# .onUnload <- function(libname, pkgname) {
#
# }




# [END]
