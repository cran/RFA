".First.lib" <-
function(lib, pkg)
{
  library.dynam("RFA", package = pkg, lib.loc = lib)
  return(invisible(0))
}

