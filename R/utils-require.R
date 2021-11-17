requireHere <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf('Package "%s" needed for this function to work. Please install it.',
                 pkg),
         call. = FALSE)
  }
}
