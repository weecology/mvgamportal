# Downloads the latest outputs of the sliding window runs
# and the forecast panels in its figures/ subdirectory)
# so the skill score app and other visualizations can be used
# without refitting anything locally.
#
# Existing files in dest_dir with the same names are replaced, including any
# outputs from local runs.

download_outputs <- function(dest_dir = ".",
                             url = "https://data.rc.ufl.edu/pub/ewhite/portal/",
                             figures = TRUE) {
  url <- sub("/?$", "/", url)
  files <- list_remote_files(url, "rds")
  if (length(files) == 0) {
    stop("No .rds files found at ", url)
  }

  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
  # summaries.rds is a few hundred MB, beyond R's default 60 second timeout
  old_timeout <- options(timeout = max(3600, getOption("timeout")))
  on.exit(options(old_timeout))

  for (file in files) {
    dest <- file.path(dest_dir, file)
    # Download to a temporary name first, so a failed download doesn't leave a
    # truncated file in place of a good one
    tmp <- paste0(dest, ".part")
    message("Downloading ", file)
    status <- tryCatch(
      download.file(paste0(url, file), tmp, mode = "wb", quiet = TRUE),
      error = function(e) {
        unlink(tmp)
        stop("Failed to download ", file, ": ", conditionMessage(e))
      }
    )
    if (status != 0) {
      unlink(tmp)
      stop("Failed to download ", file)
    }
    file.rename(tmp, dest)
  }
  downloaded <- file.path(dest_dir, files)

  if (figures) {
    downloaded <- c(downloaded, download_figures(dest_dir, url))
  }

  invisible(downloaded)
}

# There are thousands of small forecast panels, so these are fetched in
# batches of simultaneous downloads rather than one at a time
download_figures <- function(dest_dir, url, batch_size = 50) {
  figures_url <- paste0(url, "figures/")
  files <- list_remote_files(figures_url, "png")
  if (length(files) == 0) {
    warning("No figures found at ", figures_url)
    return(character())
  }

  figures_dir <- file.path(dest_dir, "figures")
  dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
  message("Downloading ", length(files), " figures")

  failed <- character()
  for (batch in split(files, ceiling(seq_along(files) / batch_size))) {
    dests <- file.path(figures_dir, batch)
    # Failures in a batch are warnings, not errors, so they can be collected
    # and reported together at the end
    suppressWarnings(tryCatch(
      download.file(paste0(figures_url, batch), dests,
                    method = "libcurl", mode = "wb", quiet = TRUE),
      error = function(e) NULL
    ))
    # A failed download can leave an empty file behind
    missing <- !file.exists(dests) | file.size(dests) == 0
    unlink(dests[missing])
    failed <- c(failed, batch[missing])
  }
  if (length(failed) > 0) {
    warning(length(failed), " figures failed to download: ",
            paste(head(failed, 10), collapse = ", "),
            if (length(failed) > 10) ", ...")
  }

  file.path(figures_dir, setdiff(files, failed))
}

# Reads the file names out of an Apache directory listing
list_remote_files <- function(url, extension) {
  listing <- paste(readLines(url, warn = FALSE), collapse = "\n")
  pattern <- paste0('href="[^"/?]+\\.', extension, '"')
  files <- regmatches(listing, gregexpr(pattern, listing))[[1]]
  unique(gsub('^href="|"$', "", files))
}
