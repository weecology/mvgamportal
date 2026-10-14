# Copies the outputs of a sliding window run to the public web directory that
# download_outputs() reads from: every .rds in the root and the forecast panels
# in figures/. Run at the end of sliding_window_slurm.sh.
#
# The outputs are copied rather than moved, because the next run reads the
# existing .rds files at startup and adds to them.

publish_outputs <- function(src_dir = ".",
                            dest_dir = "/orange/ewhite/web/public/portal") {
  files <- list.files(src_dir, pattern = "\\.rds$")
  # The input data are prepared locally from portalr, not outputs of a run
  files <- files[!grepl("^data_", files)]
  if (length(files) == 0) {
    stop("No outputs found in ", normalizePath(src_dir))
  }
  # Plots from test_model.R are labelled apart from real runs and not published
  figures <- list.files(file.path(src_dir, "figures"), pattern = "\\.png$")
  figures <- figures[!grepl("_test_[0-9]+\\.png$", figures)]

  copy_files(file.path(src_dir, files), file.path(dest_dir, files))
  copy_files(file.path(src_dir, "figures", figures),
             file.path(dest_dir, "figures", figures))

  message("Published ", length(files), " outputs and ", length(figures),
          " figures to ", dest_dir)
  invisible(c(file.path(dest_dir, files),
              file.path(dest_dir, "figures", figures)))
}

copy_files <- function(from, to) {
  # Copy to a temporary name first, so anyone downloading mid-copy gets the
  # previous version rather than a truncated file
  tmp <- paste0(to, ".part")
  copied <- file.copy(from, tmp, overwrite = TRUE)
  renamed <- file.rename(tmp[copied], to[copied])
  if (!all(copied) || !all(renamed)) {
    unlink(tmp)
    stop("Failed to publish: ",
         paste(basename(c(from[!copied], from[copied][!renamed])),
               collapse = ", "))
  }
}
