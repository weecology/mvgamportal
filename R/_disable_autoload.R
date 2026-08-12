# Shiny automatically sources every file in the R/ directory of an app when the
# app starts. app.R sits at the repo root, which makes this project's R/ that
# directory, so without this file launching the skill score app would run the
# analysis scripts in here — several of which read .rds files that are not in
# the repo, or write results. app.R is self-contained and needs none of them.
#
# Shiny only checks that this file exists; its contents are ignored.
