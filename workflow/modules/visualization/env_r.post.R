# Packages conda does not have, installed into rat-visualization-r after the
# conda env is built. Versions match renv.lock.
repos <- "https://cloud.r-project.org"
remotes::install_version("plotthis", version = "0.12.0", repos = repos, upgrade = "never")
remotes::install_version("scplotter", version = "0.7.4", repos = repos, upgrade = "never")
