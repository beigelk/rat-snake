replace_ext <- function(path, new_ext) {
    paste0(tools::file_path_sans_ext(path), ".", new_ext)
}