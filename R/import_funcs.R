#' @title Import spot data as produced from the `merge_incu_channels.py` and `trackmate_cellpose_segment.py` scripts
#' @description Given a set of IncuCyte data that has been produced via the 2 provided Python scripts, import into R
#' @param path A directory containing both spot and datetime data.
#' @param spot_suffix The suffix used to denote the spot data. This should not need to be modified.
#' @param time_suffix The suffix used to denote the datetime data. This should not need to be modified.
#' @return A named list containing spot data per well.
#' @importFrom dplyr bind_rows left_join mutate relocate rename slice
#' @importFrom purrr discard map map_chr map_int map2
#' @importFrom readr read_csv
#' @importFrom rlang is_null sym
#' @importFrom stringr str_extract str_remove str_replace_all str_split str_to_lower str_trim
#' 
import_trackmate_spot_data <- function(path, spot_suffix = "-spots.csv", time_suffix = "_datetime.csv") {
  # fun to process a single dir
  process_dir <- function(path) {
    message(paste0("processing:\n", path))
    spot_paths <- list.files(path, spot_suffix, full.names = TRUE)
    time_paths <- list.files(path, time_suffix, full.names = TRUE)

    if (length(spot_paths) != length(time_paths)) {
      message("no. of spot files does not equal no. of time files")
      message("is this an error in segmentation?")
    }

    # these should be ordered, but just to be safe:
    sort_order <- match(
      stringr::str_remove(basename(spot_paths), spot_suffix),
      stringr::str_remove(basename(time_paths), time_suffix)
    )
    spot_paths <- spot_paths[sort_order]
    time_paths <- time_paths[sort_order]

    spot_paths <- spot_paths[!is.na(spot_paths)]

    # read in both sets of matched data
    spot_list <- purrr::map(spot_paths, function(x) {
      # VID207_B10_1-spots.csv"
      well <- stringr::str_extract(x, "(?:VID\\d+_)([A-Z]\\d{1,2})(?:_)", 1)

      df <- readr::read_csv(x, show_col_types = FALSE, skip = 1)
      # remove the first and second line
      df <- dplyr::slice(df, -c(1:2))
      # else all columns read in as chr
      df <- auto_convert_numerics(df)
      # add in wells
      df <- dplyr::mutate(df, well = well) |> dplyr::relocate(well)

      # edge case - there is a tibble but no spots
      if (nrow(df) == 0) {
        return(NA)
      } else {
        return(df)
      }
    })

    # remove na
    spot_list <- purrr::discard(spot_list, rlang::is_na)

    # add in well names
    names(spot_list) <- purrr::map_chr(spot_list, function(x) {
      unique(x$well)
    })

    # read in datetime data
    datetime_list <- purrr::map(time_paths, function(x) {
      readr::read_csv(x, show_col_types = FALSE, col_names = c("index", "datetime")) |> dplyr::rename("T" = index)
    })
    # add in well names
    names(datetime_list) <- stringr::str_split(basename(time_paths), "_") |> purrr::map_chr(2)

    # filter for which spot data exists
    datetime_list <- datetime_list[names(datetime_list) %in% names(spot_list)]

    # add in datetimes
    spot_list <- purrr::map2(spot_list, datetime_list, function(x, y) {
      df <- dplyr::left_join(x, y, by = "T")
      # remove zero-based index
      df <- dplyr::mutate(df, !!rlang::sym("T") := !!rlang::sym("T") + 1)  
      df <- dplyr::relocate(df, datetime)
      
      return(df)
    })

    # clean up colnames
    spot_list <- purrr::map(spot_list, function(x) {
      colnames(x) <- colnames(x) |>
        stringr::str_to_lower() |>
        # non alphanumeric to underscore
        stringr::str_replace_all("[^a-z0-9]+", "_") |>
        # whitespace
        stringr::str_trim(side = "both") |>
        # lead/trail underscores
        stringr::str_replace_all("^_+|_+$", "") |>
        # consecutive underscores
        stringr::str_replace_all("_+", "_")

      return(x)
    })

    return(spot_list)
  }

  # put into a list of lists
  collated_spot_list <- purrr::map(path, process_dir)

  # if there is a mismatch in number of files across dirs
  if (isFALSE(all(purrr::map_int(collated_spot_list, length)))) {
    stop("no. of files across dirs is not equal")
  }

  # find the wells across each dir
  common_wells <- purrr::reduce(purrr::map(collated_spot_list, names), intersect)

  # bind across wells
  final_spot_list <- purrr::map(common_wells, function(x) {
    filtered_list <- purrr::map(collated_spot_list, x)
    dplyr::bind_rows(filtered_list) |> arrange(datetime)
  })
  names(final_spot_list) <- common_wells

  return(final_spot_list)
}

adjust_trackmate_multi_track <- function(list, time_col = "datetime", multi_spot_fun = list(max = max)) {
  purrr::map(list, function(x) {
    x <- group_by(x, !!rlang::ensym(time_col))
    x <- summarise(x, across(
          .cols = everything(),
          .fns = multi_spot_fun,
          .names = "{.col}"
        ))
    x <- ungroup(x)
    
    return(x)
  })
}
