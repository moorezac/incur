import_trackmate_spot_data <- function(
    paths,
    spot_suffix = "-spots.csv",
    time_suffix = "_datetime.csv") {
  # fun to process a single dir
  process_dir <- function(path) {
    message(paste0("processing:\n", path))
    spot_paths <- list.files(
      path = path,
      pattern = spot_suffix,
      full.names = TRUE
    )
    time_paths <- list.files(
      path = path,
      pattern = time_suffix,
      full.names = TRUE
    )

    if (length(spot_paths) != length(time_paths)) {
      message("no. of spot files does not equal no. of time files")
      message("is this an error in segmentation?")
    }

    # these should be ordered, but just to be safe:
    sort_order <- match(
      str_remove(basename(spot_paths), spot_suffix),
      str_remove(basename(time_paths), time_suffix)
    )
    spot_paths <- spot_paths[sort_order]
    time_paths <- time_paths[sort_order]

    spot_paths <- spot_paths[!is.na(spot_paths)]

    # read in both sets of matched data
    spot_list <- map(
      .x = spot_paths,
      .f = function(x) {
        df <- read_csv(
          file = x,
          show_col_types = FALSE,
          skip = 1
        ) |>
          # remove the first and second line
          slice(-c(1:2)) |>
          # else all columns read in as chr
          auto_convert_numerics() |>
          # add in wells
          mutate(
            well = str_extract(
              string = x,
              # VID207_B10_1-spots.csv"
              pattern = "(?:VID\\d+_)([A-Z]\\d{1,2})(?:_)",
              group = 1
            )
          ) |>
          relocate(well)

        # edge case - there is a tibble but no spots
        if (nrow(df) == 0) {
          return(NA)
        } else {
          return(df)
        }
      }
    )

    spot_list <- discard(
      .x = spot_list,
      .p = rlang::is_na
    )

    names(spot_list) <- map_chr(
      .x = spot_list,
      .f = function(x) {
        unique(x$well)
      }
    )

    datetime_list <- map(
      .x = time_paths,
      .f = function(x) {
        read_csv(
          file = x,
          show_col_types = FALSE,
          col_names = c("index", "datetime")
        )
      }
    )
    names(datetime_list) <- str_split(basename(time_paths), "_") |> map_chr(2)

    datetime_list <-
      datetime_list[names(datetime_list) %in% names(spot_list)]

    # add in datetimes
    spot_list <- map2(
      .x = spot_list,
      .y = datetime_list,
      .f = function(x, y) {
        left_join(
          x = x,
          y = y |> rename("T" = index),
          by = "T"
        ) |>
          relocate(datetime) |>
          # remove zero-based index
          mutate(!!sym("T") := !!sym("T") + 1)
      }
    )

    # clean up colnames
    spot_list <- map(
      .x = spot_list,
      .f = function(x) {
        colnames(x) <- colnames(x) |>
          str_to_lower() |>
          # non alphanumeric to underscore
          str_replace_all("[^a-z0-9]+", "_") |>
          # whitespace
          str_trim(side = "both") |>
          # lead/trail underscores
          str_replace_all("^_+|_+$", "") |>
          # consecutive underscores
          str_replace_all("_+", "_")

        x
      }
    )

    spot_list
  }

  # put into a list of lists
  collated_spot_list <- map(
    .x = paths,
    .f = process_dir
  )

  # if there is a mismatch in number of files across dirs
  if (isFALSE(all(map_int(.x = collated_spot_list, .f = length)))) {
    stop("no. of files across dirs is not equal")
  }

  # find the wells across each dir
  common_wells <- reduce(
    .x = map(.x = collated_spot_list, .f = names),
    .f = intersect
  )

  # bind across wells
  final_spot_list <- map(
    .x = common_wells,
    .f = function(x) {
      filtered_list <- map(
        .x = collated_spot_list,
        .f = x
      )
      bind_rows(filtered_list) |>
        arrange(datetime)
    }
  )
  names(final_spot_list) <- common_wells

  final_spot_list
}
