#' Create Temporary ASCAT Input Files
#'
#' @description
#' Generate temporary TSV files required for ASCAT analysis from processed
#' SNP array data, separating BAF and LogR values.
#'
#' @param data A data frame containing SNP array data with columns: name, chr, pos, baf, lrr.
#' @param filter_chr Character vector of chromosomes to exclude from analysis.
#' @param sample Character string specifying the sample identifier for file naming.
#'
#' @return Invisible. Creates temporary files in the system temp directory.
#'
#' @family copy_number_analysis
#' @importFrom dplyr distinct select filter rename
#' @importFrom purrr walk
#' @importFrom readr write_tsv
#' @importFrom stringr str_c
#' @importFrom rlang sym
#'
#' @examples
#' \dontrun{
#' # Create ASCAT input files
#' make_temp_data(
#'   data = snp_data,
#'   filter_chr = c("X", "Y"),
#'   sample = "tumor_sample_01"
#' )
#' }
make_temp_data <- function(data, filter_chr, sample) {
  purrr::walk(c("baf", "lrr"), function(x) {
    data_filt <- dplyr::distinct(data, name, .keep_all = TRUE)
    data_filt <- dplyr::select(data_filt, name, chr, pos, !!x)
    data_filt <- dplyr::filter(data_filt, !chr %in% filter_chr)
    data_filt <- dplyr::rename(
      data_filt,
      !!sample := !!rlang::sym(x)
    )

    readr::write_tsv(
      data_filt,
      stringr::str_c(tempdir(), "/", sample, "_", x, ".tsv")
    )
  })
}

#' Execute ASCAT Copy Number Analysis Pipeline
#'
#' @description
#' Run the complete ASCAT pipeline for copy number analysis including data loading,
#' optional GC correction, segmentation, and copy number calling.
#'
#' @param test_lrr_path Path to tumor LogR file.
#' @param test_baf_path Path to tumor BAF file.
#' @param germline_lrr_path Path to germline LogR file (optional).
#' @param germline_baf_path Path to germline BAF file (optional).
#' @param filter_chr Character vector of chromosomes to exclude.
#' @param output_dir Directory path for output files and plots.
#' @param opts_correctLogR List of options for LogR correction (optional).
#' @param opts_aspcf List of options for ASPCF segmentation (optional).
#' @param opts_runAscat List of options for ASCAT analysis (optional).
#' @param snp_platform Character string specifying SNP platform for germline prediction.
#'
#' @return A list containing:
#' \describe{
#'   \item{ascat_seg}{ASCAT segmentation object}
#'   \item{ascat_res}{ASCAT results object with copy number calls}
#' }
#'
#' @family copy_number_analysis
#' @importFrom rlang is_null exec
#'
#' @examples
#' \dontrun{
#' # Run ASCAT analysis
#' ascat_results <- execute_ascat(
#'   test_lrr_path = "tumor_lrr.tsv",
#'   test_baf_path = "tumor_baf.tsv",
#'   germline_lrr_path = "normal_lrr.tsv",
#'   germline_baf_path = "normal_baf.tsv",
#'   filter_chr = c("X", "Y"),
#'   output_dir = "ascat_results/"
#' )
#' }
execute_ascat <- function(
  test_lrr_path,
  test_baf_path,
  germline_lrr_path,
  germline_baf_path,
  filter_chr,
  output_dir,
  opts_correctLogR = NULL,
  opts_aspcf = NULL,
  opts_runAscat = NULL,
  snp_platform = NULL
) {
  # Load data
  ascat_obj <- ASCAT::ascat.loadData(
    Tumor_LogR_file = test_lrr_path,
    Tumor_BAF_file = test_baf_path,
    Germline_LogR_file = germline_lrr_path,
    Germline_BAF_file = germline_baf_path
  )

  # Correct GC
  if (!rlang::is_null(opts_correctLogR)) {
    ascat_obj <- rlang::exec(
      ASCAT::ascat.correctLogR,
      ASCATobj = ascat_obj,
      !!!opts_correctLogR
    )
  }

  if (rlang::is_null(germline_lrr_path)) {
    gg <- ASCAT::ascat.predictGermlineGenotypes(
      ascat_obj,
      platform = snp_platform
    )
    opts_aspcf <- append(opts_aspcf, list(ascat.gg = gg))
  }

  # Segment data
  if (!rlang::is_null(opts_aspcf)) {
    ascat_seg <- rlang::exec(
      ASCAT::ascat.aspcf,
      ASCATobj = ascat_obj,
      out.dir = output_dir,
      !!!opts_aspcf
    )
  } else {
    ascat_seg <- ASCAT::ascat.aspcf(ascat_obj, out.dir = output_dir)
  }

  # Run ASCAT
  if (!rlang::is_null(opts_runAscat)) {
    ascat_res <- rlang::exec(
      ASCAT::ascat.runAscat,
      ASCATobj = ascat_seg,
      !!!opts_aspcf
    )
  } else {
    ascat_res <- ASCAT::ascat.runAscat(ascat_seg, img.dir = output_dir)
  }

  list(ascat_seg = ascat_seg, ascat_res = ascat_res)
}

#' Count SNP Probes Within Genomic Segment
#'
#' @description
#' Count the number of SNP probes that fall within a specified genomic segment,
#' used for assessing segment reliability in copy number analysis.
#'
#' @param segment A single row data frame or list containing segment information
#'   with columns: chr, startpos, endpos.
#' @param snp_pos A data frame containing SNP position information with columns:
#'   chr, pos.
#'
#' @return Integer value representing the number of SNP probes within the segment boundaries.
#'
#' @details
#' This function filters SNP positions to the specified chromosome and then
#' counts how many fall within the segment boundaries (inclusive). Segments
#' with very few probes (< 10) may be less reliable for copy number calling.
#'
#' @family copy_number_analysis
#' @importFrom dplyr filter between
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Count probes in a specific segment
#' segment_data <- data.frame(
#'   chr = 1,
#'   startpos = 1000000,
#'   endpos = 2000000
#' )
#'
#' probe_count <- count_probes_in_segment(segment_data, snp_positions)
#' }
count_probes_in_segment <- function(segment, snp_pos) {
  snp_pos <- dplyr::filter(snp_pos, chr == segment$chr)
  snp_pos <- dplyr::filter(
    snp_pos,
    dplyr::between(
      pos,
      segment$startpos,
      segment$endpos
    )
  )

  nrow(snp_pos)
}

process_segments <- function(ascat_seg, ascat_res, sample, output_dir) {
  num_probes <- purrr::map_vec(
    seq_len(nrow(ascat_res$segments)),
    function(x) {
      count_probes_in_segment(
        ascat_res$segments[x, ],
        ascat_seg$SNPpos
      )
    }
  )

  readr::write_csv(
    dplyr::mutate(
      ascat_res$segments,
      num_probes = num_probes,
      ploid = ascat_res$ploidy,
      aberrant_fraction = ascat_res$aberrantcellfraction
    ),
    file.path(output_dir, stringr::str_c(sample, "_segments.csv"))
  )
}

#' Run Complete SNP Array Analysis Pipeline
#'
#' @description
#' Execute ASCAT analysis across multiple test samples with comprehensive
#' error handling, file management, and results processing.
#'
#' @param data_test_list Named list of data frames, each containing test sample
#'   SNP data with required columns: name, chr, pos, baf, lrr.
#' @param data_germline_list Optional named list with germline control data.
#' @param filter_chr Character vector of chromosomes to exclude from analysis.
#' @param output_dir Directory path for output files and plots.
#' @param opts_correctLogR List of LogR correction options.
#' @param opts_aspcf List of ASPCF segmentation options.
#' @param opts_runAscat List of ASCAT execution options.
#' @param snp_platform SNP array platform identifier for when germline is not provided.
#' @param return_ascat Logical indicating whether to return ASCAT objects.
#'
#' @return Either invisible (default) or list of ASCAT results (if return_ascat = TRUE).
#'   Processed segment files are saved to the output directory.
#'
#' @details
#' This comprehensive pipeline function:
#' \itemize{
#'   \item Validates input data structure and required columns
#'   \item Creates temporary input files for ASCAT
#'   \item Executes ASCAT analysis with error handling
#'   \item Processes and saves segment results with metadata
#'   \item Handles multiple samples in batch processing mode
#'   \item Provides detailed error reporting for failed analyses
#' }
#'
#' @family copy_number_analysis
#' @importFrom purrr imap pluck
#' @importFrom stringr str_c str_glue str_flatten_comma
#' @export
#'
#' @examples
#' \dontrun{
#' # Run pipeline on multiple samples
#' run_snp_pipeline(
#'   data_test_list = list(
#'     sample1 = tumor_data_1,
#'     sample2 = tumor_data_2
#'   ),
#'   data_germline_list = list(normal = normal_data),
#'   filter_chr = c("X", "Y"),
#'   output_dir = "copy_number_results/",
#'   return_ascat = FALSE
#' )
#' }
run_snp_pipeline <- function(
  data_test_list,
  data_germline_list = NULL,
  filter_chr = c("X", "Y", "MT", "XY", "0"),
  output_dir,
  opts_correctLogR = NULL,
  opts_aspcf = NULL,
  opts_runAscat = NULL,
  snp_platform = NULL,
  return_ascat = FALSE
) {
  # Required column names
  required <- c("name", "chr", "pos", "baf", "lrr")

  if (!rlang::is_null(data_germline_list)) {
    # Pluck the data
    data_germline <- purrr::pluck(data_germline_list, 1)

    # Check names
    if (!any(required %in% colnames(data_germline))) {
      missing <- required[!required %in% colnames(data_germline)]
      stop(str_c(
        "Required column(s) not found in `data_germline`: ",
        stringr::str_flatten_comma(missing)
      ))
    }

    # All samples to be tested against this
    # Write the germline data
    make_temp_data(
      data = data_germline,
      filter_chr = filter_chr,
      sample = names(data_germline_list)
    )

    germline_lrr_path <- stringr::str_c(
      tempdir(),
      "/",
      names(data_germline_list),
      "_lrr.tsv"
    )
    germline_baf_path <- stringr::str_c(
      tempdir(),
      "/",
      names(data_germline_list),
      "_baf.tsv"
    )
  } else {
    germline_lrr_path <- NULL
    germline_baf_path <- NULL
  }

  # Iterate across each sample to test
  result_list <- purrr::imap(data_test_list, function(x, i) {
    # Check names
    if (!any(required %in% colnames(x))) {
      missing <- required[!required %in% colnames(x)]
      stop(str_c(
        stringr::str_glue("Required column(s) not found in `{i}`: "),
        stringr::str_flatten_comma(missing)
      ))
    }

    # Write the test data
    make_temp_data(
      data = x,
      filter_chr = filter_chr,
      sample = i
    )

    test_lrr_path <- stringr::str_c(tempdir(), "/", i, "_lrr.tsv")
    test_baf_path <- stringr::str_c(tempdir(), "/", i, "_baf.tsv")

    # Execute the pipeline
    result <- try({
      execute_ascat(
        test_lrr_path,
        test_baf_path,
        germline_lrr_path,
        germline_baf_path,
        filter_chr,
        output_dir,
        opts_correctLogR,
        opts_aspcf,
        opts_runAscat,
        snp_platform = snp_platform
      )
    })

    # Remove tempfiles
    unlink(test_lrr_path)
    unlink(test_baf_path)

    if (inherits(result, "try-error")) {
      message(stringr::str_glue("Error processing: {i}"))
    }

    # Process and save segments
    try({
      process_segments(result[[1]], result[[2]], i, output_dir)
    })

    result
  })

  # Remove tempfiles
  unlink(germline_lrr_path)
  unlink(germline_baf_path)

  if (return_ascat) {
    result_list
  } else {
    invisible()
  }
}

#' Calculate Pairwise SNP Correlations with Clustering Visualization
#'
#' @description
#' Calculate pairwise correlations between SNP array samples based on B-allele
#' frequency (BAF) values, perform hierarchical clustering, and create a
#' publication-ready heatmap with dendrograms for sample authentication and
#' quality control.
#'
#' @param data_list A named list of data frames, each representing a sample
#'   with SNP array data. Each data frame must contain a `baf` column with
#'   B-allele frequency values.
#'
#' @return A list containing:
#' \describe{
#'   \item{data}{Data frame with pairwise correlation values between all samples}
#'   \item{heatmap}{ggplot object showing correlation heatmap with dendrograms}
#' }
#'
#' @details
#' This function is essential for SNP array quality control and sample authentication:
#' \itemize{
#'   \item Calculates Pearson correlations between all sample pairs using BAF values
#'   \item Handles missing values with pairwise complete observations
#'   \item Performs hierarchical clustering using (1 - correlation) as distance
#'   \item Creates symmetric correlation matrix with diagonal values = 1
#'   \item Generates clustered heatmap with row and column dendrograms
#'   \item Orders samples by clustering results for optimal visualization
#' }
#'
#' High correlations (> 0.95) typically indicate:
#' \itemize{
#'   \item Same biological sample analyzed multiple times
#'   \item Closely related samples (e.g., tumor/normal pairs)
#' }
#'
#' @family copy_number_analysis
#' @importFrom utils combn
#' @importFrom tibble as_tibble tibble column_to_rownames
#' @importFrom purrr map2_dbl
#' @importFrom dplyr mutate bind_rows rename
#' @importFrom tidyr pivot_wider
#' @importFrom stats cor hclust dist
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient guides guide_colorbar scale_y_discrete theme_minimal theme element_text element_blank unit
#' @importFrom ggdendro ggdendrogram
#' @importFrom patchwork wrap_plots plot_spacer
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate correlations between SNP samples
#' correlation_result <- correlate_snp(snp_data_list)
#'
#' # View correlation data
#' head(correlation_result$data)
#'
#' # Display heatmap
#' print(correlation_result$heatmap)
#'
#' # Identify highly correlated samples (potential duplicates)
#' high_corr <- correlation_result$data[
#'   correlation_result$data$corr > 0.95 &
#'   correlation_result$data$a != correlation_result$data$b,
#' ]
#' }
correlate_snp <- function(data_list) {
  # Get unique combinations
  combinations <- utils::combn(names(data_list), 2)
  combinations <- tibble::as_tibble(t(combinations), .name_repair = \(x) {
    c("a", "b")
  })

  # Correlate
  corr_values <- purrr::map2_dbl(
    combinations$a,
    combinations$b,
    function(a, b) {
      # Required column names
      required <- "baf"
      if (!any(required %in% colnames(data_list[[a]]))) {
        missing <- required[!required %in% colnames(a)]
        stop(str_c(
          stringr::str_glue(
            "Required column(s) not found in {a}: "
          ),
          stringr::str_flatten_comma(missing)
        ))
      }
      if (!any(required %in% colnames(data_list[[b]]))) {
        missing <- required[!required %in% colnames(b)]
        stop(str_c(
          stringr::str_glue(
            "Required column(s) not found in {b}: "
          ),
          stringr::str_flatten_comma(missing)
        ))
      }

      stats::cor(
        data_list[[a]]$baf,
        data_list[[b]]$baf,
        use = "pairwise.complete.obs",
        method = "pearson"
      )
    }
  )
  combinations <- dplyr::mutate(combinations, corr = corr_values)

  # Diagonal
  combinations <- combinations |>
    dplyr::bind_rows(combinations |> dplyr::rename(a = b, b = a)) |>
    dplyr::bind_rows(tibble::tibble(
      a = unique(c(combinations$a, combinations$b)),
      b = unique(c(combinations$a, combinations$b)),
      corr = 1
    ))

  # Matrix
  mat <- combinations |>
    tidyr::pivot_wider(names_from = b, values_from = corr) |>
    tibble::column_to_rownames("a") |>
    as.matrix()

  # Cluster
  hc <- hclust(dist(1 - mat))
  order <- hc$labels[hc$order]

  # Long and order based on cluster
  combinations_l <- as.data.frame(as.table(mat)) |>
    dplyr::rename(a = Var1, b = Var2, corr = Freq) |>
    dplyr::mutate(a = factor(a, levels = order), b = factor(b, levels = order))

  # Plots
  gg_heatmap <- ggplot2::ggplot(
    combinations_l,
    ggplot2::aes(a, b, fill = corr)
  ) +
    ggplot2::geom_tile(colour = "black", linewidth = ggplot2::unit(0.5, "pt")) +
    ggplot2::scale_fill_gradient(low = "white") +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        title = "Correlation",
        frame.colour = "black",
        ticks.colour = "black"
      )
    ) +
    ggplot2::scale_y_discrete(position = "right") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 10),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
      axis.title = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(c(0, 0, 0, 0), "pt")
    )
  gg_dendro_col <- ggdendro::ggdendrogram(
    hc,
    rotate = FALSE,
    theme_dendro = FALSE
  ) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::unit(c(0, 0, 0, 0), "pt"))
  suppressMessages({
    gg_dendro_row <- ggdendro::ggdendrogram(
      hc,
      rotate = TRUE,
      theme_dendro = FALSE
    ) +
      ggplot2::scale_y_reverse() +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = ggplot2::unit(c(0, 0, 0, 0), "pt"))
  })
  gg_final <- patchwork::wrap_plots(
    patchwork::plot_spacer(),
    gg_dendro_col,
    gg_dendro_row,
    gg_heatmap,
    nrow = 2,
    ncol = 2,
    widths = c(0.5, 1, 0.5, 1),
    heights = c(0.5, 1, 1, 0.5)
  )

  list(data = combinations, heatmap = gg_final)
}

#' Create Copy Number Profile Plots from ASCAT Segments
#'
#' @description
#' Generate publication-ready copy number profile plots from ASCAT segmentation
#' results, showing absolute copy number across all chromosomes with event
#' classification and sample metadata.
#'
#' @param path Character string specifying the directory containing ASCAT segment
#'   files (files ending in "_segments.csv").
#'
#' @return A named list of ggplot objects, one for each sample, showing:
#' \itemize{
#'   \item Absolute copy number segments across chromosomes
#'   \item Color-coded copy number events (amplifications, deletions, LOH, etc.)
#'   \item Sample ploidy and tumor purity information
#'   \item Chromosome-wise layout with appropriate scaling
#'   \item Amplification labels for high-level events (> 6 copies)
#' }
#'
#' @details
#' The function automatically:
#' \itemize{
#'   \item Reads all "_segments.csv" files from the specified directory
#'   \item Calculates total copy number (nMajor + nMinor)
#'   \item Classifies copy number events:
#'     \itemize{
#'       \item Neutral: 1+1 (normal diploid)
#'       \item Loss: 1+0 or 0+1 (hemizygous deletion)
#'       \item Deletion: 0+0 (homozygous deletion)
#'       \item LOH: 2+0 or 0+2 (loss of heterozygosity)
#'       \item Gain: total > 2 (low-level amplification)
#'       \item Balanced Gain: equal alleles, total > 2
#'       \item Amplification: total > 6 (high-level amplification)
#'     }
#'   \item Applies appropriate color scheme for each event type
#'   \item Handles amplification labeling with text annotations
#'   \item Scales plots appropriately for chromosome lengths
#' }
#'
#' @family copy_number_analysis
#' @importFrom purrr imap map
#' @importFrom readr read_csv
#' @importFrom dplyr group_by filter mutate row_number n distinct ungroup case_when lead
#' @importFrom stringr str_remove str_c
#' @importFrom ggplot2 ggplot geom_segment scale_x_continuous scale_y_continuous scale_colour_manual labs ggtitle theme element_blank element_rect element_text element_line unit facet_wrap vars
#' @importFrom ggrepel geom_text_repel
#' @export
#'
#' @examples
#' \dontrun{
#' # Create copy number plots from ASCAT results
#' cn_plots <- plot_snp("path/to/ascat/results/")
#'
#' # Display plot for first sample
#' print(cn_plots[[1]])
#'
#' # Save plots
#' purrr::iwalk(cn_plots, function(plot, sample_name) {
#'   ggsave(
#'     filename = paste(sample_name, "copy_number.pdf"),
#'     plot = plot,
#'     width = 12,
#'     height = 6
#'   )
#' })
#' }
plot_snp <- function(path) {
  # Get list of segment files
  paths <- dir(
    path = path,
    pattern = "_segments.csv",
    full.names = TRUE
  )

  # Read in
  data_list <- purrr::map(paths, read_csv, show_col_type = FALSE)
  names(data_list) <- stringr::str_remove(basename(paths), "_segments.csv")

  # Iterate over each
  purrr::imap(data_list, function(x, i) {
    # Abs. number
    x <- dplyr::mutate(x, total = nMajor + nMinor)

    # Size of each chromosome
    x <- dplyr::group_by(x, chr)
    chr_size <- dplyr::filter(x, dplyr::row_number() %in% c(1, dplyr::n()))
    chr_size <- dplyr::mutate(chr_size, start = startpos[1], end = endpos[2])
    # Single chromosome segment
    chr_size <- dplyr::mutate(
      chr_size,
      end = dplyr::case_when(
        is.na(end) ~ endpos[1],
        .default = end
      )
    )
    chr_size <- dplyr::distinct(chr_size, chr, start, end)
    chr_size <- dplyr::ungroup(chr_size)
    chr_size <- dplyr::mutate(
      chr_size,
      length = end - start,
      # Relative to chromosome 1
      ratio = length / length[1]
    )

    # Labels for amplifications
    x <- dplyr::mutate(
      x,
      label = dplyr::case_when(
        nMajor + nMinor > 6 ~ as.character(nMajor + nMinor),
        .default = NA
      )
    )
    x <- dplyr::mutate(
      x,
      prev_y = dplyr::lead(total),
      prev_y = dplyr::case_when(
        !is.na(label) ~ prev_y,
        .default = NA
      )
    )

    # Labels for events
    x <- dplyr::mutate(
      x,
      Event = dplyr::case_when(
        total > 6 ~ "Amplification",
        nMajor == 1 & nMinor == 1 ~ "Neutral",
        nMajor == 2 & nMinor == 0 | nMajor == 0 & nMinor == 2 ~ "LOH",
        nMajor == 1 & nMinor == 0 | nMajor == 0 & nMinor == 1 ~ "Loss",
        nMajor + nMinor == 0 ~ "Deletion",
        nMajor == nMinor & total > 2 ~ "Balanced Gain",
        total > 2 ~ "Gain"
      )
    )
    # Change this for ease in plots
    x <- dplyr::mutate(
      x,
      total_plot = dplyr::case_when(
        total > 6 ~ 2,
        .default = total
      )
    )

    # Colours
    event_colours <- c(
      "#ff924c",
      "#1982c4",
      "#ffca3a",
      "#ff595e",
      "#8ac926",
      "#013c58",
      "black"
    )
    names(event_colours) <- c(
      "Amplification",
      "Deletion",
      "Balanced Gain",
      "Gain",
      "LOH",
      "Loss",
      "Neutral"
    )
    x$Event <- factor(x$Event, levels = names(event_colours))
    unique_events <- unique(x$Event)

    # Plot
    ggplot2::ggplot() +
      ggplot2::geom_segment(
        x,
        mapping = ggplot2::aes(
          startpos,
          xend = endpos,
          y = total_plot,
          colour = Event
        ),
        linewidth = 2.5
      ) +
      # Add lines for amplification events - difficult to plot otherwise
      ggrepel::geom_text_repel(
        dplyr::filter(x, total > 6),
        mapping = ggplot2::aes(
          label = label,
          x = (startpos + endpos) / 2,
          y = prev_y,
          colour = Event
        ),
        colour = event_colours["Amplification"],
        force_pull = 0,
        nudge_x = 1,
        nudge_y = 2,
        segment.size = 0.25,
        segment.alpha = 0.5,
        size = 8 / .pt,
        na.rm = TRUE
      ) +
      ggplot2::scale_x_continuous(expand = c(0, 0), position = "top") +
      ggplot2::scale_y_continuous(
        breaks = 0:6,
        limits = c(c(-0.2, 6.2))
      ) +
      ggplot2::scale_colour_manual(values = event_colours) +
      ggplot2::labs(x = "Chromosome", y = "Abs. CN") +
      ggplot2::ggtitle(
        stringr::str_c(
          stringr::str_c("Sample: ", i),
          stringr::str_c("Ploid: ", round(unique(x$ploid), 2)),
          stringr::str_c("Tumour Purity: ", unique(x$aberrant_fraction)),
          sep = "\n"
        )
      ) +
      ggplot2::theme(
        axis.ticks.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        legend.key = ggplot2::element_rect(fill = NA),
        legend.key.spacing.y = ggplot2::unit(-10, "pt"),
        legend.text = ggplot2::element_text(size = 8),
        legend.title = ggplot2::element_blank(),
        panel.spacing.x = ggplot2::unit(0, "pt"),
        panel.spacing.y = ggplot2::unit(0, "pt"),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(
          colour = "black",
          linetype = "dashed",
          linewidth = 0.125
        ),
        panel.border = ggplot2::element_rect(
          color = "black",
          fill = NA,
          linewidth = 0.5
        ),
        plot.title = ggplot2::element_text(size = 8),
        strip.background = ggplot2::element_rect(
          color = "black",
          linewidth = 0.5
        ),
        strip.text = ggplot2::element_text(size = 8),
        text = ggplot2::element_text(size = 8)
      ) +
      ggplot2::facet_wrap(
        ggplot2::vars(chr),
        nrow = 1,
        scales = "free_x"
      )
  })
}
