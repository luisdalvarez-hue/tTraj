#' getPopulationWaves: Visualize celltypes along pseudotemporal coordinates
#'
#' @description
#' This function isolates each defined celltype along pseudotemporal coordinates to visualize their location
#' and gene expression(s) along a trajectory. Celltype coordinates are applied gaussian smoothing and individually rescaled
#' to follow expression values along runPlotFun()'s pseudotemporal locations.
#' Must first run runPlotFun() to run this function.
#'
#'
#' @param sce A SingleCellExperiment object containing selected cells from the isolated trajectory with scaled pseudotemporal coordinates.
#' @param pseudotime_sce A column of numeric pseudotemporal coordinate values in SingleCellExperiment object.
#' @param celltype A column of character values representing celltypes within the trajectory.
#' @param plot_as_function_values A data frame containing: pseudotemporal locations, features, and smoothed expression values.
#' @param population_smooth_factor Smooth factor to visualize population peaks along the trajectory. Smaller values gives sharper populations peaks and narrower transitions. Larger values gives wider populations peaks and smoother transitions. Default to 0.5.
#' @param matlab_version A character string of active Matlab version (e.g. "R2026a").
#' @param cyt3_script Path to Matlab folder containing cyt3 scripts. Default to "~/Documents/MATLAB".
#' @param output_dir Path to store output files from running plot_as_function. Default to "~/Downloads".
#'
#' @author Luis D. Alvarez
#' @return Attaches columns containing celltype coordinates to plot_as_function data frame.
#' @export

getPopulationWaves <- function(
    sce,
    pseudotime_sce,
    celltype,
    plot_as_function_values,
    population_smooth_factor = 0.5,
    matlab_version,
    cyt3_script = "~/Documents/MATLAB",
    output_dir = "~/Downloads") {

  # Output directory.
  output_directory <- path.expand(output_dir)
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  output_directory <- normalizePath(
    output_directory,
    winslash = "/",
    mustWork = TRUE
  )

  # Define all output files.
  population_input_csv <- file.path(
    output_directory,
    "population_indicator_plot_as_function_input.csv"
  )

  population_raw_csv <- file.path(
    output_directory,
    "population_gaussian_values_raw.csv"
  )

  population_scaled_csv <- file.path(
    output_directory,
    "population_gaussian_values_scaled.csv"
  )

  population_driver_file <- file.path(
    output_directory,
    "run_population_gaussian_smoothing.m"
  )

  population_matlab_log <- file.path(
    output_directory,
    "population_gaussian_smoothing_log.txt"
  )

  # Check if SingleCellExperiment object exist.
  if (missing(sce)) {
    stop(paste0("The object ", as.character(sce), " was not found."))
  }

  # Check if object belongs to SingleCellExperiment class.
  if(!methods::is(sce, "SingleCellExperiment")) {
    stop(paste0("Object ", as.character(sce), " is not SingleCellExperiment object. ",
                "Run class() or str() to check data type."))
  }

  # Check if pseudotime values are scaled between 0 to 1.
  if ((range(sce[[pseudotime_sce]])[1] < -1e-10) || (range(sce[[pseudotime_sce]])[2] > 1 + 1e-10)) {
    stop(base::paste0("Improper pseudotime range:\n",
                      "\n",
                      paste("Min.:", as.character(range(sce[[pseudotime_sce]])[1]), "\n"),
                      paste("Max.:", as.character(range(sce[[pseudotime_sce]])[2]), "\n"),
                      "\n",
                      "Please scale your vector from 0 to 1.",
                      sep = "\n"))
  }

  # Check if pseudotime contains non-numeric values. Stop if so.
  if (any(!is.finite(sce[[pseudotime_sce]]))) {
    stop("The pseudotime vector contains NA, NaN, Inf, or -Inf.")
  }

  # Check if plot_as_function values exist.
  if (missing(plot_as_function_values)) {

    # If user provides missing object, check if the output CSV created
    # from runPlotFun() is still available to use instead.
    plot_as_function_values_csv <- file.path(
      output_directory,
      "plot_as_function_values_from_sce.csv"
    )

    if (!file.exists(plot_as_function_values_csv)) {
      stop(
        "plot_as_function_values was not found in R, ",
        "and this file does not exist:\n",
        plot_as_function_values_csv
      )
    }

    plot_as_function_values <- data.table::fread(
      plot_as_function_values_csv,
      header = TRUE
    )
  }

  #----

  # SingleCellExperiment's pseudotime.
  scaled_pseudotime_vector <- sce[[pseudotime_sce]]

  scaled_pseudotime_vector <- as.numeric(scaled_pseudotime_vector)

  # MATLAB executable.
  matlab_exe <- paste0(
    "C:/Program Files/MATLAB/",
    matlab_version,
    "/bin",
    "/matlab.exe"
  )

  # Parent directory containing cyt3-master.
  matlab_root <- cyt3_script

  # Population curves after independent 0-1 scaling.
  population_scaled_csv <- file.path(
    output_directory,
    "population_gaussian_values_scaled.csv"
  )

  # Apply compatibility patches to the local CYT file.
  patch_cyt_compatibility <- TRUE

  # Extract celltypes from SingleCellExperiment object.
  cell_populations <- as.character(SummarizedExperiment::colData(sce)[[celltype]])

  # Check if all cells contain a celltype.
  cell_ids <- colnames(sce)

  if (length(cell_ids) != length(cell_populations)) {
    stop("Some cells do not contain specified celltypes. Please double-check before running it again.")
  }

  cat(
    "Celltype(s) used for Gaussian population waves:",
    unique(celltype),
    "\n\n"
  )

  gene_plot_values <- as.data.frame(plot_as_function_values, check.names = FALSE)

  # Must have pseudotime values in the CSV.
  if (!"pseudotime" %in% colnames(gene_plot_values)) {
    stop(
      "plot_as_function_values must contain a column ",
      "named 'pseudotime'."
    )
  }

  # CSV's pseudotime values.
  gene_trajectory_grid <- suppressWarnings(as.numeric(gene_plot_values$pseudotime))

  # Check that SingleCellExperiment's pseudotime matches CSV's pseudotime for plotting.
  sce_range <- range(scaled_pseudotime_vector)

  csv_range <- range(gene_trajectory_grid)

  ranges_overlap <-
    max(sce_range[1], csv_range[1]) <=
    min(sce_range[2], csv_range[2])

  if (!ranges_overlap) {
    stop(
      "plot_as_function's pseudotime does not match with SingleCellExperiment's pseudotime. ",
      "Make sure you are using the correct 'plot_as_function_values_from_sce.csv' or SCE object.\n"
    )
  }

  # Same locations to the plot_as_function CSV.
  num_locs <- length(gene_trajectory_grid)

  # Reorder cells based on their position in pseudotime.
  gene_grid_order <- order(gene_trajectory_grid)

  gene_trajectory_grid <- gene_trajectory_grid[gene_grid_order]

  gene_plot_values <- gene_plot_values[gene_grid_order, , drop = FALSE]

  # Sort populations from early to late.
  # Calculate population median pseudotime for ordering.
  population_median_pseudotime <- tapply(
    scaled_pseudotime_vector,
    cell_populations,
    median,
    na.rm = TRUE
  )

  populations_to_plot <- names(sort(population_median_pseudotime))

  population_counts <- table(
    factor(
      cell_populations,
      levels = populations_to_plot
    )
  )

  cat("Populations included:\n")

  print(population_counts)

  # sce assay:
  # rows    = genes
  # columns = cells
  #
  # plot_as_function:
  # rows    = cells
  # columns = genes

  # 1 = cell belongs to that population
  # 0 = cell belongs to another population

  population_indicator_matrix <- vapply(
    populations_to_plot,
    function(current_population) {
      as.numeric(
        cell_populations == current_population
      )
    },
    numeric(length(cell_populations))
  )

  internal_population_names <- sprintf(
    "POP%03d",
    seq_along(populations_to_plot)
  )

  colnames(population_indicator_matrix) <- internal_population_names

  rownames(population_indicator_matrix) <- cell_ids

  # Reorder cells based on their position in pseudotime
  cell_order <- order(scaled_pseudotime_vector)

  scaled_pseudotime_vector <- scaled_pseudotime_vector[cell_order]

  population_indicator_matrix <- population_indicator_matrix[cell_order, , drop = FALSE]

  ordered_cell_ids <- cell_ids[cell_order]

  population_input <- data.frame(
    pseudotime = scaled_pseudotime_vector,
    as.data.frame(
      population_indicator_matrix,
      check.names = FALSE
    ),
    check.names = FALSE
  )

  # Validate paths.
  # Path to Matlab executable.
  matlab_exe <- normalizePath(
    matlab_exe,
    winslash = "/",
    mustWork = TRUE
  )

  # ----

  matlab_root <- normalizePath(
    path.expand(cyt3_script),
    winslash = "/",
    mustWork = TRUE
  )

  plot_candidates <- list.files(
    path = matlab_root,
    pattern = "^plot_as_function\\.m$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(plot_candidates) == 0L) {
    stop(
      "Could not find plot_as_function.m under:\n",
      matlab_root
    )
  }

  # Prefer the copy located in cyt3-master/src.
  normalized_candidates <- gsub(
    "\\\\",
    "/",
    plot_candidates
  )

  preferred_candidates <- plot_candidates[
    grepl(
      "/src/plot_as_function\\.m$",
      normalized_candidates
    )
  ]

  if (length(preferred_candidates) == 1L) {

    cyt_plot_file <- preferred_candidates

  } else if (length(plot_candidates) == 1L) {

    cyt_plot_file <- plot_candidates

  } else {

    stop(
      "Multiple plot_as_function.m files were found.\n",
      "Provide a narrower cyt3_script directory or remove ",
      "duplicate copies:\n\n",
      paste(plot_candidates, collapse = "\n")
    )
  }

  cyt_plot_file <- normalizePath(
    cyt_plot_file,
    winslash = "/",
    mustWork = TRUE
  )

  cat(
    "\nUsing CYT function:\n",
    cyt_plot_file,
    "\n\n"
  )

  patch_cyt_compatibility <- TRUE

  if (patch_cyt_compatibility) {

    original_cyt_code <- readLines(
      cyt_plot_file,
      warn = FALSE
    )

    patched_cyt_code <- original_cyt_code

    patched_cyt_code <- gsub(
      pattern = paste0(
        "z\\s*=\\s*smootht\\(",
        "diffs\\(:,\\s*mi\\),\\s*span\\);"
      ),
      replacement = paste0(
        "z = smoothdata(",
        "diffs(:, mi), 'movmean', span);"
      ),
      x = patched_cyt_code,
      perl = TRUE
    )

    patched_cyt_code <- gsub(
      pattern = paste0(
        "z\\s*=\\s*smooth\\(",
        "diffs\\(:,\\s*mi\\),\\s*span\\);"
      ),
      replacement = paste0(
        "z = smoothdata(",
        "diffs(:, mi), 'movmean', span);"
      ),
      x = patched_cyt_code,
      perl = TRUE
    )

    patched_cyt_code <- gsub(
      pattern = paste0(
        "matColors\\s*=\\s*",
        "distinguishable_colors\\(",
        "size\\(Y,\\s*2\\)\\);"
      ),
      replacement = "matColors = lines(size(Y, 2));",
      x = patched_cyt_code,
      perl = TRUE
    )


    code_changed <- !identical(
      original_cyt_code,
      patched_cyt_code
    )

    if (code_changed) {

      backup_file <- paste0(
        cyt_plot_file,
        ".backup"
      )

      # Create the backup only once so it preserves the original.
      if (!file.exists(backup_file)) {

        backup_created <- file.copy(
          from = cyt_plot_file,
          to = backup_file,
          overwrite = FALSE
        )

        if (!backup_created) {
          stop(
            "Could not create a backup of plot_as_function.m:\n",
            backup_file
          )
        }
      }

      writeLines(
        text = patched_cyt_code,
        con = cyt_plot_file,
        useBytes = TRUE
      )

      cat(
        "CYT compatibility patches were written to:\n",
        cyt_plot_file,
        "\n\nBackup file:\n",
        backup_file,
        "\n\n"
      )

    } else {

      cat(
        "No CYT changes were required; ",
        "plot_as_function.m is already compatible.\n\n"
      )
    }

    verified_cyt_code <- readLines(
      cyt_plot_file,
      warn = FALSE
    )

    remaining_smootht <- any(
      grepl(
        "smootht\\s*\\(",
        verified_cyt_code,
        perl = TRUE
      )
    )

    remaining_distinguishable_colors <- any(
      grepl(
        "distinguishable_colors\\s*\\(",
        verified_cyt_code,
        perl = TRUE
      )
    )

    if (remaining_smootht) {
      stop(
        "The compatibility patch failed: ",
        "smootht() is still present in plot_as_function.m."
      )
    }

    if (remaining_distinguishable_colors) {
      stop(
        "The compatibility patch failed: ",
        "distinguishable_colors() is still present in ",
        "plot_as_function.m."
      )
    }
  }
  # ----

  # Output path.
  output_directory <- normalizePath(
    output_directory,
    winslash = "/",
    mustWork = TRUE
  )

  cyt_repository <- dirname(
    dirname(cyt_plot_file)
  )

  cat(
    "\nUsing CYT function:\n",
    cyt_plot_file,
    "\n\n"
  )

  # Write input for plot_as_function.m.
  data.table::fwrite(
    population_input,
    population_input_csv,
    row.names = FALSE,
    col.names = TRUE
  )

  if (!file.exists(population_input_csv)) {
    stop(
      "The MATLAB population input CSV was not created."
    )
  }

  to_matlab_path <- function(x) {
    gsub(
      "\\\\",
      "/",
      x
    )
  }

  escape_matlab_string <- function(x) {
    gsub(
      "'",
      "''",
      x,
      fixed = TRUE
    )
  }

  population_input_csv_matlab <- escape_matlab_string(
    to_matlab_path(population_input_csv)
  )

  population_raw_csv_matlab <- escape_matlab_string(
    to_matlab_path(population_raw_csv)
  )

  cyt_repository_matlab <- escape_matlab_string(
    to_matlab_path(cyt_repository)
  )

  # ----

  compat_directory <- file.path(matlab_root, "compat")
  compat_path_command <- character(0)

  if (dir.exists(compat_directory)) {
    compat_directory <- normalizePath(
      compat_directory,
      winslash = "/",
      mustWork = TRUE
    )

    compat_path_command <- paste0(
      "addpath('",
      escape_matlab_string(
        to_matlab_path(compat_directory)
      ),
      "', '-begin');"
    )
  }

  # ----

  # Create Matlab driver.
  matlab_driver <- c(

    "restoredefaultpath;",

    compat_path_command,

    paste0(
      "addpath(genpath('",
      cyt_repository_matlab,
      "'), '-begin');"
    ),

    "clear plot_as_function;",
    "rehash toolboxcache;",

    "try",

    paste0(
      "    T = readtable('",
      population_input_csv_matlab,
      "', 'VariableNamingRule', 'preserve');"
    ),

    "    x = T{:, 1};",
    "    Y = T{:, 2:end};",

    paste0(
      "    labels = cellstr(string(",
      "T.Properties.VariableNames(2:end)));"
    ),

    "    labels = labels(:)';",

    paste0(
      "    assert(numel(x) == size(Y, 1), ",
      "'Pseudotime and indicator dimensions do not match.');"
    ),

    "    f = figure('Visible','off','Color','w','Position',[100 100 1400 800]);",

    sprintf(
      paste0(
        "    plot_as_function(",
        "x, Y, ",
        "'avg_type','gaussian', ",
        "'num_locs',%d, ",
        "'smooth',%.16g, ",
        "'normalize',false, ",
        "'svGolay',false, ",
        "'show_error',false);"
      ),
      as.integer(num_locs),
      as.numeric(population_smooth_factor)
    ),

    "    ax = gca;",
    "    curve_handles = findobj(ax, 'Type', 'line');",
    "    curve_handles = flipud(curve_handles(:));",
    "    nPopulations = size(Y, 2);",

    paste0(
      "    assert(numel(curve_handles) == nPopulations, ",
      "'Expected one Gaussian curve per population.');"
    ),

    "    X_plot = get(curve_handles(1), 'XData');",
    "    X_plot = X_plot(:);",
    "    Y_plot = nan(numel(X_plot), nPopulations);",

    "    for j = 1:nPopulations",

    "        xj = get(curve_handles(j), 'XData');",
    "        yj = get(curve_handles(j), 'YData');",

    "        xj = xj(:);",
    "        yj = yj(:);",

    paste0(
      "        assert(numel(xj) == numel(X_plot), ",
      "'Gaussian curves have different lengths.');"
    ),

    paste0(
      "        assert(all(abs(xj - X_plot) < 1e-12), ",
      "'Gaussian curves use different pseudotime positions.');"
    ),

    "        Y_plot(:, j) = yj;",

    "    end",

    "    header = [{'pseudotime'}, labels];",
    "    output_cells = [header; num2cell([X_plot, Y_plot])];",

    paste0(
      "    writecell(output_cells, '",
      population_raw_csv_matlab,
      "');"
    ),

    "    close(f);",

    "catch ME",

    "    fprintf(2, '%s\\n', getReport(ME, 'extended', 'hyperlinks', 'off'));",

    "    if exist('f', 'var') && isgraphics(f)",
    "        close(f);",
    "    end",

    "    exit(1);",

    "end",

    "exit(0);"
  )

  unlink(
    c(
      population_raw_csv,
      population_scaled_csv,
      population_matlab_log
    ),
    force = TRUE
  )

  writeLines(
    matlab_driver,
    population_driver_file,
    useBytes = TRUE
  )

  population_driver_file <- normalizePath(
    population_driver_file,
    winslash = "/",
    mustWork = TRUE
  )

  batch_expression <- paste0(
    "run('",
    escape_matlab_string(population_driver_file),
    "')"
  )

  matlab_log <- suppressWarnings(
    system2(
      command = matlab_exe,
      args = c(
        "-wait",
        "-batch",
        shQuote(batch_expression)
      ),
      stdout = TRUE,
      stderr = TRUE
    )
  )

  matlab_status <- attr(
    matlab_log,
    "status"
  )

  if (is.null(matlab_status)) {
    matlab_status <- 0L
  }

  writeLines(
    matlab_log,
    population_matlab_log,
    useBytes = TRUE
  )

  if (matlab_status != 0L) {
    stop(
      "MATLAB returned a nonzero exit status.\n",
      "Review:\n",
      population_matlab_log
    )
  }

  if (!file.exists(population_raw_csv)) {
    stop(
      "MATLAB did not create the population Gaussian CSV."
    )
  }

  # Import Matlab output.
  population_gaussian_values_raw <- data.table::fread(
    population_raw_csv,
    header = TRUE
  )

  population_gaussian_values_raw[] <- lapply(
    population_gaussian_values_raw,
    as.numeric
  )

  # Rename old labels with population names.
  data.table::setnames(
    population_gaussian_values_raw,
    old = internal_population_names,
    new = populations_to_plot
  )

  population_trajectory_grid <- as.numeric(
    population_gaussian_values_raw$pseudotime
  )

  population_wave_matrix_matlab <- as.matrix(
    population_gaussian_values_raw[ , populations_to_plot, with = FALSE]
  )

  storage.mode(population_wave_matrix_matlab) <- "double"

  colnames(population_wave_matrix_matlab) <- populations_to_plot

  population_grid_order <- order(population_trajectory_grid)

  population_trajectory_grid <- population_trajectory_grid[population_grid_order]

  population_wave_matrix_matlab <- population_wave_matrix_matlab[population_grid_order, , drop = FALSE]

  population_wave_matrix_raw <- vapply(
    seq_along(populations_to_plot),
    function(population_index) {

      stats::approx(
        x = population_trajectory_grid,
        y = population_wave_matrix_matlab[
          ,
          population_index
        ],
        xout = gene_trajectory_grid,
        method = "linear",
        rule = 2,
        ties = "ordered"
      )$y
    },
    numeric(length(gene_trajectory_grid))
  )

  if (is.null(dim(population_wave_matrix_raw))) {
    population_wave_matrix_raw <- matrix(
      population_wave_matrix_raw,
      ncol = 1L
    )
  }

  colnames(population_wave_matrix_raw) <-
    populations_to_plot

  rownames(population_wave_matrix_raw) <- sprintf(
    "%.5f",
    gene_trajectory_grid
  )

  if (any(!is.finite(population_wave_matrix_raw))) {
    stop(
      "The aligned population-wave matrix contains invalid values."
    )
  }

  population_wave_matrix_raw[
    population_wave_matrix_raw < 0 &
      population_wave_matrix_raw > -1e-12
  ] <- 0

  if (
    min(
      population_wave_matrix_raw,
      na.rm = TRUE
    ) < -1e-8
  ) {
    stop(
      "The Gaussian population matrix contains negative values."
    )
  }

  # ----

  # Min-max scaling independently to each population wave.
  population_peaks <- apply(
    population_wave_matrix_raw,
    2,
    max,
    na.rm = TRUE
  )

  population_wave_matrix <- sweep(
    population_wave_matrix_raw,
    MARGIN = 2,
    STATS = population_peaks,
    FUN = "/"
  )

  # Correct negligible floating-point deviations only.
  population_wave_matrix[
    population_wave_matrix < 0 &
      population_wave_matrix > -1e-12
  ] <- 0

  population_wave_matrix[
    population_wave_matrix > 1 &
      population_wave_matrix < 1 + 1e-12
  ] <- 1

  # ----

  final_wave_range <- range(
    population_wave_matrix,
    finite = TRUE
  )

  if (
    final_wave_range[1] < -1e-8 ||
    final_wave_range[2] > 1 + 1e-8
  ) {
    stop(
      "The globally scaled curves fall outside 0-1."
    )
  }


  cat(
    "Maximum scaled value by population:\n"
  )

  print(
    apply(
      population_wave_matrix,
      2,
      max,
      na.rm = TRUE
    )
  )

  scaled_population_output <- data.frame(
    pseudotime = gene_trajectory_grid,
    population_wave_matrix,
    check.names = FALSE
  )

  data.table::fwrite(
    scaled_population_output,
    population_scaled_csv,
    row.names = FALSE,
    col.names = TRUE
  )

  # Confirm pseudotime alignment.
  if (!isTRUE(
    all.equal(
      as.numeric(plot_as_function_values$pseudotime),
      as.numeric(scaled_population_output$pseudotime),
      tolerance = 1e-10
    )
  )) {
    stop(
      "The gene and population pseudotime grids do not match."
    )
  }
  plot_as_function_values <- cbind(
    plot_as_function_values,
    scaled_population_output[populations_to_plot]
  )

  unlink(c(
    population_input_csv,
    population_raw_csv,
    population_scaled_csv,
    population_driver_file,
    population_matlab_log),
    force = TRUE
  )

  return(plot_as_function_values)

}
