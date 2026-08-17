#' runPlotFun wrapper: Run plot_as_function algorithm to apply gaussian smoothing on log-normalized single cell data
#'
#' @description
#' User interface to run plot_as_function.m from cyt (Pe'er Lab) in Matlab. Function applies
#' Gaussian smoothing (per gene) on normalized values along scaled pseudotemporal coordinates.
#' Function modifies plot_as_function.m to handle deprecated functions and errors
#' to run on current Matlab versions.
#'
#' @param sce A SingleCellExperiment object containing selected cells from the isolated trajectory with scaled pseudotemporal coordinates.
#' @param pseudotime A column of numeric pseudotemporal coordinate values in SingleCellExperiment object.
#' @param selected_features A character vector of features to apply smoothing.
#' @param num_locs Number of locations to represent expression values across pseudotemporal coordinates. Every location calculates a Gaussian-weighted average of nearby cells' normalized expression. Default to 100L.
#' @param smooth_factor Width of the Gaussian neighborhood used to estimate expression at each location. Smaller values gives sharper peaks and narrower transitions. Larger values gives wider peaks and smoother transitions. Default to 0.5.
#' @param assay Assay name containing normalized data in SingleCellExperiment object. Default to "logcounts".
#' @param matlab_version A character string of active Matlab version (e.g. "R2026a").
#' @param cyt3_script Path to Matlab folder containing cyt3 scripts. Default to "~/Documents/MATLAB".
#' @param output_dir Path to store output files from running plot_as_function. Default to "~/Downloads".
#'
#' @author Luis D. Alvarez
#' @return Input CSV file containing the pseudotime and normalized gene expression values that was used as input for plot_as.function.m in Matlab. Output CSV file containing pseudotime locations and smoothed gene expression values. Log file of plot_as_function.m console output in Matlab.
#' @export


runPlotFun <- function(
    sce,
    pseudotime,
    selected_features,
    num_locs = 100L,
    smooth_factor = 0.5,
    assay = "logcounts",
    matlab_version,
    cyt3_script = "~/Documents/MATLAB",
    output_dir = "~/Downloads") {

  # Check if object exist in global env.
  if (missing(sce)) {
    stop(paste0("The object ", as.character(sce), " was not found."))
  }

  # Check if object belongs to SingleCellExperiment class.
  if(!methods::is(sce, "SingleCellExperiment")) {
    stop(paste0("Object ", as.character(sce), " is not SingleCellExperiment object. ",
                "Run class() or str() to check data type."))
  }

  # Check if pseudotime values are scaled between 0 to 1.
  if ((range(sce[[pseudotime]])[1] < -1e-10) || (range(sce[[pseudotime]])[2] > 1 + 1e-10)) {
    stop(base::paste0("Improper pseudotime range:\n",
                     "\n",
                     paste("Min.:", as.character(range(sce[[pseudotime]])[1]), "\n"),
                     paste("Max.:", as.character(range(sce[[pseudotime]])[2]), "\n"),
                     "\n",
                     "Please scale your vector from 0 to 1.",
                     sep = "\n"))
  }

  # Check if pseudotime contains non-numeric values. Stop if so.
  if (any(!is.finite(sce[[pseudotime]]))) {
    stop("The pseudotime vector contains NA, NaN, Inf, or -Inf.")
  }

  available_assays <- SummarizedExperiment::assayNames(sce)

  # Check if object contains normalized data set.
  if (!assay %in% available_assays) {
    stop(
      "Assay '", assay, "' was not found.\n",
      "Available assays: ",
      paste(available_assays, collapse = ", ")
    )
  }

  # Pseudotime vector, scaled from 0 to 1 (min-max normalization).
  scaled_pseudotime_vector <- sce[[pseudotime]]

  # Name of the new cell-metadata column in SCE object.
  pseudotime_column <- "scaled_pseudotime"

  # Selected genes to smooth along the trajectory.
  genes_to_plot <- selected_features

  # MATLAB executable.
  matlab_exe <- paste0(
    "C:/Program Files/MATLAB/",
    matlab_version,
    "/bin",
    "/matlab.exe"
  )

  # Parent directory containing cyt3-master.
  matlab_root <- cyt3_script

  # Output directory.
  output_directory <- output_dir

  # Apply compatibility patches to the local CYT file.
  patch_cyt_compatibility <- TRUE

  scaled_pseudotime_vector <- as.numeric(scaled_pseudotime_vector)

  # Unnamed pseudotime vector must already be in the same order as sce.
  if (length(scaled_pseudotime_vector) != ncol(sce)) {
    stop(
      "The unnamed pseudotime vector contains ",
      length(scaled_pseudotime_vector),
      " values, but sce contains ",
      ncol(sce),
      " cells."
      )

  # Adding cellular IDs to pseudotime vector.
  names(scaled_pseudotime_vector) <- colnames(sce)

  print("No cellular IDs attached to pseudotime vector. Attaching cellular IDs
        from SingleCellExperiment object to vector.")
  }

  SummarizedExperiment::colData(sce)[[pseudotime_column]] <- scaled_pseudotime_vector

  # Validate paths
  matlab_exe <- normalizePath(
    matlab_exe,
    winslash = "/",
    mustWork = TRUE
  )

  matlab_root <- normalizePath(
    matlab_root,
    winslash = "/",
    mustWork = TRUE
  )

  output_directory <- normalizePath(
    output_directory,
    winslash = "/",
    mustWork = TRUE
  )

  cat(
    "\nMATLAB executable:\n",
    matlab_exe
  )

  plot_candidates <- list.files(
    matlab_root,
    pattern = "^plot_as_function\\.m$",
    recursive = TRUE,
    full.names = TRUE
  )

  # Locate plot_as_function.m
  if (length(plot_candidates) == 0) {
    stop(
      "Could not find plot_as_function.m under:\n",
      matlab_root
    )
  }

  plot_candidates <- normalizePath(
    plot_candidates,
    winslash = "/",
    mustWork = TRUE
  )

  src_candidates <- plot_candidates[
    grepl(
      "/src/plot_as_function\\.m$",
      plot_candidates,
      ignore.case = TRUE
    )]

  if (length(src_candidates) == 1) {

    cyt_plot_file <- src_candidates

  } else if (length(plot_candidates) == 1) {

    cyt_plot_file <- plot_candidates

  } else {

    stop(
      "Multiple plot_as_function.m files were found:\n\n",
      paste(plot_candidates, collapse = "\n"),
      "\n\nSet cyt_plot_file manually or remove older CYT copies."
    )
  }

  cyt_src <- dirname(cyt_plot_file)
  cyt_repository <- dirname(cyt_src)

  cat(
    "\n\nplot_as_function.m:\n",
    cyt_plot_file,
    "\n\n"
  )

  # Patch and update old Cyt3 code to run on current Matlab versions.
  if (patch_cyt_compatibility) {

    original_cyt_code <- readLines(
      cyt_plot_file,
      warn = FALSE
    )

    patched_cyt_code <- original_cyt_code

    # Correct old smootht typo.
    patched_cyt_code <- sub(
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

    # Replace old smooth dependency when present.
    patched_cyt_code <- sub(
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

    # Replace distinguishable_colors dependency.
    patched_cyt_code <- sub(
      pattern = paste0(
        "matColors\\s*=\\s*",
        "distinguishable_colors\\(",
        "size\\(Y,\\s*2\\)\\);"
      ),
      replacement = "matColors = lines(size(Y, 2));",
      x = patched_cyt_code,
      perl = TRUE
    )
  }

  # Verify genes are present in SingleCellExperiment object.
  missing_genes <- setdiff(
    genes_to_plot,
    rownames(sce)
  )

  if (length(missing_genes) > 0) {
    stop(
      "These genes were not found in SingleCellExperiment object:\n",
      paste(missing_genes, collapse = ", ")
    )
  }

  # sce assay:
  # rows    = genes
  # columns = cells
  #
  # plot_as_function:
  # rows    = cells
  # columns = genes

  # Transpose for plot_as_function formatting and extract data frame with only selected genes.
  # "normalized_data" should still retain pseudotime column for reordering.
  normalized_data <- t(
    as.matrix(
      SummarizedExperiment::assay(
        sce,
        assay, # "logcounts"
      )[genes_to_plot, , drop = FALSE])
  )

  storage.mode(normalized_data) <- "double"

  pseudotime_for_plot <- as.numeric(SummarizedExperiment::colData(sce)[[pseudotime_column]])

  names(pseudotime_for_plot) <- colnames(sce)

  # Ensure the expression matrix follows the exact sce cell order and length.
  normalized_data <- normalized_data[colnames(sce), , drop = FALSE]

  stopifnot(
    length(pseudotime_for_plot) == nrow(normalized_data)
  )

  # Remove cells with invalid values if any.
  valid_cells <- is.finite(pseudotime_for_plot) &
    rowSums(!is.finite(normalized_data)) == 0

  expression_for_plot <- normalized_data[valid_cells, , drop = FALSE]

  pseudotime_for_plot <- pseudotime_for_plot[valid_cells]

  # Sort cells by pseudotime.
  cell_order <- order(pseudotime_for_plot)

  pseudotime_for_plot <- pseudotime_for_plot[cell_order]

  expression_for_plot <- expression_for_plot[cell_order, , drop = FALSE]

  cat(
    "Input for plot_as_function:\n",
    nrow(expression_for_plot), " cells x ",
    ncol(expression_for_plot), " genes\n\n"
  )

  # Define input and output files.
  input_csv <- file.path(
    output_directory,
    "plot_as_function_input_from_sce.csv"
  )

  output_values_csv <- file.path(
    output_directory,
    "plot_as_function_values_from_sce.csv"
  )

  output_png <- file.path(
    output_directory,
    "plot_as_function_from_sce.png"
  )

  matlab_driver_file <- file.path(
    output_directory,
    "run_plot_as_function_from_sce.m"
  )

  matlab_log_file <- file.path(
    output_directory,
    "plot_as_function_from_sce_log.txt"
  )

  # Delete old files if present before running plot_as_function in Matlab.
  unlink(
    c(
      output_values_csv,
      output_png,
      matlab_log_file
    ),
    force = TRUE
  )

  # Write input table for Matlab.
  plot_input <- data.frame(
    pseudotime = pseudotime_for_plot,
    as.data.frame(
      expression_for_plot,
      check.names = FALSE
    ),
    check.names = FALSE
  )

  data.table::fwrite(
    plot_input,
    input_csv,
    row.names = FALSE,
    col.names = TRUE
  )

  # Validate output paths.
  input_csv <- normalizePath(
    input_csv,
    winslash = "/",
    mustWork = TRUE
  )

  output_values_csv <- normalizePath(
    output_values_csv,
    winslash = "/",
    mustWork = FALSE
  )

  output_png <- normalizePath(
    output_png,
    winslash = "/",
    mustWork = FALSE
  )

  matlab_driver_file <- normalizePath(
    matlab_driver_file,
    winslash = "/",
    mustWork = FALSE
  )

  matlab_log_file <- normalizePath(
    matlab_log_file,
    winslash = "/",
    mustWork = FALSE
  )

  escape_matlab_string <- function(x) {
    gsub(
      "'",
      "''",
      x,
      fixed = TRUE
    )
  }

  cyt_repository_matlab <- escape_matlab_string(cyt_repository)

  input_csv_matlab <- escape_matlab_string(input_csv)

  output_values_matlab <- escape_matlab_string(output_values_csv)

  output_png_matlab <- escape_matlab_string(output_png)

  # Create Matlab driver.
  matlab_driver <- c(

    "restoredefaultpath;",

    paste0(
      "addpath(genpath('",
      cyt_repository_matlab,
      "'), '-begin');"
    ),

    "clear plot_as_function;",
    "rehash toolboxcache;",

    "try",

    "    resolved_function = which('plot_as_function');",
    "    disp(['plot_as_function: ' resolved_function]);",

    paste0(
      "    T = readtable('",
      input_csv_matlab,
      "', 'VariableNamingRule', 'preserve');"
    ),

    "    x = T{:, 1};",
    "    Y = T{:, 2:end};",

    paste0(
      "    labels = cellstr(string(",
      "T.Properties.VariableNames(2:end)));"
    ),

    paste0(
      "    assert(numel(x) == size(Y, 1), ",
      "'Pseudotime and expression cell counts do not match.');"
    ),

    # Keep the complete figure command on one MATLAB line.
    paste0(
      "    f = figure(",
      "'Visible','off',",
      "'Color','w',",
      "'Position',[100 100 1400 800]);"
    ),

    # Run the actual CYT plot_as_function.m.
    sprintf(
      paste0(
        "    plot_as_function(",
        "x, Y, ",
        "'avg_type','gaussian', ",
        "'num_locs',%d, ",
        "'smooth',%.16g, ",
        "'normalize',false, ",
        "'svGolay',false, ",
        "'show_error',false, ",
        "'labels',labels);"
      ),
      as.integer(num_locs),
      as.numeric(smooth_factor)
    ),

    "    ax = gca;",

    # Extract the exact curves drawn by plot_as_function.
    "    curve_handles = findobj(ax, 'Type', 'line');",
    "    curve_handles = flipud(curve_handles(:));",
    "    nGenes = size(Y, 2);",

    paste0(
      "    assert(numel(curve_handles) == nGenes, ",
      "'Expected one plotted line per gene.');"
    ),

    "    X_plot = get(curve_handles(1), 'XData');",
    "    X_plot = X_plot(:);",

    "    Y_plot = nan(numel(X_plot), nGenes);",

    "    for j = 1:nGenes",

    "        xj = get(curve_handles(j), 'XData');",
    "        yj = get(curve_handles(j), 'YData');",

    "        xj = xj(:);",
    "        yj = yj(:);",

    paste0(
      "        assert(numel(xj) == numel(X_plot), ",
      "'The plotted curves have different lengths.');"
    ),

    paste0(
      "        assert(all(abs(xj - X_plot) < 1e-12), ",
      "'The plotted curves use different pseudotime locations.');"
    ),

    "        Y_plot(:, j) = yj;",

    "    end",

    "    labels = labels(:)';",
    "    header = [{'pseudotime'}, labels];",

    paste0(
      "    output_cells = ",
      "[header; num2cell([X_plot, Y_plot])];"
    ),

    paste0(
      "    writecell(output_cells, '",
      output_values_matlab,
      "');"
    ),

    "    xlabel('Scaled pseudotime');",
    "    ylabel('Normalized expression');",
    "    title('Gaussian-smoothed expression along pseudotime');",
    "    box on;",
    "    grid on;",

    paste0(
      "    exportgraphics(f, '",
      output_png_matlab,
      "', 'Resolution', 300);"
    ),

    "    close(f);",

    "catch ME",

    paste0(
      "    fprintf(2, '%s\\n', ",
      "getReport(ME, 'extended', ",
      "'hyperlinks', 'off'));"
    ),

    "    if exist('f', 'var') && isgraphics(f)",
    "        close(f);",
    "    end",

    "    exit(1);",

    "end",

    "exit(0);"
  )

  # Write Matlab driver file.
  writeLines(
    matlab_driver,
    matlab_driver_file,
    useBytes = TRUE
  )

  if (!file.exists(matlab_driver_file)) {
    stop("The MATLAB driver file was not created.")
  }

  # Run Matlab.
  driver_file_matlab <- escape_matlab_string(
    matlab_driver_file
  )

  batch_expression <- paste0(
    "run('",
    driver_file_matlab,
    "')"
  )

  cat(
    "Running Matlab...",
    sep = "\n"
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
    matlab_log_file,
    useBytes = TRUE
  )

  cat(
    "\nMATLAB exit status:",
    matlab_status,
    "\n\n"
  )

  # Verify Matlab output.
  if (matlab_status != 0L) {
    stop(
      "MATLAB returned a nonzero exit status.\n",
      "Review:\n",
      matlab_log_file
    )
  }

  if (!file.exists(output_values_csv)) {
    stop(
      "MATLAB did not create the numerical output CSV.\n",
      "Review:\n",
      matlab_log_file
    )
  }

  if (!file.exists(output_png)) {
    stop(
      "MATLAB did not create the PNG figure.\n",
      "Review:\n",
      matlab_log_file
    )
  }

  if (file.info(output_values_csv)$size == 0) {
    stop("The numerical output CSV is empty.")
  }

  # Import plot_as_function.m values into R.
  plot_as_function_values <- data.table::fread(
    output_values_csv,
    header = TRUE,
    check.names = FALSE
  )

  plot_as_function_values[] <- lapply(
    plot_as_function_values,
    as.numeric
  )

  cat("Preview of plot_as_function.m output:\n")

  print(
    head(plot_as_function_values[, 1:12], n = 5)
  )

  unlink(c(
    output_png,
    matlab_driver_file),
    force = TRUE
  )

  cat(
    "\n\nInput CSV:\n",
    input_csv,
    "\n\nOutput CSV:\n",
    output_values_csv,
    "\n\nMATLAB log:\n",
    matlab_log_file,
    "\n"
  )

  return(plot_as_function_values)

}
