#!/usr/bin/env Rscript
# build_dataset_flow.R
# ---------------------------------------------------------------------------
# Exact dataset accounting for every species / ion / source-field analysis.
# Produces a sequential flow that reconciles exactly:
#
#   raw rows = retained + excl_species + excl_ion + excl_missing_source
#              + excl_nonpositive + duplicates
#
# and a set of INDEPENDENT diagnostic counts (missing observed / Ritz /
# wavenumber) that may overlap and therefore are reported separately so that
# overlapping exclusions never create a false total.
#
# Writes:
#   tables_r/statistical_audit/dataset_flow.csv
#   tables_r/statistical_audit/dataset_flow_by_species.csv
#   figures_r/statistical_audit/dataset_flow.png   (and fig01_dataset_flow.png)
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
if (!exists("emp_p", mode = "function")) source(file.path(if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this), "audit_utils.R"))

# data_path_for(): map a species to its CSV in the data root.
data_path_for <- function(species, data_root) {
  file.path(data_root, sprintf("%s_lines.csv", species))
}

# flow_for(): compute the flow record for one species/ion/source.
flow_for <- function(species, ion, source, data_root) {
  path <- data_path_for(species, data_root)
  if (!file.exists(path)) return(NULL)
  raw <- read_nist_csv(path)
  cl <- clean_lines_source(raw, species = species, ion = ion, source = source)
  as.data.frame(cl$flow, stringsAsFactors = FALSE)
}

# build_dataset_flow(): iterate over a data.frame of (species, ion, source).
build_dataset_flow <- function(specs, data_root) {
  rows <- lapply(seq_len(nrow(specs)), function(i) {
    flow_for(specs$species[i], specs$ion[i], specs$source[i], data_root)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  do.call(rbind, rows)
}

# plot_dataset_flow(): a colour-blind-safe stacked accounting bar per species
# (wavenumber source), making the sequential exclusions visible.
plot_dataset_flow <- function(flow) {
  fe_wn <- flow[flow$source == "wavenumber", , drop = FALSE]
  comp_cols <- c("retained", "excl_ion", "excl_missing_source",
                 "excl_nonpositive", "n_duplicates", "excl_species")
  long <- do.call(rbind, lapply(seq_len(nrow(fe_wn)), function(i) {
    data.frame(species = fe_wn$species[i],
               component = comp_cols,
               count = as.numeric(unlist(fe_wn[i, comp_cols])),
               stringsAsFactors = FALSE)
  }))
  long$component <- factor(long$component, levels = rev(comp_cols))
  pal <- c(retained = "#117733", excl_ion = "#999999",
           excl_missing_source = "#DDCC77", excl_nonpositive = "#CC6677",
           n_duplicates = "#88CCEE", excl_species = "#AA4499")
  ggplot2::ggplot(long, ggplot2::aes(x = stats::reorder(species, -count, sum),
                                     y = count, fill = component)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(values = pal, name = "disposition") +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::labs(
      title = "Dataset accounting: raw rows reconciled to retained lines",
      subtitle = "Sequential, mutually exclusive dispositions (wavenumber source). Stacks sum to raw row count.",
      x = "species (ion II)", y = "rows",
      caption = "NIST is the data provider only; no endorsement implied.") +
    theme_audit()
}

main <- function(argv = commandArgs(TRUE)) {
  cfg <- default_audit_config()
  root <- audit_repo_root()
  data_root <- file.path(root, "data")

  species_all <- c("Fe", "Cr", "Mn", "Co", "Ni", "Ti")
  # Fe is reported under all three source representations; neighbours under wn.
  specs <- rbind(
    expand.grid(species = "Fe", ion = 2L,
                source = c("wavenumber", "observed", "ritz"), stringsAsFactors = FALSE),
    expand.grid(species = setdiff(species_all, "Fe"), ion = 2L,
                source = "wavenumber", stringsAsFactors = FALSE)
  )

  flow <- build_dataset_flow(specs, data_root)
  if (is.null(flow) || nrow(flow) == 0L) stop("No dataset flow computed (missing data?)")

  write_table(flow, file.path(root, "tables_r/statistical_audit/dataset_flow.csv"))

  # by-species summary (wavenumber source) for quick reference
  by_sp <- flow[flow$source == "wavenumber",
                c("species", "ion", "n_raw", "n_valid_species_ion",
                  "retained", "n_duplicates", "excl_ion",
                  "missing_observed", "missing_ritz", "missing_wavenumber",
                  "wavenumber_min_cm", "wavenumber_max_cm", "ell_min", "ell_max")]
  write_table(by_sp, file.path(root, "tables_r/statistical_audit/dataset_flow_by_species.csv"))

  p <- plot_dataset_flow(flow)
  save_fig(p, file.path(root, "figures_r/statistical_audit/dataset_flow.png"), width = 9, height = 6)
  save_fig(p, file.path(root, "figures_r/statistical_audit/fig01_dataset_flow.png"), width = 9, height = 6)

  cat(sprintf("[dataset_flow] %d (species,source) rows; all reconcile: %s\n",
              nrow(flow), all(flow$reconciles)))
  if (!all(flow$reconciles)) {
    bad <- flow[!flow$reconciles, c("species", "source")]
    cat("[WARN] non-reconciling rows:\n"); print(bad)
  }
  invisible(flow)
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("build_dataset_flow\\.R$", .invoked_file)) {
  main()
}
