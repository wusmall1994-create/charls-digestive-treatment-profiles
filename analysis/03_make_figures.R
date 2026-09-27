options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "analysis/03_make_figures.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
res_dir <- file.path(repo_root, "source_data")
fig_dir <- file.path(repo_root, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Figure 1: transparent sample construction and adjacent-wave eligibility.
tiff(file.path(fig_dir, "Fig1.tif"), width = 3300, height = 2100, res = 300,
     compression = "lzw", bg = "white")
par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))

box <- function(x, y, w, h, text, fill = "#F3F6F9", border = "#52708F", cex = 0.88, bold = FALSE) {
  rect(x - w/2, y - h/2, x + w/2, y + h/2, col = fill, border = border, lwd = 1.8)
  text(x, y, text, cex = cex, family = "sans", font = if (bold) 2 else 1)
}
arrow <- function(x0, y0, x1, y1) arrows(x0, y0, x1, y1, length = 0.08, angle = 22, lwd = 1.5, col = "#5F6B76")

text(0.5, 0.965, "Construction of wave-specific profile samples and adjacent-wave transition samples", cex = 1.25, font = 2, family = "sans")

yrs <- c("2011", "2013", "2015", "2018")
disease <- c(3613, 4133, 5348, 5858)
profile <- c(2787, 3431, 4012, 5184)
x <- c(0.14, 0.38, 0.62, 0.86)
for (i in seq_along(x)) {
  box(x[i], 0.79, 0.185, 0.105, paste0(yrs[i], "\nSelf-reported digestive disease\nn = ", format(disease[i], big.mark = ",")), fill = "#E8F1F8", bold = TRUE)
  arrow(x[i], 0.735, x[i], 0.655)
  excluded <- disease[i] - profile[i]
  box(x[i], 0.60, 0.185, 0.11, paste0("Eligible questionnaire profile\nn = ", format(profile[i], big.mark = ","), "\nExcluded n = ", format(excluded, big.mark = ",")), fill = "#F8FAFC")
}

interval_x <- c(0.26, 0.50, 0.74)
intervals <- c("2011-2013", "2013-2015", "2015-2018")
observed <- c(2114, 2708, 3454)
for (i in 1:3) {
  arrow(x[i] + 0.05, 0.535, interval_x[i], 0.43)
  arrow(x[i+1] - 0.05, 0.535, interval_x[i], 0.43)
  box(interval_x[i], 0.35, 0.205, 0.14,
      paste0(intervals[i], " transition sample\nEligible profiles at both waves\nn = ", format(observed[i], big.mark = ",")),
      fill = "#EAF5EE", border = "#4D8060", bold = TRUE)
}

text(0.5, 0.15,
     "Wave-specific samples required age >=45 years, a positive survey weight, a positive self-reported physician diagnosis,\nand a mutually exclusive eligible questionnaire profile. Transition samples additionally required an eligible profile at the next wave.",
     cex = 0.86, family = "sans", col = "#333333")
dev.off()

# Figure 2: interval-specific five-profile transition matrices.
files <- c("transition_matrix_2011-2013.csv", "transition_matrix_2013-2015.csv", "transition_matrix_2015-2018.csv")
labs <- c("2011-2013", "2013-2015", "2015-2018")
canonical <- c("None listed", "TCM without Western", "Western without TCM", "Combined TCM and Western", "Other only")
short <- c("None listed", "TCM without WM", "WM without TCM", "Combined TCM + WM", "Other only")

tiff(file.path(fig_dir, "Fig2.tif"), width = 4200, height = 1900, res = 300,
     compression = "lzw", bg = "white")
par(mfrow = c(1, 3), mar = c(7.2, 7.3, 3.2, 1.4), oma = c(1.0, 0.5, 1.7, 0.5), family = "sans")
pal <- colorRampPalette(c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"))(101)
for (k in seq_along(files)) {
  d <- read.csv(file.path(res_dir, files[k]), check.names = FALSE)
  row_name <- names(d)[1]
  rownames(d) <- d[[row_name]]
  d[[row_name]] <- NULL
  m <- as.matrix(d[canonical, canonical])
  m <- apply(m, 2, as.numeric)
  rownames(m) <- canonical; colnames(m) <- canonical
  image(1:5, 1:5, t(m[5:1, ]), zlim = c(0, 60), col = pal, axes = FALSE, xlab = "", ylab = "", main = labs[k], cex.main = 1.15)
  axis(1, at = 1:5, labels = short, las = 2, cex.axis = 0.76, tick = FALSE)
  axis(2, at = 1:5, labels = rev(short), las = 2, cex.axis = 0.76, tick = FALSE)
  for (i in 1:5) for (j in 1:5) {
    val <- m[6-j, i]
    text(i, j, sprintf("%.1f", val), cex = 0.82, font = ifelse(val >= 35, 2, 1), col = ifelse(val >= 38, "white", "#1A1A1A"))
  }
  graphics::box(col = "#B7C0C8")
  mtext("Later-wave profile", side = 1, line = 5.9, cex = 0.88)
  if (k == 1) mtext("Earlier-wave profile", side = 2, line = 6.0, cex = 0.88)
}
mtext("Weighted row percentages", side = 3, outer = TRUE, cex = 1.25, font = 2, line = 0.25)
mtext("Each cell is the percentage of respondents in an earlier-wave row who reported the later-wave column profile.", side = 1, outer = TRUE, cex = 0.78, line = -0.2)
dev.off()
