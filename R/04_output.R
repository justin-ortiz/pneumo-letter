# 04_output.R
# The Table, the Figure, eFigure 1 and the numbers quoted in the text.
# The Table is tab-delimited so it pastes into Word.

rates <- read.csv("output/rates.csv")
by_period <- read.csv("output/florida_periods.csv")
deaths <- read.csv("output/florida_deaths.csv")
deaths_by_year <- read.csv("output/florida_deaths_by_year.csv")
groups <- read.csv("output/serotype_groups.csv")
cov <- read.csv("output/nis_florida_pooled.csv")
cov_us <- read.csv("output/nis_us_pooled.csv")

dash <- "–"
label <- c("50-64" = "50 to 64 years", "65+" = "65 years or older", "50+" = "50 years or older")
period_label <- c("2011-2013" = paste0("Early childhood PCV13 period, 2011", dash, "2013"),
                  "2023-2024" = paste0("Postpandemic period, 2023", dash, "2024"))

# An estimate with its 95% CI. Rates have 1 decimal, cases are rounded to 10
# and deaths to whole numbers.
est <- function(x, lo, hi, f) paste0(f(x), " (", f(lo), dash, f(hi), ")")
rate1 <- function(x) formatC(x, format = "f", digits = 1)
by10 <- function(x) formatC(round(x, -1), format = "d", big.mark = ",")
by1 <- function(x) formatC(round(x), format = "d", big.mark = ",")


# --- Table -------------------------------------------------------------------
# florida_deaths.csv has the same row order as florida_periods.csv.

p <- by_period
stopifnot(all(deaths$band == p$band), all(deaths$period == p$period))
cells <- data.frame(
  "Age group" = label[p$band],
  "IPD rate reduction per 100,000 (95% CI)" =
    ifelse(is.na(p$reduction), "", est(p$reduction, p$reduction_lo, p$reduction_hi, rate1)),
  "IPD cases averted (95% CI)" = est(p$ipd, p$ipd_lo, p$ipd_hi, by10),
  "IPD deaths averted (95% CI)" = est(deaths$ipd_deaths, deaths$ipd_deaths_lo, deaths$ipd_deaths_hi, by1),
  "IPD and nonbacteremic pneumococcal pneumonia cases averted (95% CI)" =
    est(p$total, p$total_lo, p$total_hi, by10),
  "IPD and nonbacteremic pneumococcal pneumonia deaths averted (95% CI)" =
    est(deaths$total_deaths, deaths$total_deaths_lo, deaths$total_deaths_hi, by1),
  check.names = FALSE)

# Each period starts with a row that holds only its label, as in the Table.
table1 <- data.frame()
for (k in names(period_label)) {
  heading <- cells[1, ]
  heading[1, ] <- ""
  heading[1, 1] <- period_label[k]
  table1 <- rbind(table1, heading, cells[p$period == k, ])
}
write.table(table1, "output/table1.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# --- Figure ------------------------------------------------------------------
# Panel A. Cases averted by serotype group. Groups whose rate fell stack above
# zero and groups whose rate rose stack below zero. The black line is the net.
# Panel B. Coverage, plotted cov_shift years before the survey year.

cov_shift <- 2
years <- 2000:2024
group_label <- c(pcv7 = "7 PCV7 serotypes",
                 pcv13_added = "6 serotypes added in PCV13",
                 comp_pcv20 = "7 comparison serotypes in PCV20",
                 comp_other = "15 other comparison serotypes",
                 other = "All other IPD")
# Okabe-Ito colors, distinct for readers with color vision deficiency.
col_groups <- c(pcv7 = "#0072B2", pcv13_added = "#56B4E9", comp_pcv20 = "#D55E00",
                comp_other = "#E69F00", other = "#BBBBBB")
col_us_band <- "#7FCFB6"

# Cases averted by group (rows) and year (columns), both age groups added.
g <- groups[groups$year %in% years, ]
averted <- tapply(g$averted, list(g$group, g$year), sum)[names(group_label), ]
net <- colSums(averted)

draw_figure <- function() {
  layout(matrix(1:3, 3, 1), heights = c(1.5, 0.55, 0.35))
  xr <- c(1997.5, 2024.5)

  # Panel A
  par(mar = c(2.5, 5.5, 2.5, 1))
  up_total <- colSums(pmax(averted, 0))
  down_total <- colSums(pmin(averted, 0))
  yl <- c(min(down_total) * 1.15, max(up_total) * 1.12)
  plot(years, net, type = "n", xlim = xr, ylim = yl, xlab = "", ylab = "", las = 1, bty = "l")
  title(ylab = "Estimated cases averted per year", line = 4)
  title("A   Estimated IPD cases averted, Florida adults aged 50 years or older",
        adj = 0, font.main = 1)
  rect(2019.5, yl[1], 2022.5, yl[2], col = "#FFF3C4", border = NA)
  text(2020.5, yl[1], "Pandemic\nyears", pos = 3)
  rect(1997.5, yl[1], 1999.5, yl[2], col = "grey96", border = NA)
  text(1998.5, yl[2] * 0.5, paste0("Baseline, 1998", dash, "1999"), srt = 90)
  v <- c(2000, 2010, 2022, 2023) - 0.5
  vlab <- c("PCV7", "PCV13", "PCV15", "PCV20")
  segments(v, yl[1], v, yl[2] * 0.98, lty = 2, lwd = 1.5)
  text(v, yl[2] * 0.98, vlab, pos = c(4, 4, 2, 4), offset = 0.2)
  for (j in seq_along(years)) {
    up <- 0
    down <- 0
    for (k in names(group_label)) {
      a <- averted[k, j]
      if (a >= 0) {
        rect(years[j] - 0.42, up, years[j] + 0.42, up + a, col = col_groups[k], border = NA)
        up <- up + a
      } else {
        rect(years[j] - 0.42, down + a, years[j] + 0.42, down, col = col_groups[k], border = NA)
        down <- down + a
      }
    }
  }
  abline(h = 0, col = "grey20")
  lines(years, net, lwd = 2.5)
  points(years, net, pch = 19, cex = 0.6)

  # Panel B. The US line is a wide band under the thin Florida line.
  par(mar = c(4, 5.5, 2.5, 1))
  plot(cov$survey_year - cov_shift, cov$coverage_3plus, type = "n", xlim = xr, ylim = c(0, 100),
       xlab = "", ylab = "", las = 1, bty = "l", yaxt = "n")
  title(ylab = "Coverage (%)", line = 4)
  axis(2, at = c(0, 50, 100), las = 1)
  title("B   PCV coverage in children, 3 or more doses", adj = 0, font.main = 1)
  title(xlab = paste("Coverage plotted", cov_shift, "years before the survey year"), line = 2.6)
  rect(1997.5, 0, 1999.5, 100, col = "grey96", border = NA)
  abline(v = c(2000, 2010) - 0.5, lty = 2, lwd = 1.5)
  text(2001.75, 50, "No public use data\nbefore survey year 2005")
  lines(cov_us$survey_year - cov_shift, cov_us$coverage_3plus, col = col_us_band, lwd = 6)
  lines(cov$survey_year - cov_shift, cov$coverage_3plus, col = "grey20", lwd = 1.5)

  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", bty = "n", ncol = 2, y.intersp = 1.1,
         legend = c(group_label, "Net IPD cases averted", "Florida coverage", "US coverage", "", ""),
         fill = c(col_groups, NA, NA, NA, NA, NA), border = NA,
         lty = c(rep(NA, 5), 1, 1, 1, NA, NA), lwd = c(rep(NA, 5), 2.5, 1.5, 6, NA, NA),
         pch = c(rep(NA, 5), 19, NA, NA, NA, NA),
         col = c(rep(NA, 5), "black", "grey20", col_us_band, NA, NA))
}

# The PNG is for Word. The PDF is a vector file for the journal.
png("output/figure.png", width = 10, height = 9, units = "in", res = 600)
draw_figure()
invisible(dev.off())
cairo_pdf("output/figure.pdf", width = 10, height = 9)
draw_figure()
invisible(dev.off())


# --- eFigure 1 ---------------------------------------------------------------
# IPD rates in the ABCs areas for the four named serotype groups.

col_rates <- c(pcv7 = "#1F4E79", pcv13_added = "#5B9BD5", comp_pcv20 = "#C55A11", comp_other = "#F4B183")
draw_efigure <- function() {
  layout(matrix(c(1, 0, 2, 3, 3, 3), 2, byrow = TRUE), widths = c(1, 0.12, 1), heights = c(1, 0.3))
  par(cex = 0.9)
  panels <- c("50-64" = "A", "65+" = "B")
  for (b in names(panels)) {
    par(mar = c(3, 4.5, 3, 1))
    r <- groups[groups$band == b & groups$group %in% names(col_rates), ]
    top <- max(r$rate)
    plot(r$year, r$rate, type = "n", ylim = c(0, top), xlab = "",
         ylab = "Cases per 100,000", las = 1, bty = "l")
    title(paste0(panels[b], "   Adults aged ", label[b]), adj = 0, font.main = 1)
    abline(v = c(2000, 2010, 2022), col = "grey85")
    text(c(2000, 2010, 2022), top, c("PCV7", "PCV13", "PCV15"), pos = 4, offset = 0.2, cex = 0.8, col = "grey40")
    for (k in names(col_rates)) lines(r$year[r$group == k], r$rate[r$group == k], col = col_rates[k], lwd = 2.5)
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", bty = "n", ncol = 2, lwd = 2.5, col = col_rates, cex = 0.9,
         legend = c("7 PCV7 serotypes", "6 serotypes added by PCV13",
                    "7 comparison serotypes added by PCV20", "15 other comparison serotypes"))
}
png("output/efigure1.png", width = 11, height = 6.5, units = "in", res = 300)
draw_efigure()
invisible(dev.off())
pdf("output/efigure1.pdf", width = 11, height = 6.5)
draw_efigure()
invisible(dev.off())


# --- Numbers quoted in the text ----------------------------------------------

r <- function(b, yr, col) rates[rates$band == b & rates$year == yr, col]
late <- p$period == "2023-2024"
a <- p[p$band == "50+" & late, ]
d <- deaths[deaths$band == "50+" & deaths$period == "2023-2024", ]
cfr <- deaths_by_year[deaths_by_year$year %in% c(2011:2013, 2023:2024), ]

txt <- c(
  paste("Florida population 50+, average 2023-2024:", by1(a$pop)),
  paste("Florida population 65+, average 2011-2013:", by1(p$pop[p$band == "65+" & !late]),
        " 2023-2024:", by1(p$pop[p$band == "65+" & late])),
  paste("All-serotype IPD cases averted, 50+, 2023-2024:", est(a$all_ipd, a$all_ipd_lo, a$all_ipd_hi, by10)),
  paste("All-serotype IPD deaths averted, 50+, 2023-2024:",
        est(d$all_ipd_deaths, d$all_ipd_deaths_lo, d$all_ipd_deaths_hi, by1)),
  paste("Comparison-serotype rate, 50-64: 1998-1999", rate1(r("50-64", 1998, "comparison_base")),
        " 2013", rate1(r("50-64", 2013, "comparison_rate")), " 2023", rate1(r("50-64", 2023, "comparison_rate")),
        " 2024", rate1(r("50-64", 2024, "comparison_rate"))),
  paste("Comparison-serotype rate, 65+: 1998-1999", rate1(r("65+", 1998, "comparison_base")),
        " 2013", rate1(r("65+", 2013, "comparison_rate")), " 2023", rate1(r("65+", 2023, "comparison_rate")),
        " 2024", rate1(r("65+", 2024, "comparison_rate"))),
  paste("All-serotype rate, 65+: 2019", r("65+", 2019, "all_rate"), " 2024", r("65+", 2024, "all_rate")),
  paste("Largest difference, calculated minus published PCV13 rate:",
        formatC(max(abs(rates$vaccine_rate_calc - rates$pcv13_published)), format = "f", digits = 2)),
  paste("IPD case fatality in the analysis periods, 50-64:",
        paste(rate1(100 * range(cfr$cfr[cfr$band == "50-64"])), collapse = " to "), "percent"),
  paste("IPD case fatality in the analysis periods, 65+:",
        paste(rate1(100 * range(cfr$cfr[cfr$band == "65+"])), collapse = " to "), "percent"))
writeLines(txt, "output/text_numbers.txt")

print(table1, row.names = FALSE)
cat("\n")
cat(txt, sep = "\n")
