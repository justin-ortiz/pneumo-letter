# 01_rates.R
# Invasive pneumococcal disease in adults aged 50-64 and 65 and older, 1998-2024,
# from CDC's Active Bacterial Core surveillance (ABCs).
#
# ABCs publishes two datasets. The serotype dataset counts cases by serotype but
# gives no rates. The aggregate dataset gives rates but not by single serotype.
# The PCV13-serotype rate is CDC's published rate. CDC publishes no rate for
# the comparison serotypes, so we calculated theirs as their share of typed
# cases times the published all-serotype rate for the same age group and year.
# The same calculation for the PCV13 serotypes reproduces CDC's published rate,
# which checks the method. The reduction is the 1998-1999 PCV13-serotype rate
# minus each year's rate.
#
# Untyped cases (coded MISS) are left out of the share, which is the same as
# assuming they had the serotype mix of the typed cases. Non-typeable isolates
# (coded NT) were tested and matched no serotype, so they count as typed cases
# that are neither PCV13 nor comparison serotypes.

pcv13 <- c("1", "3", "4", "5", "6A", "6B", "7F", "9V", "14", "18C", "19A", "19F", "23F")
untyped <- "MISS"
bands <- c("50-64", "65+")
baseline_years <- c(1998, 1999)
min_cases <- 50   # comparison serotypes need at least this many cases, 1998-2024

s <- read.csv("data-raw/abcs_serotype_1998_2024.csv")
a <- read.csv("data-raw/abcs_aggregate_spn_2026-11-04.csv")
names(s) <- c("age", "site", "year", "serotype", "n")
# Only the total for all surveillance sites matches the published rates, which
# are for all sites combined.
s <- s[s$site == "All_Sites", ]
a$Value <- as.numeric(gsub(",", "", a$Value))

rates <- data.frame()
comp <- list()

for (b in bands) {
  sb <- s[s$age == paste("Age", b), ]
  # Rates come from the Serotypes topic, the same table that holds CDC's PCV13
  # rate, so the calculated PCV13 rate can be checked against it. The Case Rates
  # topic differs from it by up to 0.8 per 100,000 in some years.
  # The 65+ label in the aggregate file has a >= sign, so match on the number.
  ab <- a[a$Topic == "Serotypes" & grepl(sub("\\+", "", b), a$ViewBy), ]
  
  # Comparison serotypes: outside PCV13, not 6C, and with at least min_cases
  # cases over 1998-2024. 6C is left out because protection against 6A, which
  # is in PCV13, may also reduce 6C disease. MISS and NT are not serotypes.
  totals <- tapply(sb$n, sb$serotype, sum)
  keep <- totals >= min_cases & !(names(totals) %in% c(pcv13, "6C", "MISS", "NT"))
  comp[[b]] <- names(totals)[keep]
  
  for (y in 1998:2024) {
    sy <- sb[sb$year == y, ]
    r <- ab$Value[ab$Year == y & ab$ViewBy2 == "All Serotypes"]
    p <- ab$Value[ab$Year == y & ab$ViewBy2 == "PCV13 Serotypes"]
    if (length(r) != 1 | length(p) != 1) stop("published rate not found for ", b, " ", y)
    rates <- rbind(rates, data.frame(
      band = b, year = y,
      cases = sum(sy$n),
      typed = sum(sy$n[!(sy$serotype %in% untyped)]),
      nt = sum(sy$n[sy$serotype == "NT"]),
      vaccine = sum(sy$n[sy$serotype %in% pcv13]),
      comparison = sum(sy$n[sy$serotype %in% comp[[b]]]),
      all_rate = r, pcv13_published = p))
  }
}

# The PCV13-serotype rate is CDC's published rate. The same rate calculated
# from the counts is kept as vaccine_rate_calc. It checks the calculation we
# used for the comparison serotypes, for which CDC publishes no rate, and the
# simulation in 02_florida.R draws from it.
rates$vaccine_rate <- rates$pcv13_published
rates$vaccine_rate_calc <- rates$vaccine / rates$typed * rates$all_rate
rates$comparison_rate <- rates$comparison / rates$typed * rates$all_rate

# Baseline: pool 1998 and 1999, weighting each year by its person-years in the
# surveillance areas. Person-years each year are cases divided by the
# all-serotype rate.
for (b in bands) {
  i <- rates$band == b & rates$year %in% baseline_years
  j <- rates$band == b
  py <- rates$cases[i] / rates$all_rate[i]
  pooled <- sum(rates$cases[i]) / sum(py)
  rates$vaccine_base[j] <- sum(rates$pcv13_published[i] * py) / sum(py)
  rates$vaccine_base_calc[j] <- sum(rates$vaccine[i]) / sum(rates$typed[i]) * pooled
  rates$comparison_base[j] <- sum(rates$comparison[i]) / sum(rates$typed[i]) * pooled
}
rates$reduction <- rates$vaccine_base - rates$vaccine_rate

write.csv(rates, "output/rates.csv", row.names = FALSE)

# Numbers to check against the Supplement.
for (b in bands) {
  x <- rates[rates$band == b, ]
  i <- x$year %in% baseline_years
  cat("\n", b, "\n")
  cat("comparison serotypes:", length(comp[[b]]), "\n")
  cat("baseline typed:", sum(x$typed[i]), " untyped (MISS):", sum(x$cases[i] - x$typed[i]),
      " NT:", sum(x$nt[i]), " vaccine:", sum(x$vaccine[i]), "\n")
  cat("baseline PCV13 rate, published:", round(x$vaccine_base[1], 2),
      " calculated:", round(x$vaccine_base_calc[1], 2),
      " 2024 reduction:", round(x$reduction[x$year == 2024], 2), "\n")
  cat("comparison rate, baseline:", round(x$comparison_base[1], 2),
      " 2024:", round(x$comparison_rate[x$year == 2024], 2), "\n")
  cat("largest difference, calculated minus published PCV13 rate:",
      round(max(abs(x$vaccine_rate_calc - x$pcv13_published)), 2), "\n")
}