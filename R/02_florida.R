# 02_florida.R
# Cases and deaths averted in Florida adults, each year and averaged over
# 2011-2013 and 2023-2024. Rate reductions from 01_rates.R are applied to the
# Florida population. 95% CIs come from a Monte Carlo simulation.

rates <- read.csv("output/rates.csv")
bands <- c("50-64", "65+")
years <- 2000:2024
periods <- list("2011-2013" = 2011:2013, "2023-2024" = 2023:2024)
draws <- 10000
set.seed(2026)

# Non-bacteremic pneumonia per invasive case, from the trial placebo arm in
# adults 65 and older (Table 2, per-protocol, vaccine-type first episodes):
# 60 pneumonia and 28 invasive episodes. For 50-64 the ratio is scaled by the
# Canadian study's Table 1 counts: non-bacteremic pneumonia over bacteremic
# pneumonia plus invasive disease without pneumonia, 50-64 against 65+.
pneu <- 60
inv <- 28
scale_50 <- (244 / (152 + 64)) / (439 / (190 + 105))


# --- Florida population ------------------------------------------------------

fl <- data.frame()

# 2000-2009. STATE 12 is Florida, SEX 0 is both sexes, AGE 999 is all ages.
st <- read.csv("data-raw/st-est00int-agesex.csv")
st <- st[st$STATE == 12 & st$SEX == 0 & st$AGE <= 100, ]
for (y in 2000:2009) {
  p <- st[[paste0("POPESTIMATE", y)]]
  fl <- rbind(fl, data.frame(band = bands, year = y,
    pop = c(sum(p[st$AGE >= 50 & st$AGE <= 64]), sum(p[st$AGE >= 65]))))
}

# 2010-2019. This file has 5-year age groups, so 50-64 is the sum of three
# rows. Numbers are written with commas, and each year has a total column
# under a header row of years.
raw <- read.csv("data-raw/sc-est2020int-agesex-12.csv", header = FALSE, colClasses = "character")
hdr <- which(apply(raw, 1, function(r) sum(r %in% as.character(2010:2019)) >= 2))[1]
num <- function(label, y) as.numeric(gsub(",", "", raw[raw[, 1] == label, which(raw[hdr, ] == as.character(y))[1]]))
for (y in 2010:2019) {
  fl <- rbind(fl, data.frame(band = bands, year = y,
    pop = c(num(".50 to 54 years", y) + num(".55 to 59 years", y) + num(".60 to 64 years", y),
    num(".65 years and over", y))))
}

# 2020-2024. Vintage 2025 estimates by single year of age. Ages are written
# ".50" or ".85+", and each year has a total column under a header row of years.
census <- list("data-raw/sc-est2025-syasex-12.csv" = 2020:2024)
for (f in names(census)) {
  raw <- read.csv(f, header = FALSE, colClasses = "character")
  hdr <- which(apply(raw, 1, function(r) sum(r %in% as.character(census[[f]])) >= 2))[1]
  body <- raw[grepl("^\\.[0-9]", raw[, 1]), ]
  age <- as.numeric(sub("\\+", "", sub("^\\.", "", body[, 1])))
  for (y in census[[f]]) {
    p <- as.numeric(gsub(",", "", body[, which(raw[hdr, ] == as.character(y))[1]]))
    fl <- rbind(fl, data.frame(band = bands, year = y,
      pop = c(sum(p[age >= 50 & age <= 64]), sum(p[age >= 65]))))
  }
}


# --- Simulation --------------------------------------------------------------
# Each draw redraws total cases from a Poisson, typed cases within them, and
# PCV13-serotype cases within the typed cases. The 1998-1999 baseline is drawn
# once per draw and used for every year.
# Do not move or remove a random draw. It changes every 95% CI after it.

ratio <- rpois(draws, pneu) / rpois(draws, inv)
sim <- list()

for (b in bands) {
  x <- rates[rates$band == b, ]
  py <- x$cases / x$all_rate     # surveillance population, per 100,000
  n <- sapply(x$cases, function(m) rpois(draws, m))
  typ <- matrix(rbinom(length(n), n, rep(x$typed / x$cases, each = draws)), draws)
  v <- matrix(rbinom(length(n), typ, rep(x$vaccine / x$typed, each = draws)), draws)
  # Each draw is multiplied by CDC's published rate over the calculated rate,
  # so the intervals are centered on the published rate.
  rate <- v / typ * n / rep(py, each = draws) * rep(x$vaccine_rate / x$vaccine_rate_calc, each = draws)
  i <- x$year %in% c(1998, 1999)
  base <- rowSums(v[, i]) / rowSums(typ[, i]) * rowSums(n[, i]) / sum(py[i]) * x$vaccine_base[1] / x$vaccine_base_calc[1]
  red <- (base - rate)[, match(years, x$year)]

  pop <- fl$pop[fl$band == b][match(years, fl$year[fl$band == b])]
  ipd <- red * rep(pop, each = draws) / 1e5
  r <- if (b == "65+") ratio else ratio * scale_50
  sim[[b]] <- list(red = red, ipd = ipd, total = ipd * (1 + r), n = n, typ = typ, v = v)
}

# Comparison serotypes and all serotypes. Comparison-serotype cases are drawn
# from the typed cases that are not PCV13 serotypes.

for (b in bands) {
  x <- rates[rates$band == b, ]
  s <- sim[[b]]
  py <- x$cases / x$all_rate
  i <- x$year %in% c(1998, 1999)
  rest <- s$typ - s$v
  comp <- matrix(rbinom(length(rest), rest, rep(x$comparison / (x$typed - x$vaccine), each = draws)), draws)
  comp_rate <- comp / s$typ * s$n / rep(py, each = draws)
  comp_base <- rowSums(comp[, i]) / rowSums(s$typ[, i]) * rowSums(s$n[, i]) / sum(py[i])
  all_rate <- s$n / rep(py, each = draws)
  all_base <- rowSums(s$n[, i]) / sum(py[i])
  pop <- fl$pop[fl$band == b][match(years, fl$year[fl$band == b])]
  cols <- match(years, x$year)
  sim[[b]]$comp_rise <- (comp_rate - comp_base)[, cols]
  sim[[b]]$all_red <- (all_base - all_rate)[, cols]
  sim[[b]]$all_ipd <- sim[[b]]$all_red * rep(pop, each = draws) / 1e5
}


# --- Results by year ---------------------------------------------------------

lo <- function(m) apply(m, 2, quantile, 0.025)
hi <- function(m) apply(m, 2, quantile, 0.975)

by_year <- data.frame()
for (b in bands) {
  x <- rates[rates$band == b, ][match(years, rates$year[rates$band == b]), ]
  pop <- fl$pop[fl$band == b][match(years, fl$year[fl$band == b])]
  r <- if (b == "65+") pneu / inv else pneu / inv * scale_50
  s <- sim[[b]]
  xb <- rates[rates$band == b & rates$year %in% c(1998, 1999), ]
  all_red <- sum(xb$cases) / sum(xb$cases / xb$all_rate) - x$all_rate
  by_year <- rbind(by_year, data.frame(
    band = b, year = years, pop = pop,
    reduction = x$reduction, reduction_lo = lo(s$red), reduction_hi = hi(s$red),
    ipd = x$reduction * pop / 1e5, ipd_lo = lo(s$ipd), ipd_hi = hi(s$ipd),
    total = x$reduction * pop / 1e5 * (1 + r), total_lo = lo(s$total), total_hi = hi(s$total),
    comp_rise = x$comparison_rate - x$comparison_base, comp_rise_lo = lo(s$comp_rise), comp_rise_hi = hi(s$comp_rise),
    all_red = all_red, all_red_lo = lo(s$all_red), all_red_hi = hi(s$all_red),
    all_ipd = all_red * pop / 1e5, all_ipd_lo = lo(s$all_ipd), all_ipd_hi = hi(s$all_ipd)))
}
write.csv(by_year, "output/florida_by_year.csv", row.names = FALSE)


# --- Results averaged over each period ---------------------------------------
# The rows for adults 50 and older add the two age groups within each draw.

ci <- function(v) quantile(v, c(0.025, 0.975))
by_period <- data.frame()
for (p in names(periods)) {
  cols <- match(periods[[p]], years)
  both_ipd <- 0
  both_total <- 0
  both_all_ipd <- 0
  for (b in bands) {
    y <- by_year[by_year$band == b & by_year$year %in% periods[[p]], ]
    s <- sim[[b]]
    red <- rowMeans(s$red[, cols])
    ipd <- rowMeans(s$ipd[, cols])
    total <- rowMeans(s$total[, cols])
    comp_rise <- rowMeans(s$comp_rise[, cols])
    all_red <- rowMeans(s$all_red[, cols])
    all_ipd <- rowMeans(s$all_ipd[, cols])
    both_ipd <- both_ipd + ipd
    both_total <- both_total + total
    both_all_ipd <- both_all_ipd + all_ipd
    by_period <- rbind(by_period, data.frame(
      band = b, period = p, pop = mean(y$pop),
      reduction = mean(y$reduction), reduction_lo = ci(red)[1], reduction_hi = ci(red)[2],
      ipd = mean(y$ipd), ipd_lo = ci(ipd)[1], ipd_hi = ci(ipd)[2],
      total = mean(y$total), total_lo = ci(total)[1], total_hi = ci(total)[2],
      comp_rise = mean(y$comp_rise), comp_rise_lo = ci(comp_rise)[1], comp_rise_hi = ci(comp_rise)[2],
      all_red = mean(y$all_red), all_red_lo = ci(all_red)[1], all_red_hi = ci(all_red)[2],
      all_ipd = mean(y$all_ipd), all_ipd_lo = ci(all_ipd)[1], all_ipd_hi = ci(all_ipd)[2]))
  }
  y <- by_year[by_year$year %in% periods[[p]], ]
  by_period <- rbind(by_period, data.frame(
    band = "50+", period = p, pop = sum(y$pop) / length(periods[[p]]),
    reduction = NA, reduction_lo = NA, reduction_hi = NA,
    ipd = sum(y$ipd) / length(periods[[p]]), ipd_lo = ci(both_ipd)[1], ipd_hi = ci(both_ipd)[2],
    total = sum(y$total) / length(periods[[p]]), total_lo = ci(both_total)[1], total_hi = ci(both_total)[2],
    comp_rise = NA, comp_rise_lo = NA, comp_rise_hi = NA,
    all_red = NA, all_red_lo = NA, all_red_hi = NA,
    all_ipd = sum(y$all_ipd) / length(periods[[p]]), all_ipd_lo = ci(both_all_ipd)[1], all_ipd_hi = ci(both_all_ipd)[2]))
}
row.names(by_period) <- NULL
write.csv(by_period, "output/florida_periods.csv", row.names = FALSE)


# --- Deaths ------------------------------------------------------------------
# Deaths averted are cases averted times case fatality. For invasive disease,
# case fatality is the ABCs death rate over the case rate for the same age
# group and year, all serotypes.

# Non-bacteremic pneumonia deaths by age, from the Canadian study (LeBlanc et
# al, Vaccine 2022). The authors sent the 30-day deaths by age group on
# 9 November 2026. The cases by age are from Table 1.
nbpp_episodes <- c("50-64" = 244, "65+" = 439)
nbpp_deaths <- c("50-64" = 7, "65+" = 31)

agg <- read.csv("data-raw/abcs_aggregate_spn_2026-11-04.csv")
agg$Value <- as.numeric(gsub(",", "", agg$Value))

cfr <- data.frame()
for (b in bands) {
  x <- rates[rates$band == b, ]
  age <- agg[agg$ViewBy == "Age" & grepl(sub("\\+", "", b), agg$ViewBy2), ]
  dr <- age[age$Topic == "Death Rates", ]
  cr <- age[age$Topic == "Case Rates", ]
  ratio <- dr$Value[match(x$year, dr$Year)] / cr$Value[match(x$year, cr$Year)]
  if (any(is.na(ratio))) stop("death or case rate missing for ", b)
  cfr <- rbind(cfr, data.frame(band = b, year = x$year, cfr = ratio))

  cols <- match(years, x$year)
  died <- matrix(rbinom(draws * nrow(x), rep(x$cases, each = draws), rep(ratio, each = draws)), draws)
  cfr_draw <- (died / rep(x$cases, each = draws))[, cols]
  nbpp_draw <- rbinom(draws, nbpp_episodes[[b]], nbpp_deaths[[b]] / nbpp_episodes[[b]]) / nbpp_episodes[[b]]
  s <- sim[[b]]
  sim[[b]]$ipd_deaths <- s$ipd * cfr_draw
  sim[[b]]$all_ipd_deaths <- s$all_ipd * cfr_draw
  sim[[b]]$pneu_deaths <- (s$total - s$ipd) * nbpp_draw
}

# Point estimates, each year, from the results by year above.
dy <- merge(by_year, cfr, by = c("band", "year"))
dy <- dy[order(dy$band, dy$year), ]
dy$ipd_deaths <- dy$ipd * dy$cfr
dy$all_ipd_deaths <- dy$all_ipd * dy$cfr
dy$pneu_deaths <- (dy$total - dy$ipd) * nbpp_deaths[dy$band] / nbpp_episodes[dy$band]

# Averages over each period, with the rows for adults 50 and older adding the
# two age groups within each draw, as above.
deaths <- data.frame()
for (p in names(periods)) {
  cols <- match(periods[[p]], years)
  both <- list(ipd_deaths = 0, all_ipd_deaths = 0, pneu_deaths = 0, total_deaths = 0)
  for (b in c(bands, "50+")) {
    if (b != "50+") {
      s <- sim[[b]]
      d <- list(ipd_deaths = rowMeans(s$ipd_deaths[, cols]),
                all_ipd_deaths = rowMeans(s$all_ipd_deaths[, cols]),
                pneu_deaths = rowMeans(s$pneu_deaths[, cols]))
      d$total_deaths <- d$ipd_deaths + d$pneu_deaths
      for (k in names(both)) both[[k]] <- both[[k]] + d[[k]]
      y <- dy[dy$band == b & dy$year %in% periods[[p]], ]
    } else {
      d <- both
      y <- dy[dy$year %in% periods[[p]], ]
    }
    n <- length(periods[[p]])
    pt <- c(ipd_deaths = sum(y$ipd_deaths) / n, all_ipd_deaths = sum(y$all_ipd_deaths) / n,
            pneu_deaths = sum(y$pneu_deaths) / n)
    pt["total_deaths"] <- pt[["ipd_deaths"]] + pt[["pneu_deaths"]]
    row <- data.frame(band = b, period = p)
    for (k in names(pt)) {
      row[[k]] <- pt[[k]]
      row[[paste0(k, "_lo")]] <- ci(d[[k]])[[1]]
      row[[paste0(k, "_hi")]] <- ci(d[[k]])[[2]]
    }
    deaths <- rbind(deaths, row)
  }
}
write.csv(dy[, c("band", "year", "cfr", "ipd_deaths", "all_ipd_deaths", "pneu_deaths")],
          "output/florida_deaths_by_year.csv", row.names = FALSE)
write.csv(deaths, "output/florida_deaths.csv", row.names = FALSE)


# --- Serotype groups ---------------------------------------------------------
# Rates and Florida cases averted for five serotype groups, for the Figure and
# eFigure 1. Rates are calculated as in 01_rates.R. "other" is the
# all-serotype rate minus the four named groups. A group whose rate rose above
# its 1998-1999 baseline has negative cases averted.

serotypes <- list(
  pcv7 = c("4", "6B", "9V", "14", "18C", "19F", "23F"),
  pcv13_added = c("1", "3", "5", "6A", "7F", "19A"),
  comp_pcv20 = c("8", "10A", "11A", "12F", "15B", "22F", "33F"),
  comp_other = c("7C", "9N", "13", "15A", "15C", "16F", "17F", "20", "23A", "23B",
                 "31", "34", "35B", "35F", "38"))

s <- read.csv("data-raw/abcs_serotype_1998_2024.csv")
names(s) <- c("age", "site", "year", "serotype", "n")
s <- s[s$site == "All_Sites", ]

groups <- data.frame()
for (b in bands) {
  x <- rates[rates$band == b, ]
  sb <- s[s$age == paste("Age", b), ]
  py <- x$cases / x$all_rate
  i <- x$year %in% c(1998, 1999)
  pop <- fl$pop[fl$band == b][match(x$year, fl$year[fl$band == b])]
  r <- list(other = x$all_rate)
  for (g in names(serotypes)) {
    n <- sapply(x$year, function(y) sum(sb$n[sb$year == y & sb$serotype %in% serotypes[[g]]]))
    r[[g]] <- n / x$typed * x$all_rate
    r$other <- r$other - r[[g]]
  }
  # The groups must add back to the rates in rates.csv.
  stopifnot(max(abs(r$pcv7 + r$pcv13_added - x$vaccine_rate_calc)) < 1e-9)
  stopifnot(max(abs(r$comp_pcv20 + r$comp_other - x$comparison_rate)) < 1e-9)
  for (g in c(names(serotypes), "other")) {
    base <- sum(r[[g]][i] * py[i]) / sum(py[i])
    groups <- rbind(groups, data.frame(band = b, year = x$year, group = g, rate = r[[g]],
      base_rate = base, pop = pop, averted = (base - r[[g]]) * pop / 1e5))
  }
}
write.csv(groups, "output/serotype_groups.csv", row.names = FALSE)

# The five groups must add to the all-serotype cases averted in every year.
g <- groups[groups$year >= 2000, ]
net <- tapply(g$averted, g$year, sum)
stopifnot(max(abs(net - tapply(by_year$all_ipd, by_year$year, sum))) < 1e-6)

# The PCV13 groups use calculated rates and the main estimate uses CDC's
# published rate, so they differ by a few cases a year.
g <- g[g$group %in% c("pcv7", "pcv13_added"), ]
cat("PCV13 groups minus main PCV13 estimate, largest difference in any year:",
    round(max(abs(tapply(g$averted, g$year, sum) - tapply(by_year$ipd, by_year$year, sum))), 1),
    "cases\n")
