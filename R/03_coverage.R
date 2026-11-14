# 03_coverage.R
# PCV coverage (3 or more doses) in children aged 19 to 35 months, Florida and
# US, from the National Immunization Survey public use files, 2005 to 2024.
# Each survey year needs two files in data-raw/nis: the data file
# (NISPUFyy.DAT) and CDC's program for reading it (NISPUFyy.R or NISPUFyy.sas).

nis_dir <- "data-raw/nis"


# --- Which years can be read ------------------------------------------------
# File names are matched without regard to case.

all_files <- list.files(nis_dir)

find_file <- function(name) {
  hit <- all_files[tolower(all_files) == tolower(name)]
  if (length(hit) == 0) return(NA_character_)
  file.path(nis_dir, hit[1])
}

dat_files  <- grep("^nispuf[0-9]{2}\\.dat$", all_files,
                   ignore.case = TRUE, value = TRUE)
data_years <- sort(2000 + as.numeric(substr(dat_files, 7, 8)))

if (length(data_years) == 0)
  stop("No NIS data files were found in ", nis_dir, ". Each file should be ",
       "named NISPUFyy.DAT (for example NISPUF24.DAT) and sit directly in ",
       "that folder, not in a subfolder.", call. = FALSE)

has_program <- sapply(data_years, function(y) {
  yy <- sprintf("%02d", y - 2000)
  !is.na(find_file(paste0("NISPUF", yy, ".R"))) ||
    !is.na(find_file(paste0("NISPUF", yy, ".sas")))
})

run_years <- data_years[has_program]

if (any(!has_program))
  cat("Skipped, no CDC program for:",
      paste(data_years[!has_program], collapse = " "), "\n")

if (length(run_years) == 0)
  stop("Data files were found, but none has its CDC program (NISPUFyy.R or ",
       "NISPUFyy.sas) beside it in ", nis_dir, ".", call. = FALSE)


# --- The field layout -------------------------------------------------------
# Each data file is fixed-width text, one child per line. CDC's program for the
# year gives the name and width of every field.

# CDC's R programs list each field as a quoted name and a width: "PDAT",1
layout_from_r <- function(pgm) {
  src <- readLines(pgm, warn = FALSE)
  a <- grep("^LIST\\.NAMEWIDTH", src)[1]
  b <- grep("^NAMEWIDTH <-", src)[1]
  if (is.na(a) || is.na(b))
    stop("Could not find the field layout in ", basename(pgm), call. = FALSE)
  pairs <- unlist(regmatches(src[a:b],
                             gregexpr('"[A-Za-z_0-9]+",[0-9]+', src[a:b])))
  list(name  = sub('^"([^"]+)",.*$', "\\1", pairs),
       width = as.numeric(sub('^.*",([0-9]+)$', "\\1", pairs)))
}

# CDC's SAS programs list each field as a starting column, a name and a
# width: @13 PROVWT 19.14
layout_from_sas <- function(pgm) {
  src <- readLines(pgm, warn = FALSE)
  ln  <- grep("^\\s*@[0-9]+\\s+[A-Za-z_][A-Za-z_0-9]*\\s+\\$?[0-9]+", src,
              value = TRUE)
  pos <- as.numeric(sub("^\\s*@([0-9]+).*$", "\\1", ln))
  nme <- sub("^\\s*@[0-9]+\\s+([A-Za-z_][A-Za-z_0-9]*).*$", "\\1", ln)
  wid <- as.numeric(sub("^\\s*@[0-9]+\\s+[A-Za-z_][A-Za-z_0-9]*\\s+\\$?([0-9]+).*$",
                        "\\1", ln))

  o <- order(pos)
  pos <- pos[o]
  nme <- nme[o]
  wid <- wid[o]

  # Each field should start where the previous one ended.
  expected <- c(1, head(pos, -1) + head(wid, -1))
  if (any(pos != expected))
    stop("SAS layout in ", basename(pgm), " is not contiguous at ",
         nme[which(pos != expected)[1]], call. = FALSE)

  list(name = nme, width = wid)
}


# --- Read one year ----------------------------------------------------------

read_year <- function(year) {
  yy  <- sprintf("%02d", year - 2000)
  dat <- find_file(paste0("NISPUF", yy, ".DAT"))
  r_p <- find_file(paste0("NISPUF", yy, ".R"))
  s_p <- find_file(paste0("NISPUF", yy, ".sas"))

  # Where a year has both programs, the R program is used.
  lay <- if (!is.na(r_p)) layout_from_r(r_p) else layout_from_sas(s_p)

  end   <- cumsum(lay$width)
  start <- end - lay$width + 1

  recs <- readLines(dat, warn = FALSE)

  # The line length the program implies must match the lines in the file.
  if (sum(lay$width) != nchar(recs[1]))
    stop("year ", year, ": program implies ", sum(lay$width),
         " characters, file has ", nchar(recs[1]), call. = FALSE)

  # The provider weight is PROVWT through 2010, PROVWT_D from 2011 to 2017
  # and PROVWT_C from 2018.
  wt <- NA_character_
  for (cand in c("PROVWT_C", "PROVWT_D", "PROVWT")) {
    if (cand %in% lay$name) { wt <- cand; break }
  }
  if (is.na(wt))
    stop("year ", year, ": no provider weight found among ",
         paste(grep("^PROVWT", lay$name, value = TRUE), collapse = " "),
         call. = FALSE)

  for (v in c("PDAT", "P_UTDPCV", "P_NUMPCV", "STATE"))
    if (!(v %in% lay$name))
      stop("year ", year, ": the layout has no field named ", v, call. = FALSE)

  # CDC marks a missing value with blanks or a dot.
  get_col <- function(v) {
    i <- match(v, lay$name)
    x <- trimws(substr(recs, start[i], end[i]))
    x[x %in% c("", ".")] <- NA
    as.numeric(x)
  }

  cat("read", year, " weight", wt, "\n")
  data.frame(PDAT     = get_col("PDAT"),
             P_UTDPCV = get_col("P_UTDPCV"),
             P_NUMPCV = get_col("P_NUMPCV"),
             STATE    = get_col("STATE"),
             weight   = get_col(wt))
}


# --- Children who count toward coverage -------------------------------------
# Children with adequate provider data (PDAT equal to 1), a value for P_UTDPCV
# and a provider weight. Territory records have no provider weight. P_NUMPCV
# is the number of pneumococcal doses. Florida is STATE 12.

fl_children <- data.frame()
us_children <- data.frame()

for (y in run_years) {
  d <- read_year(y)
  us <- d[d$PDAT %in% 1 & !is.na(d$P_UTDPCV) & !is.na(d$weight),
          c("P_NUMPCV", "STATE", "weight")]
  us$survey_year <- rep(y, nrow(us))
  us_children <- rbind(us_children, us)
  fl_children <- rbind(fl_children, us[us$STATE %in% 12, ])
}


# --- Coverage pooled over adjacent survey years -----------------------------
# Each estimate pools the survey year with the year before and the year after.
# The last estimate pools two years. The series starts at 2006.

pool <- function(children) {
  out <- data.frame()
  for (y in run_years[run_years > min(run_years)]) {
    use <- run_years[run_years >= y - 1 & run_years <= y + 1]
    s   <- children[children$survey_year %in% use, ]
    out <- rbind(out, data.frame(
      survey_year    = y,
      first_year     = min(use),
      last_year      = max(use),
      coverage_3plus = 100 * sum(s$weight[s$P_NUMPCV %in% 3:99]) / sum(s$weight),
      n_children     = nrow(s)))
  }
  out
}

write.csv(pool(fl_children), "output/nis_florida_pooled.csv", row.names = FALSE)
write.csv(pool(us_children), "output/nis_us_pooled.csv", row.names = FALSE)
