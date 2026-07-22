# The version must be one number, not three that agree by hand.
#
# Ported from rHRG, where CITATION.cff sat two releases behind DESCRIPTION,
# and from GeoPalette, where __version__ was two behind pyproject.toml. Both
# times the bump edited the places someone remembered instead of searching
# for the old number, and one of the two fixes searched for a string the file
# had never contained, so it silently changed nothing. The CITATION is what a
# DOI resolves to, so a stale one misattributes the work.
#
# This lives in tools/ and runs from the workflow rather than from tests/,
# because CITATION.cff is .Rbuildignore'd: under R CMD check it is not in
# the tarball at all, so a test there could only skip.
#
#   Rscript tools/check_versions.R
#
# Copyright (C) 2025 Igor Pawelec. Licence: GPLv3.

fail <- character(0)

desc <- read.dcf("DESCRIPTION")
version <- unname(desc[1, "Version"])
if (!grepl("^[0-9]+\\.[0-9]+\\.[0-9]+$", version))
  stop("DESCRIPTION Version is '", version, "', which is not X.Y.Z. ",
       "Everything below compares against it, so a nonsense value here ",
       "would make the rest of this script agree with nothing.")

if (file.exists("CITATION.cff")) {
  cff <- readLines("CITATION.cff", warn = FALSE)
  hit <- grep('^version:', cff, value = TRUE)
  if (!length(hit)) {
    fail <- c(fail, "CITATION.cff has no version: field")
  } else {
    got <- gsub('^version:\\s*"?|"?\\s*$', "", hit[1])
    if (got != version)
      fail <- c(fail, sprintf(
        "CITATION.cff says %s, DESCRIPTION says %s -- the CITATION is what a DOI cites",
        got, version))
  }
}

if (file.exists("NEWS.md")) {
  news <- readLines("NEWS.md", warn = FALSE)
  if (!any(grepl(paste0("^#+\\s.*", gsub(".", "\\.", version, fixed = TRUE), "\\s*$"),
                 news)))
    fail <- c(fail, sprintf("NEWS.md has no heading for %s", version))
}

if (length(fail)) {
  for (f in fail) cat("  FAIL: ", f, "\n", sep = "")
  stop("version declarations disagree", call. = FALSE)
}
cat("version ", version, " is consistent across DESCRIPTION, CITATION.cff and NEWS.md\n",
    sep = "")
