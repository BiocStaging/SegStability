# Package Name Collision Check

Date: 2026-10-06

Candidate package name: `SegStability`

Independently checked official sources (2026-10-06 UTC):

- Current Bioconductor 3.23 release and 3.24 devel `PACKAGES.gz`
  indexes at `https://bioconductor.org/packages/{version}/{category}/src/contrib/PACKAGES.gz`,
  for `bioc`, `data/annotation`, `data/experiment`, and `workflows`.
  Release counts: 2384, 928, 434, 25; devel counts: 2412, 926, 433, 26.
- Complete current [3.23 catalog](https://bioconductor.org/packages/3.23/BiocViews.html)
  and [3.24 catalog](https://bioconductor.org/packages/3.24/BiocViews.html)
  responses also contain zero case-insensitive candidate occurrences
  (2,729,450 and 2,758,531 bytes), covering listed packages even when a source
  tarball is absent from an installable-package index.
- [Removed Bioconductor packages](https://bioconductor.org/about/removed-packages/).
- Historical releases linked from the official
  [release archive](https://bioconductor.org/about/release-announcements/):
  package-name JavaScript data referenced by `BiocViews.html` for 2.5–3.22
  (`https://bioconductor.org/packages/json/{version}/{category}/packages.js`).
  For 3.6, where workflows JavaScript returned 404, all four `PACKAGES.gz`
  indexes were checked instead. For 1.8–2.4, the `Software.html`,
  `AnnotationData.html`, and `ExperimentData.html` listings were checked;
  the 2.0 software `PACKAGES.gz` index supplied its full 212-package census.
  Historical 1.5–1.7 HTML directories linked from the archive return 404;
  these unavailable records are not treated as empty namespaces. The official
  removed-package record supplies historical/removed-name coverage.
- [Current CRAN package names](https://cran.r-project.org/web/packages/available_packages_by_name.html)
  (complete response: 4,192,948 bytes).
- [CRAN archived package directories](https://cran.r-project.org/src/contrib/Archive/)
  (complete response: 5,562,115 bytes).

Method: extract `Package:` names from the decompressed current indexes and
compare each lowercased name with `segstability`. For historical package-name
data, removed-package records, and CRAN listings, search the complete responses
case-insensitively for `segstability`; zero occurrences also rules out an exact
case-insensitive package-name match. Older non-UTF-8 pages were searched as
ASCII-lowercased bytes. Page shells alone were not used as namespace evidence.
No `BiocManager::install("SegStability")` attempt was used.

Result:
- No current package named `SegStability` found.
- No historical/archived package named `SegStability` found.
- Case-insensitive matches: zero across all available sources checked.

Conclusion:
`SegStability` passes the authoritative current/historical package-name collision check as of 2026-10-06.
