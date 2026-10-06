# Lossless-ID laboratory — Gate D

Environment: 2026-10-06, R 4.6.1, Bioconductor 3.24, S4Vectors 0.51.10,
bit64 4.8.6, testthat 3.3.2. This records validated experimental paths, not
implemented public API support. Helpers exist only under `tests/testthat/`.

| Representation | Exact operations demonstrated | Unsupported or conditional paths |
|---|---|---|
| Integer | Equality, duplicates, matching, order/sort, subset, DataFrame extraction/alignment and RDS round trips | NA rejected; signed 32-bit range only |
| Character | Same operations; UUIDs, opaque spellings and decimal values above signed 64-bit range retained | NA/empty IDs rejected; ordering uses radix to avoid locale dependence |
| LLint | Same-class equality, base matching, base order/sort, integer subsetting, DataFrame extraction/alignment, exact decimal conversion and RDS round trips | Native duplicated() fails; use validated exact decimal keys. NA subscripts fail. Mixed-double operations are forbidden |
| integer64 | Same-class equality, S3 duplicate detection, base matching, explicit bit64::order, S3 sort/subset, DataFrame extraction/alignment, exact decimal conversion and RDS round trips | base::order misorders payloads; use bit64::order. Default expect_equal tolerates distinct IDs; use exact assertions. Preserve the class; never interpret its buffer as ordinary doubles |

Both 64-bit classes preserve adjacent values above 2^53 and near both signed
boundaries. Their usable interval is -9223372036854775807 through
9223372036854775807: -9223372036854775808 is reserved for missing values.
Overflow and the reserved minimum produce missing values and are rejected.
The prototype's 9223372036854775809/9223372036854775810 case is retained as
character IDs, not narrowed into either signed 64-bit class.

Ordinary numeric/double, logical, factor, complex, raw, unsupported classes,
missing IDs, NaN and Inf are rejected by the test-only acceptance probe.
Integer and character representations remain distinct. The probe neither repairs
IDs nor defines production universe/alignment utilities.
Matching was validated within each family only. Base matching across integer and
character families conflates `1L` and `"1"`; that unguarded mixed-family path is
unsupported. Later alignment must preserve the frozen typed identity distinction.

Evidence combines adversarial tests, guarded numeric conversion methods, exact
decimal expectations, and bit-exact `identical` checks. Serialized wire formats
are not compared as scientific identity: ALTREP can serialize identical integer
values differently. Class and bit-exact identity are checked after RDS restoration.

Source inspection:

- [S4Vectors 0.51.10 official source](https://bioconductor.org/packages/3.24/bioc/src/contrib/S4Vectors_0.51.10.tar.gz):
  `R/LLint-class.R`, `R/subsetting-utils.R`, `src/LLint_class.c`; raw-byte storage,
  integer decimal parsing/formatting, same-class C comparisons and byte subsetting.
- [bit64 4.8.6 official source](https://cran.r-project.org/src/contrib/bit64_4.8.6.tar.gz):
  `src/integer64.c`, `src/sort64.h`, sorting/hash routines and registered R methods;
  long-long interpretation of its buffer, strtoll parsing, integer formatting,
  explicit integer ordering. Its double-sized buffer is not a numeric conversion.
- [R 4.6.1 sorting source](https://svn.r-project.org/R/tags/R-4-6-1/src/main/sort.c)
  and installed `rank`, `xtfrm.default`, `mtfrm.default` methods: LLint ranks compare
  objects through their comparison methods; only positional ranks are numeric.
  Matching uses exact character conversion. integer64's `mtfrm` method likewise
  returns exact character keys; base ordering's class stripping is unsafe.

Historical identity cases come from the read-only prototype's
`scripts/144_validate_final_portability_cases.py` and its recorded
`results/final_pre_api/portability_validation.json`. R container, serialization,
boundary, dispatch and assertion tests are new package-native coverage.

Later implementation may accept LLint/integer64 only through these audited paths,
with explicit dispatch and exact assertions. Native LLint duplication and base
integer64 ordering remain unsupported. No ordinary-double fallback is permitted.
