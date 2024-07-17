# 0.3.0.9000

- Updated default `poly` to 2 so that limited non-linearities get discovered by default,
  as advertised (thank you to the anonymous reviewer who noticed this)
- Fixed formatting issue in documentation for `sparseR`

# 0.3.0

- Fixed bug in the `effect_plot` function that was causing issues with 
  numeric `by` variables
- Add control statement for plotting residuals (thanks to Max McGrath)
- Improving test coverage
- Allow moving legend location for `effect_plot`

# 0.2.3

- Fixed bug that caused error when making predictions 
  without outcome variable in `newdata` in `predict.sparseR`
  ((issue 3)[https://github.com/petersonR/sparseR/issues/3])
- Fixed issue in recipes revdepcheck. 

# 0.2.2

- update readme
- Fix bad test in test_sparseR_surv.R
- Fix citation style for CRAN check

# 0.2.0

- Resubmitted CRAN submission
- Updated documentation
- Changes to remain compatible with recipes 1.0.0+

# 0.1.0

- Initial CRAN submission

# 0.0.2.900

- improved coverage and documentation, fixing checks 

# 0.0.1.903

- Updated documentation, added experimental message regarding RBIC
- Update citation to new paper in AStA
- prepping for CRAN
- Improve test coverage, begin cleaning up stepwise selection code
- fix bug occurring for larger n and p problems

# 0.0.1.901

- Updated documentation
- Fixed all R CMD CHECK issues
- Updated to only include publicly-available data
- Added `pkgdown` functionality and github actions

# 0.0.0.901

- Initial private alpha version
