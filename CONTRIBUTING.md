# Contributing to OpenClade

Bug reports, test matrices and feature requests are welcome as GitHub issues.

**Reporting a problem with an analysis.** Please attach the matrix (or a minimal version of it), the settings you used (gap treatment, outgroup, intensity, seed, number of CPU cores) and what you expected, for example a tree length reported by TNT or TreeSearch.

**Code changes.** Fork the repository, make your change with a test, and open a pull request. Run `npm test` before submitting; changes to scoring or search should also pass `npm run test:data` (see README for the `DATASETS` variable).

**Support.** This is maintained by one person. I aim to respond to issues within a week.
