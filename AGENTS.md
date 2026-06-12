# Instructions for Codex working on FLightR

This is an R package for reconstructing animal paths from solar geolocation logger data.

## Privacy rules

- Do not read, summarize, copy, move, upload, or commit files in:
  - data-raw/local_validation/private_data/
  - data-raw/local_validation/outputs/
- Do not commit files matching:
  - *.RData
  - *.rds
  - *.Rds
  - *.local.yml
  - *.local.yaml
- Real GPS and GLS/TAGS data are private and must remain local.
- Use fake or simulated data for tests.
- Public code may refer to private file paths only through example config files with placeholder paths.

## Development rules

- Work in small, reviewable changes.
- Do not modify the particle filter yet.
- First create only a local validation framework.
- Add tests using fake data only.
- Keep dependencies minimal.
- Prefer base R and existing package imports.
- Do not use `git add .`.
- Do not commit unless explicitly asked.

## R package rules

- Preserve existing exported function names unless explicitly asked.
- Use roxygen comments for new exported or user-facing functions.
- Add tests under tests/testthat/.
- Run `devtools::test()` if available.
- If a test cannot be run, explain why.