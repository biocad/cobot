# Changelog for cobot

## Unreleased changes

## [0.1.1.3] - 2020-04-30
### Fixed
- Change Hackage upload command to `--pvp-bounds upper` instead of `both` to fix plan construction
  for older GHCs. This is needed to fix haddocks building on Hackage.

## [0.1.1.2] - 2020-04-15
### Changed
- Resolver version up.
- Move from `less-wrong` to `biocad`.

## [0.1.1.1] - 2020-03-24
### Fixed
- Affine alignment matrix construction and traceback.
### Added
- `similarity'` family of functions to calculate similarity directly on `AlignnmentResult`.
### Changed
- Generalize affine gap matrix construction;
- Improve alignment documentation.

## [0.1.1.0] - 2019-06-17
### Added
- Typeclass `IsGap`.
- Now you can align two sequences using different gaps for each one of them.
### Changed
- Removed `affine` function from `SequenceAlignment` typeclass.
