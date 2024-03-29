# Changelog for cobot

## Unreleased changes

## [0.1.1.8] - 2023-03-23
### Changed
- Speedup alignment of `Chain Int Char`.

## [0.1.1.7] - 2021-10-06
### Added
- `Complementary` instance for `Char`.

## [0.1.1.6] - 2021-09-14
### Added
- `ChainLike` instance for `Vector`

## [0.1.1.5] - 2021-07-03
### Changed
- `*` -> `Type` for GHC-9.

## [0.1.1.4] - 2020-04-30
### Fixed
- Tweak `package.yaml` to make `stack-2.3.1` happy.

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
