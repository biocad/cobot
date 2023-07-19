{-# LANGUAGE BangPatterns #-}
module Bio.Chain.Alignment.Scoring.Loader where

-- https://github.com/haskell/core-libraries-committee/blob/main/guides/export-lifta2-prelude.md
import Prelude hiding (Applicative(..))
import Control.Applicative (Applicative(..))

class ScoringMatrix a where
    scoring :: a -> Char -> Char -> Int

loadMatrix :: String -> [((Char, Char), Int)]
loadMatrix txt = concatMap (lineMap . words) (tail txtlns)
    where
    -- Lines of matrix
    txtlns :: [String]
    !txtlns  = filter (liftA2 (||) null ((/= '#') . head)) (strip <$> lines txt)

    -- Letters of matrix
    letters :: [Char]
    !letters = (map head . words . head) txtlns

    -- Strip spaces
    strip :: String -> String
    strip = reverse . dropWhile (== ' ') . reverse . dropWhile (== ' ')

    -- Map one line to a matrix
    lineMap :: [String] -> [((Char, Char), Int)]
    lineMap []       = []
    lineMap (x : xs) =
        let !hx = head x
            g (n, c) = ((hx, c), n)
        in  g <$> (read <$> xs) `zip` letters
