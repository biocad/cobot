module Bio.Chain.Alignment
  ( AlignmentResult (..)
  , EditDistance (..)
  , GlobalAlignment (..), LocalAlignment (..), SemiglobalAlignment (..)
  , align
  , viewAlignment
  ) where

import           Control.Lens                   ( Index
                                                , Ixed (..)
                                                , (^?!) 
                                                )
import           Data.Array                     ( array
                                                , range 
                                                )

import           Bio.Utils.Monomer              ( Symbol(..) )
import           Bio.Chain
import           Bio.Chain.Alignment.Type
import           Bio.Chain.Alignment.Algorithms

-- | Align chains using specifed algorithm
--
align :: forall algo m m'.(SequenceAlignment algo, Alignable m, Alignable m') => algo -> m -> m' -> AlignmentResult m m'
align algo s t = AlignmentResult alignmentScore alignmentResult s t
  where
    -- Bounds of chains specify bounds of alignment matrix
    (m, _M) = bounds s
    (n, _N) = bounds t
    -- Number of matricies one Match in case of simple alignment
    -- and three (Insert, Delete, Match) in case of affine
    (a, _A) = if affine algo then (Insert, Match) else (Match, Match)
    -- Bounds of alignment matrix
    bounds' :: ((Index m, Index m', EditOp), (Index m, Index m', EditOp))
    bounds' = ((m, n, a), (succ _M, succ _N, _A))
    -- Fill the matrix
    mat :: Matrix m m'
    mat = array bounds' [(ijk, dist algo mat s t ijk) | ijk <- range bounds' ]
    -- Result coordinates
    coords :: (Index m, Index m')
    coords = traceStart algo mat s t
    -- Score of alignment
    alignmentScore :: Int
    alignmentScore = let (x, y) = coords in mat ! (x, y, Match)
    -- Traceback function
    traceback :: Index m -> Index m' -> [Operation (Index m) (Index m')] -> [Operation (Index m) (Index m')]
    traceback i j ar | isStop  (cond algo) mat s t i j = ar
                     | isVert  (cond algo) mat s t i j = traceback (pred i)       j  (DELETE (pred i):ar)
                     | isHoriz (cond algo) mat s t i j = traceback       i  (pred j) (INSERT (pred j):ar)
                     | isDiag  (cond algo) mat s t i j = traceback (pred i) (pred j) (MATCH (pred i) (pred j):ar)
                     | otherwise                       = error "Alignment traceback: you cannot be here"
    -- Resulting alignment should contain additional deletions/insertions in case of semiglobal alignment 
    alignmentResult :: [Operation (Index m) (Index m')]
    alignmentResult = let preResult = uncurry traceback coords []
                      in  if not (semi algo)
                             then preResult
                             else case last preResult of
                                    MATCH i j -> preResult ++ if i /= _M then DELETE <$> [succ i .. _M]
                                                                         else INSERT <$> [succ j .. _N]
                                    _         -> error "Semiglobal alignment should always end with MATCH"

-- | View alignment results as simple strings with gaps
--
viewAlignment :: forall m m'.(Alignable m, Alignable m') => AlignmentResult m m' -> (String, String)
viewAlignment ar = unzip (toChars <$> alignment ar) 
  where
    (s, t) = (sequence1 ar, sequence2 ar)

    toChars :: Operation (Index m) (Index m') -> (Char, Char)
    toChars (MATCH i j) = (symbol (s ^?! ix i), symbol (t ^?! ix j))
    toChars (DELETE i)  = (symbol (s ^?! ix i), '-')
    toChars (INSERT j)  = ('-', symbol (t ^?! ix j))