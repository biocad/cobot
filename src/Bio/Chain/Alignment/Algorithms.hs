{-# LANGUAGE TupleSections #-}
module Bio.Chain.Alignment.Algorithms where

import           Control.Lens                   ( (^?!)
                                                , Index
                                                , ix
                                                )
import           Data.List                      ( maximumBy )
import qualified Data.Array                as A ( bounds
                                                , range
                                                )
import           Data.Array                     ( Ix (..) )

import           Bio.Utils.Monomer              ( Symbol(..) )
import           Bio.Chain
import           Bio.Chain.Alignment.Type

-- | Alignnment methods
--
data EditDistance          = EditDistance
data GlobalAlignment a     = GlobalAlignment     Scoring a
data LocalAlignment a      = LocalAlignment      Scoring a
data SemiglobalAlignment a = SemiglobalAlignment Scoring a

-- Common functions

-- | Lift simple substitution function to a ChainLike collection
--
substitute :: (Alignable m, Alignable m') => (Char -> Char -> Int) -> m -> m' -> Index m -> Index m' -> Int
substitute f s t i j =
    f (symbol $ s ^?! ix (pred i)) (symbol $ t ^?! ix (pred j))

-- | Simple substitution function for edit distance
--
substituteED :: Char -> Char -> Int
substituteED c d = if c == d then 0 else 1

defStop :: (Alignable m, Alignable m') => Matrix m m' -> m -> m' -> Index m -> Index m' -> Bool
defStop _ s t i j = let (m, _) = bounds s
                        (n, _) = bounds t
                    in  i == m && j == n

localStop :: (Alignable m, Alignable m') => Matrix m m' -> m -> m' -> Index m -> Index m' -> Bool
localStop m' s t i j = let (m, _) = bounds s
                           (n, _) = bounds t
                       in  i == m || j == n || m' ! (i, j, Match) == 0

defVert :: (Alignable m, Alignable m') => Int -> Matrix m m' -> m -> m' -> Index m -> Index m' -> Bool
defVert gap m _ _ i j = m ! (pred i, j, Match) + gap == m ! (i, j, Match)

defHoriz :: (Alignable m, Alignable m') => Int -> Matrix m m' -> m -> m' -> Index m -> Index m' -> Bool
defHoriz gap m _ _ i j = m ! (i, pred j, Match) + gap == m ! (i, j, Match)

defDiag :: (Alignable m, Alignable m') => (Char -> Char -> Int) -> Matrix m m' -> m -> m' -> Index m -> Index m' -> Bool
defDiag sub' m s t i j = let sub = substitute sub' s t
                         in  m ! (pred i, pred j, Match) + sub i j == m ! (i, j, Match)

defStart :: (Alignable m, Alignable m') => Matrix m m' -> m -> m' -> (Index m, Index m')
defStart m _ _ = let ((_, _, _), (_M, _N, _)) = A.bounds m in (_M, _N)

localStart :: (Alignable m, Alignable m') => Matrix m m' -> m -> m' -> (Index m, Index m')
localStart m _ _ = (\(a, b, _) -> (a, b)) $ maximumBy (\a b -> (m ! a) `compare` (m ! b)) $ A.range (A.bounds m)

semiStart :: (Alignable m, Alignable m') => Matrix m m' -> m -> m' -> (Index m, Index m')
semiStart m _ _ = let ((_m, _n, _), (_M, _N, _)) = A.bounds m
                      lastCol = (, _N, Match) <$> [_m .. _M]
                      lastRow = (_M, , Match) <$> [_n .. _N]
                  in  (\(a, b, _) -> (a, b)) $ maximumBy (\a b -> (m ! a) `compare` (m ! b)) $ lastCol ++ lastRow

-- Alignment algorithm instances

instance SequenceAlignment EditDistance where
    -- Conditions of traceback are described below
    cond = const $ Conditions defStop (defDiag substituteED) (defVert 1) (defHoriz 1)
    -- Start from botton right corner
    traceStart = const defStart
    -- Next cell = max (d_i-1,j + 1, d_i,j-1 + 1, d_i-1,j-1 + 1 if different else 0)
    dist = const $ \mat s t  (i, j, k) -> let sub     = substitute substituteED s t
                                              (m, _M) = bounds s
                                              (n, _N) = bounds t
                                          in  if | i == m    -> index (n, succ _N) j
                                                 | j == n    -> index (m, succ _M) i
                                                 | otherwise -> minimum [ mat ! (pred i, pred j, k) + sub i j
                                                                        , mat ! (pred i,      j, k) + 1
                                                                        , mat ! (i,      pred j, k) + 1
                                                                        ]

instance SequenceAlignment (GlobalAlignment SimpleGap) where
    -- Conditions of traceback are described below
    cond (GlobalAlignment subC gap) = Conditions defStop (defDiag subC) (defVert gap) (defHoriz gap)
    -- Start from botton right corner
    traceStart = const defStart
    -- Next cell = max (d_i-1,j + gap, d_i,j-1 + gap, d_i-1,j-1 + s(i,j))
    dist (GlobalAlignment subC gap) mat s t (i, j, k) = let sub     = substitute subC s t
                                                            (m, _M) = bounds s
                                                            (n, _N) = bounds t
                                                        in  if | i == m    -> gap * index (n, succ _N) j
                                                               | j == n    -> gap * index (m, succ _M) i
                                                               | otherwise -> maximum [ mat ! (pred i, pred j, k) + sub i j
                                                                                      , mat ! (pred i,      j, k) + gap
                                                                                      , mat ! (i,      pred j, k) + gap
                                                                                      ]

instance SequenceAlignment (LocalAlignment SimpleGap) where
    -- Conditions of traceback are described below
    cond (LocalAlignment subC gap) = Conditions localStop (defDiag subC) (defVert gap) (defHoriz gap)
    -- Start from botton right corner
    traceStart = const localStart
    -- Next cell = max (d_i-1,j + gap, d_i,j-1 + gap, d_i-1,j-1 + s(i,j))
    dist (LocalAlignment subC gap) mat s t (i, j, k) = let sub     = substitute subC s t
                                                           (m, _M) = bounds s
                                                           (n, _N) = bounds t
                                                       in  if | i == m    -> 0
                                                              | j == n    -> 0
                                                              | otherwise -> maximum [ mat ! (pred i, pred j, k) + sub i j
                                                                                     , mat ! (pred i,      j, k) + gap
                                                                                     , mat ! (i,      pred j, k) + gap
                                                                                     , 0
                                                                                     ]

instance SequenceAlignment (SemiglobalAlignment SimpleGap) where
    -- The alignment is semiglobal, so we have to perform some additional operations
    semi = const True
    -- This is not a affine alignment, so we don't need multiple matricies
    affine = const False
    -- Conditions of traceback are described below
    cond (SemiglobalAlignment subC gap) = Conditions defStop (defDiag subC) (defVert gap) (defHoriz gap)
    -- Start from botton right corner
    traceStart = const semiStart
    -- Next cell = max (d_i-1,j + gap, d_i,j-1 + gap, d_i-1,j-1 + s(i,j))
    dist (SemiglobalAlignment subC gap) mat s t (i, j, k) = let sub     = substitute subC s t
                                                                (m, _M) = bounds s
                                                                (n, _N) = bounds t
                                                            in  if | i == m    -> 0
                                                                   | j == n    -> 0
                                                                   | otherwise -> maximum [ mat ! (pred i, pred j, k) + sub i j
                                                                                          , mat ! (pred i,      j, k) + gap
                                                                                          , mat ! (i,      pred j, k) + gap
                                                                                          ]