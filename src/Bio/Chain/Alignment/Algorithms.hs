{-# LANGUAGE InstanceSigs  #-}
{-# LANGUAGE TupleSections #-}
module Bio.Chain.Alignment.Algorithms where

import           Control.Lens             (Index, IxValue)
import qualified Data.Array.Unboxed       as A (bounds, range)
import           Data.List                (maximumBy)

import           Bio.Chain                hiding ((!))
import           Bio.Chain.Alignment.Type
import           Control.Monad            (forM_)
import           Data.Ord                 (comparing)

import           Control.Monad.ST         (ST)
import           Data.Array.Base          (readArray, writeArray)
import           Data.Array.ST            (MArray (..), STUArray, newArray,
                                           runSTUArray)
import           Data.Array.Unboxed       (Ix (..), UArray, (!))


-- | Alignnment methods
--
newtype EditDistance e1 e2       = EditDistance        (e1 -> e2 -> Bool)
data GlobalAlignment a e1 e2     = GlobalAlignment     (Scoring e1 e2) a
data LocalAlignment a e1 e2      = LocalAlignment      (Scoring e1 e2) a
data SemiglobalAlignment a e1 e2 = SemiglobalAlignment (Scoring e1 e2) a

-- Common functions

-- | Lift simple substitution function to a ChainLike collection
--
{-# SPECIALISE substitute :: (Char -> Char -> Int) -> Chain Int Char -> Chain Int Char -> Int -> Int -> Int #-}
{-# INLINE substitute #-}
substitute :: (Alignable m, Alignable m') => (IxValue m -> IxValue m' -> Int) -> m -> m' -> Index m -> Index m' -> Int
substitute f s t i j = f (s `unsafeRead` (pred i)) (t `unsafeRead` (pred j))

-- | Simple substitution function for edit distance
--
substituteED :: EditDistance e1 e2 -> (e1 -> e2 -> Int)
substituteED (EditDistance genericEq) x y = if x `genericEq` y then 1 else 0

-- | Default traceback stop condition.
--
{-# SPECIALISE defStop :: Matrix (Chain Int Char) (Chain Int Char) -> Chain Int Char -> Chain Int Char -> Int -> Int -> Bool #-}
{-# INLINE defStop #-}
defStop :: (Alignable m, Alignable m') => Matrix m m' -> m -> m' -> Index m -> Index m' -> Bool
defStop _ s t i j = let (lowerS, _) = bounds s
                        (lowerT, _) = bounds t
                    in  i == lowerS && j == lowerT

-- | Traceback stop condition for the local alignment.
--
{-# SPECIALISE localStop :: Matrix (Chain Int Char) (Chain Int Char) -> Chain Int Char -> Chain Int Char -> Int -> Int -> Bool #-}
{-# INLINE localStop #-}
localStop :: (Alignable m, Alignable m') => Matrix m m' -> m -> m' -> Index m -> Index m' -> Bool
localStop m' s t i j = let (lowerS, _) = bounds s
                           (lowerT, _) = bounds t
                       in  i == lowerS || j == lowerT || m' ! (i, j, Match) == 0

horiz :: (Alignable m, Alignable m', IsGap g) => g -> Matrix m m' -> m -> m' -> Index m -> Index m' -> Bool
horiz g m s t i j = (j > lowerT) && ((i == lowerS) || (m ! (i, pred j, Match) + add == m ! (i, j, Match)))
  where
    add | isAffine g = m ! (i, pred j, Insert)
        | otherwise  = insertCostOpen g

    (lowerT, _) = bounds t
    (lowerS, _) = bounds s

vert :: (Alignable m, Alignable m', IsGap g) => g -> Matrix m m' -> m -> m' -> Index m -> Index m' -> Bool
vert g m s t i j = (i > lowerS) && ((lowerT == j) || (m ! (pred i, j, Match) + add == m ! (i, j, Match)))
  where
    add | isAffine g = m ! (pred i, j, Delete)
        | otherwise  = deleteCostOpen g

    (lowerT, _) = bounds t
    (lowerS, _) = bounds s

-- | Default condition of moving diagonally in traceback.
--
{-# SPECIALISE defDiag :: (Char -> Char -> Int) -> Matrix (Chain Int Char) (Chain Int Char) -> Chain Int Char -> Chain Int Char -> Int -> Int -> Bool #-}
{-# INLINE defDiag #-}
defDiag :: (Alignable m, Alignable m') => (IxValue m -> IxValue m' -> Int) -> Matrix m m' -> m -> m' -> Index m -> Index m' -> Bool
defDiag sub' m s t i j = let sub = substitute sub' s t
                         in  m ! (pred i, pred j, Match) + sub i j == m ! (i, j, Match)

-- | Default start condition for traceback.
--
{-# SPECIALISE defStart :: Matrix (Chain Int Char) (Chain Int Char) -> Chain Int Char -> Chain Int Char -> (Int, Int) #-}
{-# INLINE defStart #-}
defStart :: (Alignable m, Alignable m') => Matrix m m' -> m -> m' -> (Index m, Index m')
defStart m _ _ = let ((_, _, _), (upperS, upperT, _)) = A.bounds m in (upperS, upperT)

-- | Default start condition for traceback in local alignment.
--
{-# SPECIALISE localStart :: Matrix (Chain Int Char) (Chain Int Char) -> Chain Int Char -> Chain Int Char -> (Int, Int) #-}
{-# INLINE localStart #-}
localStart :: (Alignable m, Alignable m') => Matrix m m' -> m -> m' -> (Index m, Index m')
localStart m _ _ = let ((lowerS, lowerT, _), (upperS, upperT, _)) = A.bounds m
                       range' = A.range ((lowerS, lowerT, Match), (upperS, upperT, Match))
                   in  (\(a, b, _) -> (a, b)) $ maximumBy (comparing (m !)) range'

-- | Default start condition for traceback in semiglobal alignment.
--
{-# SPECIALISE semiStart :: Matrix (Chain Int Char) (Chain Int Char) -> Chain Int Char -> Chain Int Char -> (Int, Int) #-}
{-# INLINE semiStart #-}
semiStart :: (Alignable m, Alignable m') => Matrix m m' -> m -> m' -> (Index m, Index m')
semiStart m _ _ = let ((lowerS, lowerT, _), (upperS, upperT, _)) = A.bounds m
                      lastCol = (, upperT, Match) <$> [lowerS .. upperS]
                      lastRow = (upperS, , Match) <$> [lowerT .. upperT]
                  in  (\(a, b, _) -> (a, b)) $ maximumBy (comparing (m !)) $ lastCol ++ lastRow

-- Alignment algorithm instances

-------------------
  --
  --                        About affine gaps
  --
  -- There are three matrices used in all the algorithms below:
  -- 1) One stores the resulting scores for each prefix pair;
  -- 2) One stores insertion costs in the first sequence for each prefix pair;
  -- 3) One stores insertion costs in the second sequence for each prefix pair.
  --
  -- Matrices 2 and 3 are used in affine penalty calculation:
  -- Let `gapOpen` and `gapExtend` be affine gap penalties (for opening and extending a gap correspondingly).
  -- M2[i, j] = gapOpen   if neither of sequences has insertion or deletion at prefix (i, j);
  -- M2[i, j] = gapExtend if sequence1 has insertion or, equivalently, sequence2 has deletion at prefix (i, j);
  --
  -- gap penalty in the first sequence for the prefix (i, j) is gapOpen + M2[i, j]
  -- So, if there are no gaps in the sequence before, the penalty will be `gapOpen`.
  -- Otherwise it will be `gapExtend`.
  --
  -- So, M2 holds insertion costs for the first sequence, and M3 holds insertion costs for the second sequence.
  --
  -- The resulting score is the same as in plain gap penalty:
  -- the biggest one between substitution, insertion and deletion scores.
  --
-------------------

instance IsGap g => SequenceAlignment (GlobalAlignment g) where

    -- Conditions of traceback are described below
    --
    {-# INLINE cond #-}
    cond (GlobalAlignment subC gap) = Conditions defStop (defDiag subC) (vert gap) (horiz gap)

    -- Start from bottom right corner
    --
    {-# INLINE traceStart #-}
    traceStart = const defStart

    scoreMatrix :: forall m m' . (Alignable m, Alignable m')
                => GlobalAlignment g (IxValue m) (IxValue m')
                -> m
                -> m'
                -> Matrix m m'
    scoreMatrix (GlobalAlignment subC g) s t | isAffine g = uMatrixAffine
                                             | otherwise  = uMatrixSimple
      where
        uMatrixSimple :: UArray (Index m, Index m', EditOp) Int
        uMatrixSimple = runSTUArray $ do
          matrix <- newArray ((lowerS, lowerT, Match), (nilS, nilT, Match)) 0 :: ST s (STUArray s (Index m, Index m', EditOp) Int)

          forM_ [lowerS .. nilS] $ \ixS ->
            forM_ [lowerT .. nilT] $ \ixT ->

              -- Next cell = max (d_i-1,j + gap, d_i,j-1 + gap, d_i-1,j-1 + s(i,j))
              --
              if | ixS == lowerS -> writeArray matrix (ixS, ixT, Match) $ (insertCostOpen g) * index (lowerT, nilT) ixT
                 | ixT == lowerT -> writeArray matrix (ixS, ixT, Match) $ (deleteCostOpen g) * index (lowerS, nilS) ixS
                 | otherwise -> do
                   predDiag <- matrix `readArray` (pred ixS, pred ixT, Match)
                   predS    <- matrix `readArray` (pred ixS,      ixT, Match)
                   predT    <- matrix `readArray` (     ixS, pred ixT, Match)
                   writeArray matrix (ixS, ixT, Match) $ maximum [ predDiag + sub ixS ixT
                                                                 , predS + deleteCostOpen g
                                                                 , predT + insertCostOpen g
                                                                 ]
          pure matrix

        uMatrixAffine :: UArray (Index m, Index m', EditOp) Int
        uMatrixAffine = runSTUArray $ do
          matrix <- newArray ((lowerS, lowerT, Insert), (nilS, nilT, Match)) 0 :: ST s (STUArray s (Index m, Index m', EditOp) Int)
          forM_ [lowerS .. nilS] $ \ixS ->
            forM_ [lowerT .. nilT] $ \ixT ->

              -- Next cell = max (d_i-1,j + gap, d_i,j-1 + gap, d_i-1,j-1 + s(i,j))
              -- Matrices with gap costs are also filled as follows:
              -- gepMatrix[i, j] <- gapExtend if one of strings has gap at this position else gapOpen
              --
              if | ixS == lowerS && ixT == lowerT -> do
                   writeArray matrix (ixS, ixT, Match) 0
                   writeArray matrix (ixS, ixT, Insert) $ insertCostOpen g
                   writeArray matrix (ixS, ixT, Delete) $ deleteCostOpen g
                 | ixS == lowerS -> do
                   writeArray matrix (ixS, ixT, Match) $ insertCostOpen g + (insertCostExtend g) * pred (index (lowerT, nilT) ixT)
                   writeArray matrix (ixS, ixT, Insert) $ insertCostExtend g
                   writeArray matrix (ixS, ixT, Delete) $ deleteCostOpen g
                 | ixT == lowerT -> do
                   writeArray matrix (ixS, ixT, Match) $ deleteCostOpen g + (deleteCostExtend g) * pred (index (lowerS, nilS) ixS)
                   writeArray matrix (ixS, ixT, Delete) $ deleteCostExtend g
                   writeArray matrix (ixS, ixT, Insert) $ insertCostOpen g
                 | otherwise -> do
                   predDiag <- matrix `readArray` (pred ixS, pred ixT, Match)
                   predS    <- matrix `readArray` (pred ixS,      ixT, Match)
                   predT    <- matrix `readArray` (     ixS, pred ixT, Match)

                   delCost  <- matrix `readArray` (pred ixS,      ixT, Delete)
                   insCost  <- matrix `readArray` (     ixS, pred ixT, Insert)

                   let maxScore = maximum [ predDiag + sub ixS ixT
                                          , predS + delCost
                                          , predT + insCost
                                          ]

                   writeArray matrix (ixS, ixT, Delete) $ if predS + delCost == maxScore then deleteCostExtend g else deleteCostOpen g
                   writeArray matrix (ixS, ixT, Insert) $ if predT + insCost == maxScore then insertCostExtend g else insertCostOpen g
                   writeArray matrix (ixS, ixT, Match) maxScore
          pure matrix

        (lowerS, upperS) = bounds s
        (lowerT, upperT) = bounds t
        nilS = succ upperS
        nilT = succ upperT

        sub :: Index m -> Index m' -> Int
        sub = substitute subC s t

instance IsGap g => SequenceAlignment (LocalAlignment g) where

    -- Conditions of traceback are described below
    --
    {-# INLINE cond #-}
    cond (LocalAlignment subC gap) = Conditions localStop (defDiag subC) (vert gap) (horiz gap)

    -- Start from bottom right corner
    --
    {-# INLINE traceStart #-}
    traceStart = const localStart

    scoreMatrix :: forall m m' . (Alignable m, Alignable m')
                => LocalAlignment g (IxValue m) (IxValue m')
                -> m
                -> m'
                -> Matrix m m'
    scoreMatrix (LocalAlignment subC g) s t | isAffine g = uMatrixAffine
                                            | otherwise  = uMatrixSimple
      where
        uMatrixSimple :: UArray (Index m, Index m', EditOp) Int
        uMatrixSimple = runSTUArray $ do
          matrix <- newArray ((lowerS, lowerT, Match), (nilS, nilT, Match)) 0  :: ST s (STUArray s (Index m, Index m', EditOp) Int)
          forM_ [lowerS .. nilS] $ \ixS ->
            forM_ [lowerT .. nilT] $ \ixT ->

              -- Next cell = max (d_i-1,j + gap, d_i,j-1 + gap, d_i-1,j-1 + s(i,j))
              --
              if | ixS == lowerS -> writeArray matrix (ixS, ixT, Match) 0
                 | ixT == lowerT -> writeArray matrix (ixS, ixT, Match) 0
                 | otherwise -> do
                   predDiag <- matrix `readArray` (pred ixS, pred ixT, Match)
                   predS    <- matrix `readArray` (pred ixS,      ixT, Match)
                   predT    <- matrix `readArray` (     ixS, pred ixT, Match)
                   writeArray matrix (ixS, ixT, Match) $ maximum [ predDiag + sub ixS ixT
                                                                 , predS + deleteCostOpen g
                                                                 , predT + insertCostOpen g
                                                                 , 0
                                                                 ]
          pure matrix

        uMatrixAffine :: UArray (Index m, Index m', EditOp) Int
        uMatrixAffine = runSTUArray $ do
          matrix <- newArray ((lowerS, lowerT, Insert), (nilS, nilT, Match)) 0 :: ST s (STUArray s (Index m, Index m', EditOp) Int)
          forM_ [lowerS .. nilS] $ \ixS ->
            forM_ [lowerT .. nilT] $ \ixT ->

              -- Next cell = max (d_i-1,j + gap, d_i,j-1 + gap, d_i-1,j-1 + s(i,j))
              -- Matrices with gap costs are also filled as follows:
              -- gepMatrix[i, j] <- gapExtend if one of strings has gap at this position else gapOpen
              --
              if | ixS == lowerS || ixT == lowerT -> do
                   writeArray matrix (ixS, ixT, Match)  0
                   writeArray matrix (ixS, ixT, Insert) $ insertCostOpen g
                   writeArray matrix (ixS, ixT, Delete) $ deleteCostOpen g
                 | otherwise -> do
                   predDiag <- matrix `readArray` (pred ixS, pred ixT, Match)
                   predS    <- matrix `readArray` (pred ixS,      ixT, Match)
                   predT    <- matrix `readArray` (     ixS, pred ixT, Match)

                   delCost  <- matrix `readArray` (pred ixS,      ixT, Delete)
                   insCost  <- matrix `readArray` (     ixS, pred ixT, Insert)

                   let maxScore = maximum [ predDiag + sub ixS ixT
                                          , predS + delCost
                                          , predT + insCost
                                          , 0
                                          ]

                   writeArray matrix (ixS, ixT, Delete) $ if predS + delCost == maxScore then deleteCostExtend g else deleteCostOpen g
                   writeArray matrix (ixS, ixT, Insert) $ if predT + insCost == maxScore then insertCostExtend g else insertCostOpen g
                   writeArray matrix (ixS, ixT, Match) maxScore
          pure matrix

        (lowerS, upperS) = bounds s
        (lowerT, upperT) = bounds t
        nilS = succ upperS
        nilT = succ upperT

        sub :: Index m -> Index m' -> Int
        sub = substitute subC s t

instance IsGap g => SequenceAlignment (SemiglobalAlignment g) where

    -- The alignment is semiglobal, so we have to perform some additional operations
    --
    {-# INLINE semi #-}
    semi = const True

    -- Conditions of traceback are described below
    --
    {-# INLINE cond #-}
    cond (SemiglobalAlignment subC gap) = Conditions defStop (defDiag subC) (vert gap) (horiz gap)

    -- Start from bottom right corner
    --
    {-# INLINE traceStart #-}
    traceStart = const semiStart

    scoreMatrix :: forall m m' . (Alignable m, Alignable m')
                => SemiglobalAlignment g (IxValue m) (IxValue m')
                -> m
                -> m'
                -> Matrix m m'
    scoreMatrix (SemiglobalAlignment subC g) s t | isAffine g = uMatrixAffine
                                                 | otherwise  = uMatrixSimple
      where
        uMatrixSimple :: UArray (Index m, Index m', EditOp) Int
        uMatrixSimple = runSTUArray $ do
          matrix <- newArray ((lowerS, lowerT, Match), (nilS, nilT, Match)) 0  :: ST s (STUArray s (Index m, Index m', EditOp) Int)
          forM_ [lowerS .. nilS] $ \ixS ->
            forM_ [lowerT .. nilT] $ \ixT ->

              -- Next cell = max (d_i-1,j + gap, d_i,j-1 + gap, d_i-1,j-1 + s(i,j))
              --
              if | ixS == lowerS -> writeArray matrix (ixS, ixT, Match) 0
                 | ixT == lowerT -> writeArray matrix (ixS, ixT, Match) 0
                 | otherwise -> do
                   predDiag <- matrix `readArray` (pred ixS, pred ixT, Match)
                   predS    <- matrix `readArray` (pred ixS,      ixT, Match)
                   predT    <- matrix `readArray` (     ixS, pred ixT, Match)
                   writeArray matrix (ixS, ixT, Match) $ maximum [ predDiag + sub ixS ixT
                                                                 , predS + deleteCostOpen g
                                                                 , predT + insertCostOpen g
                                                                 ]
          pure matrix

        uMatrixAffine :: UArray (Index m, Index m', EditOp) Int
        uMatrixAffine = runSTUArray $ do
          matrix <- newArray ((lowerS, lowerT, Insert), (nilS, nilT, Match)) 0 :: ST s (STUArray s (Index m, Index m', EditOp) Int)
          forM_ [lowerS .. nilS] $ \ixS ->
            forM_ [lowerT .. nilT] $ \ixT ->

              -- Next cell = max (d_i-1,j + gap, d_i,j-1 + gap, d_i-1,j-1 + s(i,j))
              -- Matrices with gap costs are also filled as follows:
              -- gepMatrix[i, j] <- gapExtend if one of strings has gap at this position else gapOpen
              --
              if | ixS == lowerS || ixT == lowerT -> do
                   writeArray matrix (ixS, ixT, Match)  0
                   writeArray matrix (ixS, ixT, Insert) $ insertCostOpen g
                   writeArray matrix (ixS, ixT, Delete) $ deleteCostOpen g
                 | otherwise -> do
                   predDiag <- matrix `readArray` (pred ixS, pred ixT, Match)
                   predS    <- matrix `readArray` (pred ixS,      ixT, Match)
                   predT    <- matrix `readArray` (     ixS, pred ixT, Match)

                   delCost  <- matrix `readArray` (pred ixS,      ixT, Delete)
                   insCost  <- matrix `readArray` (     ixS, pred ixT, Insert)

                   let maxScore = maximum [ predDiag + sub ixS ixT
                                          , predS + delCost
                                          , predT + insCost
                                          ]

                   writeArray matrix (ixS, ixT, Delete) $ if predS + delCost == maxScore then deleteCostExtend g else deleteCostOpen g
                   writeArray matrix (ixS, ixT, Insert) $ if predT + insCost == maxScore then insertCostExtend g else insertCostOpen g
                   writeArray matrix (ixS, ixT, Match) maxScore
          pure matrix

        (lowerS, upperS) = bounds s
        (lowerT, upperT) = bounds t
        nilS = succ upperS
        nilT = succ upperT

        sub :: Index m -> Index m' -> Int
        sub = substitute subC s t
