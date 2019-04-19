{-# LANGUAGE InstanceSigs  #-}
{-# LANGUAGE TupleSections #-}
module Bio.Chain.Alignment.Algorithms where

import           Control.Lens             (Index, IxValue, ix, (^?!))
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
substitute f s t i j = f (s ^?! ix (pred i)) (t ^?! ix (pred j))

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

-- | Default condition of moving vertically in traceback.
--
{-# SPECIALISE defVert :: Int -> Matrix (Chain Int Char) (Chain Int Char) -> Chain Int Char -> Chain Int Char -> Int -> Int -> Bool #-}
{-# INLINE defVert #-}
defVert :: (Alignable m, Alignable m') => Int -> Matrix m m' -> m -> m' -> Index m -> Index m' -> Bool
defVert gap m s t i j = (i > lowerS) && ((lowerT == j) || (m ! (pred i, j, Match) + gap == m ! (i, j, Match)))
  where
    (lowerT, _) = bounds t
    (lowerS, _) = bounds s

-- | Default condition of moving vertically in traceback with affine gap penalty.
--
{-# SPECIALISE affVert :: AffineGap -> Matrix (Chain Int Char) (Chain Int Char) -> Chain Int Char -> Chain Int Char -> Int -> Int -> Bool #-}
{-# INLINE affVert #-}
affVert :: (Alignable m, Alignable m') => AffineGap -> Matrix m m' -> m -> m' -> Index m -> Index m' -> Bool
affVert _ m s t i j =  (i > lowerS) && ((lowerT == j) || (m ! (pred i, j, Match) + deletions == m ! (i, j, Match)))
  where
    deletions  = m ! (pred i, j, Delete)
    (lowerT, _) = bounds t
    (lowerS, _) = bounds s

-- | Default condition of moving horizontally in traceback.
--
{-# SPECIALISE defHoriz :: Int -> Matrix (Chain Int Char) (Chain Int Char) -> Chain Int Char -> Chain Int Char -> Int -> Int -> Bool #-}
{-# INLINE defHoriz #-}
defHoriz :: (Alignable m, Alignable m') => Int -> Matrix m m' -> m -> m' -> Index m -> Index m' -> Bool
defHoriz gap m s t i j = (j > lowerT) && ((i == lowerS) || (m ! (i, pred j, Match) + gap == m ! (i, j, Match)))
  where
    (lowerT, _) = bounds t
    (lowerS, _) = bounds s

-- | Default condition of moving horizontally in traceback with affine gap penalty.
--
{-# SPECIALISE affHoriz :: AffineGap -> Matrix (Chain Int Char) (Chain Int Char) -> Chain Int Char -> Chain Int Char -> Int -> Int -> Bool #-}
{-# INLINE affHoriz #-}
affHoriz :: (Alignable m, Alignable m') => AffineGap -> Matrix m m' -> m -> m' -> Index m -> Index m' -> Bool
affHoriz _ m s t i j = (j > lowerT) && ((i == lowerS) || (m ! (i, pred j, Match) + insertions == m ! (i, j, Match)))
  where
    insertions  = m ! (i, pred j, Insert)
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

instance SequenceAlignment EditDistance where

    -- Conditions of traceback are described below
    --
    {-# INLINE cond #-}
    cond ed = Conditions defStop (defDiag (substituteED ed)) (defVert 1) (defHoriz 1)

    -- Start from bottom right corner
    --
    {-# INLINE traceStart #-}
    traceStart = const defStart

    scoreMatrix :: forall m m' . (Alignable m, Alignable m')
                => EditDistance (IxValue m) (IxValue m')
                -> m
                -> m'
                -> Matrix m m'
    scoreMatrix ed s t = uMatrix
      where
        uMatrix :: UArray (Index m, Index m', EditOp) Int
        uMatrix = runSTUArray $ do
          matrix <- newArray ((lowerS, lowerT, Match), (nilS, nilT, Match)) 0 :: ST s (STUArray s (Index m, Index m', EditOp) Int)

          forM_ [lowerS .. nilS] $ \ixS ->
            forM_ [lowerT .. nilT] $ \ixT ->

              -- Next cell = max (d_i-1,j + 1, d_i,j-1 + 1, d_i-1,j-1 + 1 if different else 0)
              --
              if | ixS == lowerS -> writeArray matrix (ixS, ixT, Match) $ index (lowerT, nilT) ixT
                 | ixT == lowerT -> writeArray matrix (ixS, ixT, Match) $ index (lowerS, nilS) ixS
                 | otherwise -> do
                   predDiag <- matrix `readArray` (pred ixS, pred ixT, Match)
                   predS    <- matrix `readArray` (pred ixS,      ixT, Match)
                   predT    <- matrix `readArray` (     ixS, pred ixT, Match)
                   writeArray matrix (ixS, ixT, Match) $ maximum [ predDiag + sub ixS ixT
                                                                 , predS + 1
                                                                 , predT + 1
                                                                 ]
          pure matrix

        (lowerS, upperS) = bounds s
        (lowerT, upperT) = bounds t
        nilS = succ upperS
        nilT = succ upperT

        sub :: Index m -> Index m' -> Int
        sub = substitute (substituteED ed) s t

instance SequenceAlignment (GlobalAlignment SimpleGap) where

    -- Conditions of traceback are described below
    --
    {-# INLINE cond #-}
    cond (GlobalAlignment subC gap) = Conditions defStop (defDiag subC) (defVert gap) (defHoriz gap)

    -- Start from bottom right corner
    --
    {-# INLINE traceStart #-}
    traceStart = const defStart

    scoreMatrix :: forall m m' . (Alignable m, Alignable m')
                => GlobalAlignment SimpleGap (IxValue m) (IxValue m')
                -> m
                -> m'
                -> Matrix m m'
    scoreMatrix (GlobalAlignment subC gap) s t = uMatrix
      where
        uMatrix :: UArray (Index m, Index m', EditOp) Int
        uMatrix = runSTUArray $ do
          matrix <- newArray ((lowerS, lowerT, Match), (nilS, nilT, Match)) 0 :: ST s (STUArray s (Index m, Index m', EditOp) Int)

          forM_ [lowerS .. nilS] $ \ixS ->
            forM_ [lowerT .. nilT] $ \ixT ->

              -- Next cell = max (d_i-1,j + gap, d_i,j-1 + gap, d_i-1,j-1 + s(i,j))
              --
              if | ixS == lowerS -> writeArray matrix (ixS, ixT, Match) $ gap * index (lowerT, nilT) ixT
                 | ixT == lowerT -> writeArray matrix (ixS, ixT, Match) $ gap * index (lowerS, nilS) ixS
                 | otherwise -> do
                   predDiag <- matrix `readArray` (pred ixS, pred ixT, Match)
                   predS    <- matrix `readArray` (pred ixS,      ixT, Match)
                   predT    <- matrix `readArray` (     ixS, pred ixT, Match)
                   writeArray matrix (ixS, ixT, Match) $ maximum [ predDiag + sub ixS ixT
                                                                 , predS + gap
                                                                 , predT + gap
                                                                 ]
          pure matrix

        (lowerS, upperS) = bounds s
        (lowerT, upperT) = bounds t
        nilS = succ upperS
        nilT = succ upperT

        sub :: Index m -> Index m' -> Int
        sub = substitute subC s t



instance SequenceAlignment (LocalAlignment SimpleGap) where

    -- Conditions of traceback are described below
    --
    {-# INLINE cond #-}
    cond (LocalAlignment subC gap) = Conditions localStop (defDiag subC) (defVert gap) (defHoriz gap)

    -- Start from bottom right corner
    --
    {-# INLINE traceStart #-}
    traceStart = const localStart

    scoreMatrix :: forall m m' . (Alignable m, Alignable m')
                => LocalAlignment SimpleGap (IxValue m) (IxValue m')
                -> m
                -> m'
                -> Matrix m m'
    scoreMatrix (LocalAlignment subC gap) s t = uMatrix
      where
        uMatrix :: UArray (Index m, Index m', EditOp) Int
        uMatrix = runSTUArray $ do
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
                                                                 , predS + gap
                                                                 , predT + gap
                                                                 , 0
                                                                 ]
          pure matrix

        (lowerS, upperS) = bounds s
        (lowerT, upperT) = bounds t
        nilS = succ upperS
        nilT = succ upperT

        sub :: Index m -> Index m' -> Int
        sub = substitute subC s t


instance SequenceAlignment (SemiglobalAlignment SimpleGap) where

    -- The alignment is semiglobal, so we have to perform some additional operations
    --
    {-# INLINE semi #-}
    semi = const True

    -- This is not a affine alignment, so we don't need multiple matricies
    --
    {-# INLINE affine #-}
    affine = const False

    -- Conditions of traceback are described below
    --
    {-# INLINE cond #-}
    cond (SemiglobalAlignment subC gap) = Conditions defStop (defDiag subC) (defVert gap) (defHoriz gap)

    -- Start from bottom right corner
    --
    {-# INLINE traceStart #-}
    traceStart = const semiStart

    scoreMatrix :: forall m m' . (Alignable m, Alignable m')
                => SemiglobalAlignment SimpleGap (IxValue m) (IxValue m')
                -> m
                -> m'
                -> Matrix m m'
    scoreMatrix (SemiglobalAlignment subC gap) s t = uMatrix
      where
        uMatrix :: UArray (Index m, Index m', EditOp) Int
        uMatrix = runSTUArray $ do
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
                                                                 , predS + gap
                                                                 , predT + gap
                                                                 ]
          pure matrix

        (lowerS, upperS) = bounds s
        (lowerT, upperT) = bounds t
        nilS = succ upperS
        nilT = succ upperT

        sub :: Index m -> Index m' -> Int
        sub = substitute subC s t

-------------------
  --
  --                        Affine gaps
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

instance SequenceAlignment (GlobalAlignment AffineGap) where

    -- The alignment uses affine gap penalty
    --
    {-# INLINE affine #-}
    affine = const True

    -- Conditions of traceback are described below
    --
    {-# INLINE cond #-}
    cond (GlobalAlignment subC gap) = Conditions defStop (defDiag subC) (affVert gap) (affHoriz gap)

    -- Start from bottom right corner
    --
    {-# INLINE traceStart #-}
    traceStart = const defStart

    scoreMatrix :: forall m m' . (Alignable m, Alignable m')
                => GlobalAlignment AffineGap (IxValue m) (IxValue m')
                -> m
                -> m'
                -> Matrix m m'
    scoreMatrix (GlobalAlignment subC AffineGap {..}) s t = uMatrix
      where
        uMatrix :: UArray (Index m, Index m', EditOp) Int
        uMatrix = runSTUArray $ do
          matrix <- newArray ((lowerS, lowerT, Insert), (nilS, nilT, Match)) 0 :: ST s (STUArray s (Index m, Index m', EditOp) Int)
          forM_ [lowerS .. nilS] $ \ixS ->
            forM_ [lowerT .. nilT] $ \ixT ->

              -- Next cell = max (d_i-1,j + gap, d_i,j-1 + gap, d_i-1,j-1 + s(i,j))
              -- Matrices with gap costs are also filled as follows:
              -- gepMatrix[i, j] <- gapExtend if one of strings has gap at this position else gapOpen
              --
              if | ixS == lowerS && ixT == lowerT -> do
                   writeArray matrix (ixS, ixT, Match) 0
                   writeArray matrix (ixS, ixT, Insert) gapOpen
                   writeArray matrix (ixS, ixT, Delete) gapOpen
                 | ixS == lowerS -> do
                   writeArray matrix (ixS, ixT, Match) $ gapOpen + gapExtend * pred (index (lowerT, nilT) ixT)
                   writeArray matrix (ixS, ixT, Insert) gapExtend
                   writeArray matrix (ixS, ixT, Delete) gapOpen
                 | ixT == lowerT -> do
                   writeArray matrix (ixS, ixT, Match) $ gapOpen + gapExtend * pred (index (lowerS, nilS) ixS)
                   writeArray matrix (ixS, ixT, Delete) gapExtend
                   writeArray matrix (ixS, ixT, Insert) gapOpen
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

                   writeArray matrix (ixS, ixT, Delete) $ if predS + delCost == maxScore then gapExtend else gapOpen
                   writeArray matrix (ixS, ixT, Insert) $ if predT + insCost == maxScore then gapExtend else gapOpen
                   writeArray matrix (ixS, ixT, Match) maxScore
          pure matrix

        (lowerS, upperS) = bounds s
        (lowerT, upperT) = bounds t
        nilS = succ upperS
        nilT = succ upperT

        sub :: Index m -> Index m' -> Int
        sub = substitute subC s t


instance SequenceAlignment (LocalAlignment AffineGap) where

    -- The alignment uses affine gap penalty
    --
    {-# INLINE affine #-}
    affine = const True

    -- Conditions of traceback are described below
    --
    {-# INLINE cond #-}
    cond (LocalAlignment subC gap) = Conditions localStop (defDiag subC) (affVert gap) (affHoriz gap)

    -- Start from bottom right corner
    --
    {-# INLINE traceStart #-}
    traceStart = const localStart

    scoreMatrix :: forall m m' . (Alignable m, Alignable m')
                => LocalAlignment AffineGap (IxValue m) (IxValue m')
                -> m
                -> m'
                -> Matrix m m'
    scoreMatrix (LocalAlignment subC AffineGap {..}) s t = uMatrix
      where
        uMatrix :: UArray (Index m, Index m', EditOp) Int
        uMatrix = runSTUArray $ do
          matrix <- newArray ((lowerS, lowerT, Insert), (nilS, nilT, Match)) 0 :: ST s (STUArray s (Index m, Index m', EditOp) Int)
          forM_ [lowerS .. nilS] $ \ixS ->
            forM_ [lowerT .. nilT] $ \ixT ->

              -- Next cell = max (d_i-1,j + gap, d_i,j-1 + gap, d_i-1,j-1 + s(i,j))
              -- Matrices with gap costs are also filled as follows:
              -- gepMatrix[i, j] <- gapExtend if one of strings has gap at this position else gapOpen
              --
              if | ixS == lowerS || ixT == lowerT -> do
                   writeArray matrix (ixS, ixT, Match)  0
                   writeArray matrix (ixS, ixT, Insert) gapOpen
                   writeArray matrix (ixS, ixT, Delete) gapOpen
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

                   writeArray matrix (ixS, ixT, Delete) $ if predS + delCost == maxScore then gapExtend else gapOpen
                   writeArray matrix (ixS, ixT, Insert) $ if predT + insCost == maxScore then gapExtend else gapOpen
                   writeArray matrix (ixS, ixT, Match) maxScore
          pure matrix

        (lowerS, upperS) = bounds s
        (lowerT, upperT) = bounds t
        nilS = succ upperS
        nilT = succ upperT

        sub :: Index m -> Index m' -> Int
        sub = substitute subC s t

instance SequenceAlignment (SemiglobalAlignment AffineGap) where

    -- The alignment uses affine gap penalty
    --
    {-# INLINE affine #-}
    affine = const True

    -- The alignment is semiglobal, so we have to perform some additional operations
    --
    {-# INLINE semi #-}
    semi = const True

    -- Conditions of traceback are described below
    --
    {-# INLINE cond #-}
    cond (SemiglobalAlignment subC gap) = Conditions defStop (defDiag subC) (affVert gap) (affHoriz gap)

    -- Start from bottom right corner
    --
    {-# INLINE traceStart #-}
    traceStart = const semiStart

    scoreMatrix :: forall m m' . (Alignable m, Alignable m')
                => SemiglobalAlignment AffineGap (IxValue m) (IxValue m')
                -> m
                -> m'
                -> Matrix m m'
    scoreMatrix (SemiglobalAlignment subC AffineGap {..}) s t = uMatrix
      where
        uMatrix :: UArray (Index m, Index m', EditOp) Int
        uMatrix = runSTUArray $ do
          matrix <- newArray ((lowerS, lowerT, Insert), (nilS, nilT, Match)) 0 :: ST s (STUArray s (Index m, Index m', EditOp) Int)
          forM_ [lowerS .. nilS] $ \ixS ->
            forM_ [lowerT .. nilT] $ \ixT ->

              -- Next cell = max (d_i-1,j + gap, d_i,j-1 + gap, d_i-1,j-1 + s(i,j))
              -- Matrices with gap costs are also filled as follows:
              -- gepMatrix[i, j] <- gapExtend if one of strings has gap at this position else gapOpen
              --
              if | ixS == lowerS || ixT == lowerT -> do
                   writeArray matrix (ixS, ixT, Match)  0
                   writeArray matrix (ixS, ixT, Insert) gapOpen
                   writeArray matrix (ixS, ixT, Delete) gapOpen
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

                   writeArray matrix (ixS, ixT, Delete) $ if predS + delCost == maxScore then gapExtend else gapOpen
                   writeArray matrix (ixS, ixT, Insert) $ if predT + insCost == maxScore then gapExtend else gapOpen
                   writeArray matrix (ixS, ixT, Match) maxScore
          pure matrix

        (lowerS, upperS) = bounds s
        (lowerT, upperT) = bounds t
        nilS = succ upperS
        nilT = succ upperT

        sub :: Index m -> Index m' -> Int
        sub = substitute subC s t
