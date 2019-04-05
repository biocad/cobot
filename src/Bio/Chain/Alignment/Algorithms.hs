{-# LANGUAGE InstanceSigs  #-}
{-# LANGUAGE TupleSections #-}
module Bio.Chain.Alignment.Algorithms where

import           Control.Lens             (Index, IxValue, ix, (^?!))
import           Data.Array.ST            (Ix (..), readArray, getBounds, range)
import           Data.List                (maximumBy)
import           Control.Monad
import           Bio.Chain
import           Bio.Chain.Alignment.Type
import           Data.Ord                 (comparing)
import           Control.Monad.ST
import           Data.Functor


-- | Alignnment methods
--
newtype EditDistance e1 e2       = EditDistance        (e1 -> e2 -> Bool)
data GlobalAlignment a e1 e2     = GlobalAlignment     (Scoring e1 e2) a
data LocalAlignment a e1 e2      = LocalAlignment      (Scoring e1 e2) a
data SemiglobalAlignment a e1 e2 = SemiglobalAlignment (Scoring e1 e2) a

-- Common functions

-- | Lift simple substitution function to a ChainLike collection
--
substitute :: (Alignable m, Alignable m') => (IxValue m -> IxValue m' -> Int) -> m -> m' -> Index m -> Index m' -> Int
substitute f s t i j = f (s ^?! ix (pred i)) (t ^?! ix (pred j))

-- | Simple substitution function for edit distance
--
substituteED :: EditDistance e1 e2 -> (e1 -> e2 -> Int)
substituteED (EditDistance genericEq) x y = if x `genericEq` y then 1 else 0

-- | Default traceback stop condition.
--
defStop :: (Alignable m, Alignable m') => Matrix s m m' -> m -> m' -> Index m -> Index m' -> ST s Bool
defStop _ s t i j = do
    let (lowerS, _) = bounds s
    let (lowerT, _) = bounds t
    pure $ i == lowerS && j == lowerT

-- | Traceback stop condition for the local alignment.
--
localStop :: (Alignable m, Alignable m') => Matrix s m m' -> m -> m' -> Index m -> Index m' -> ST s Bool
localStop m s t i j = do
    let (lowerS, _) = bounds s
    let (lowerT, _) = bounds t
    m' <- readArray m (i, j, Match)
    pure $ i == lowerS || j == lowerT || m' == 0

-- | Default condition of moving vertically in traceback.
--
defVert :: (Alignable m, Alignable m') => Int -> Matrix s m m' -> m -> m' -> Index m -> Index m' -> ST s Bool
defVert gap m s t i j
    | i <= lowerS = pure False
    | otherwise = do
        m' <- readArray m (pred i, j, Match)
        m'' <- readArray m (i, j, Match)
        pure $ (lowerT == j) || (m' + gap == m'')
  where
    (lowerT, _) = bounds t
    (lowerS, _) = bounds s

-- | Default condition of moving vertically in traceback with affine gap penalty.
--
affVert :: (Alignable m, Alignable m') => AffineGap -> Matrix s m m' -> m -> m' -> Index m -> Index m' -> ST s Bool
affVert AffineGap{..} m s t i j
    | i <= lowerS = pure False
    | otherwise = do
        m' <- readArray m (pred i, j, Match)
        m'' <- readArray m (i, j, Match)
        insertions <- readArray m (pred i, j, Insert)
        let gap = if insertions == 0 then gapOpen else gapExtend
        pure $ (lowerT == j) || (m' + gap == m'')
  where
    (lowerT, _) = bounds t
    (lowerS, _) = bounds s

-- | Default condition of moving horizontally in traceback.
--
defHoriz :: (Alignable m, Alignable m') => Int -> Matrix s m m' -> m -> m' -> Index m -> Index m' -> ST s Bool
defHoriz gap m s t i j
    | j <= lowerT = pure False
    | otherwise = do
        m' <- readArray m (i, pred j, Match)
        m'' <- readArray m (i, j, Match)
        pure $ (i == lowerS) || (m' + gap == m'')
  where
    (lowerT, _) = bounds t
    (lowerS, _) = bounds s

-- | Default condition of moving horizontally in traceback with affine gap penalty.
--
affHoriz :: (Alignable m, Alignable m') => AffineGap -> Matrix s m m' -> m -> m' -> Index m -> Index m' -> ST s Bool
affHoriz AffineGap{..} m s t i j
    | j <= lowerT = pure False
    | otherwise = do
        m' <- readArray m (i, pred j, Match)
        m'' <- readArray m (i, j, Match)
        deletions <- readArray m (i, pred j, Delete)
        let gap = if deletions == 0 then gapOpen else gapExtend
        pure $ (i == lowerS) || (m' + gap == m'')
  where
    (lowerT, _) = bounds t
    (lowerS, _) = bounds s

-- | Default condition of moving diagonally in traceback.
--
defDiag :: (Alignable m, Alignable m') => (IxValue m -> IxValue m' -> Int) -> Matrix s m m' -> m -> m' -> Index m -> Index m' -> ST s Bool
defDiag sub' m s t i j
    | j <= lowerT = pure False
    | i <= lowerS = pure False
    | otherwise = do
        let sub = substitute sub' s t
        m' <- readArray m (pred i, pred j, Match)
        m'' <- readArray m (i, j, Match)
        pure $ m' + sub i j == m''
  where
    (lowerT, _) = bounds t
    (lowerS, _) = bounds s

-- | Default start condition for traceback.
--
defStart :: (Alignable m, Alignable m') => Matrix s m m' -> m -> m' -> ST s (Index m, Index m')
defStart m _ _ = do
    ((_, _, _), (upperS, upperT, _)) <- getBounds m
    pure $ (upperS, upperT)

-- | Default start condition for traceback in local alignment.
--
localStart :: (Alignable m, Alignable m') => Matrix s m m' -> m -> m' -> ST s (Index m, Index m')
localStart m _ _ = do
    ((lowerS, lowerT, _), (upperS, upperT, _)) <- getBounds m
    let range' = range ((lowerS, lowerT, Match), (upperS, upperT, Match))
    rangeWithValues <- forM range' $ \i -> readArray m i >>= \x -> pure (i, x)
    pure $ (\((a, b, _), _) -> (a, b)) $ maximumBy (comparing snd) rangeWithValues

-- | Default start condition for traceback in semiglobal alignment.
--
semiStart :: (Alignable m, Alignable m') => Matrix s m m' -> m -> m' -> ST s (Index m, Index m')
semiStart m _ _ = do
    ((lowerS, lowerT, _), (upperS, upperT, _)) <- getBounds m
    let lastCol = (, upperT, Match) <$> [lowerS .. upperS]
    let lastRow = (upperS, , Match) <$> [lowerT .. upperT]
    rangeWithValues <- forM (lastCol ++ lastRow) $ \i -> readArray m i >>= \x -> pure (i, x)
    pure $ (\((a, b, _), _) -> (a, b)) $ maximumBy (comparing snd) rangeWithValues

-- Alignment algorithm instances

instance SequenceAlignment EditDistance where
    -- Conditions of traceback are described below
    cond ed = Conditions defStop (defDiag (substituteED ed)) (defVert 1) (defHoriz 1)
    -- Start from bottom right corner
    traceStart = const defStart
    -- Next cell = max (d_i-1,j + 1, d_i,j-1 + 1, d_i-1,j-1 + 1 if different else 0)
    dist :: forall s m m' . (Alignable m, Alignable m')
         => EditDistance (IxValue m) (IxValue m')
         -> Matrix s m m'
         -> m
         -> m'
         -> (Index m, Index m', EditOp)
         -> ST s Int
    dist ed mat s t (i, j, k) = result
      where
        sub :: Index m -> Index m' -> Int
        sub = substitute (substituteED ed) s t

        (lowerS, upperS) = bounds s
        (lowerT, upperT) = bounds t

        result :: ST s Int
        result = if | i == lowerS -> pure $ index (lowerT, succ upperT) j
                    | j == lowerT -> pure $ index (lowerS, succ upperS) i
                    | otherwise -> do
                        a <- readArray mat (pred i, pred j, k) <&> (+ sub i j)
                        b <- readArray mat (pred i,      j, k) <&> (+ 1)
                        c <- readArray mat (i,      pred j, k) <&> (+ 1)
                        pure (minimum [a, b, c])

instance SequenceAlignment (GlobalAlignment SimpleGap) where
    -- Conditions of traceback are described below
    cond (GlobalAlignment subC gap) = Conditions defStop (defDiag subC) (defVert gap) (defHoriz gap)
    -- Start from bottom right corner
    traceStart = const defStart
    -- Next cell = max (d_i-1,j + gap, d_i,j-1 + gap, d_i-1,j-1 + s(i,j))
    dist :: forall s m m' . (Alignable m, Alignable m')
         => GlobalAlignment SimpleGap (IxValue m) (IxValue m')
         -> Matrix s m m'
         -> m
         -> m'
         -> (Index m, Index m', EditOp)
         -> ST s Int
    dist (GlobalAlignment subC gap) mat s t (i, j, k) = result
      where
        sub :: Index m -> Index m' -> Int
        sub = substitute subC s t

        (lowerS, upperS) = bounds s
        (lowerT, upperT) = bounds t

        result :: ST s Int
        result = if | i == lowerS -> pure $ gap * index (lowerT, succ upperT) j
                    | j == lowerT -> pure $ gap * index (lowerS, succ upperS) i
                    | otherwise -> do
                        a <- readArray mat (pred i, pred j, k) <&> (+ sub i j)
                        b <- readArray mat (pred i,      j, k) <&> (+ gap)
                        c <- readArray mat (i,      pred j, k) <&> (+ gap)
                        pure $ maximum [a, b, c]

instance SequenceAlignment (LocalAlignment SimpleGap) where
    -- Conditions of traceback are described below
    cond (LocalAlignment subC gap) = Conditions localStop (defDiag subC) (defVert gap) (defHoriz gap)
    -- Start from bottom right corner
    traceStart = const localStart
    -- Next cell = max (d_i-1,j + gap, d_i,j-1 + gap, d_i-1,j-1 + s(i,j))
    dist :: forall s m m' . (Alignable m, Alignable m')
         => LocalAlignment SimpleGap (IxValue m) (IxValue m')
         -> Matrix s m m'
         -> m
         -> m'
         -> (Index m, Index m', EditOp)
         -> ST s Int
    dist (LocalAlignment subC gap) mat s t (i, j, k) = result
      where
        sub :: Index m -> Index m' -> Int
        sub = substitute subC s t

        (lowerS, _) = bounds s
        (lowerT, _) = bounds t

        result :: ST s Int
        result = if | i == lowerS -> pure 0
                    | j == lowerT -> pure 0
                    | otherwise -> do
                        a <- readArray mat (pred i, pred j, k) <&> (+ sub i j)
                        b <- readArray mat (pred i,      j, k) <&> (+ gap)
                        c <- readArray mat (i,      pred j, k) <&> (+ gap)
                        pure $ maximum [a, b, c, 0]

instance SequenceAlignment (SemiglobalAlignment SimpleGap) where
    -- The alignment is semiglobal, so we have to perform some additional operations
    semi = const True
    -- This is not a affine alignment, so we don't need multiple matricies
    affine = const False
    -- Conditions of traceback are described below
    cond (SemiglobalAlignment subC gap) = Conditions defStop (defDiag subC) (defVert gap) (defHoriz gap)
    -- Start from bottom right corner
    traceStart = const semiStart
    -- Next cell = max (d_i-1,j + gap, d_i,j-1 + gap, d_i-1,j-1 + s(i,j))
    dist :: forall s m m' . (Alignable m, Alignable m')
         => SemiglobalAlignment SimpleGap (IxValue m) (IxValue m')
         -> Matrix s m m'
         -> m
         -> m'
         -> (Index m, Index m', EditOp)
         -> ST s Int
    dist (SemiglobalAlignment subC gap) mat s t (i, j, k) = result
      where
        sub :: Index m -> Index m' -> Int
        sub = substitute subC s t

        (lowerS, _) = bounds s
        (lowerT, _) = bounds t

        result :: ST s Int
        result = if | i == lowerS -> pure 0
                    | j == lowerT -> pure 0
                    | otherwise -> do
                        a <- readArray mat (pred i, pred j, k) <&> (+ sub i j)
                        b <- readArray mat (pred i,      j, k) <&> (+ gap)
                        c <- readArray mat (i,      pred j, k) <&> (+ gap)
                        pure $ maximum [a, b, c]

-------------------
  --
  --                        Affine gaps
  --
  -- There are three matrices used in all the algorithms below:
  -- 1) One stores the resulting scores for each prefix pair;
  -- 2) One stores lengths of gaps in the first sequence for each prefix pair;
  -- 3) One stores lengths of gaps in the second sequence for each prefix pair.
  --
  -- Matrices 2 and 3 are used in affine penalty calculation:
  -- gap penalty in the first sequence for the prefix (i, j) is gapOpen + M2[i, j] * gapExtend
  --
  -- The resulting score is the same as in plain gap penalty:
  -- the biggest one between substitution, insertion and deletion scores.
  --
-------------------

instance SequenceAlignment (GlobalAlignment AffineGap) where
    -- The alignment uses affine gap penalty
    affine = const True
    -- Conditions of traceback are described below
    cond (GlobalAlignment subC gap) = Conditions defStop (defDiag subC) (affVert gap) (affHoriz gap)
    -- Start from bottom right corner
    traceStart = const defStart

    -- Next cell = max (d_i-1,j + gap, d_i,j-1 + gap, d_i-1,j-1 + s(i,j))
    dist :: forall s m m' . (Alignable m, Alignable m')
         => GlobalAlignment AffineGap (IxValue m) (IxValue m')
         -> Matrix s m m'
         -> m
         -> m'
         -> (Index m, Index m', EditOp)
         -> ST s Int
    dist (GlobalAlignment subC AffineGap {..}) mat s t (i, j, k) = result
      where
        sub :: Index m -> Index m' -> Int
        sub = substitute subC s t

        gapCost :: Int -> Int
        gapCost 0 = gapOpen
        gapCost _ = gapExtend

        (lowerS, upperS) = bounds s
        (lowerT, upperT) = bounds t

        -- Replacement cost at prefixes (i, j)
        replacement :: ST s Int
        replacement = readArray mat (pred i, pred j, Match) <&> (+ sub i j)

        -- Number of insertions in the `s` sequence on prefixes (i - 1, j)
        insertions :: ST s Int
        insertions = readArray mat (pred i,      j, Insert)

        -- Number of deletions in the `s` sequence on prefixes (i, j - 1)
        deletions :: ST s Int
        deletions = readArray mat (     i, pred j, Delete)

        -- Insertion cost at prefixes (i, j)
        insertion :: ST s Int
        insertion = do
            n <- insertions
            readArray mat (pred i,      j, Match) <&> (+ gapCost n)

        -- Deletion cost at prefixes (i, j)
        deletion :: ST s Int
        deletion = do
            n <- deletions
            readArray mat (     i, pred j, Match) <&> (+ gapCost n)

        maxIxValue :: ST s Int
        maxIxValue = maximum <$> sequence [replacement, insertion, deletion]

        result :: ST s Int
        result = if | i == lowerS -> pure $ gapOpen + gapExtend * index (lowerT, succ upperT) j
                    | j == lowerT -> pure $ gapOpen + gapExtend * index (lowerS, succ upperS) i
                    | k == Insert -> do
                        maxIxValue' <- maxIxValue
                        insertion' <- insertion
                        if maxIxValue' == insertion'
                            then succ <$> insertions
                            else pure 0
                    | k == Delete -> do
                        maxIxValue' <- maxIxValue
                        deletion' <- deletion
                        if maxIxValue' == deletion'
                            then succ <$> deletions
                            else pure 0
                    | otherwise -> maxIxValue


instance SequenceAlignment (LocalAlignment AffineGap) where
    -- The alignment uses affine gap penalty
    affine = const True
    -- Conditions of traceback are described below
    cond (LocalAlignment subC gap) = Conditions localStop (defDiag subC) (affVert gap) (affHoriz gap)
    -- Start from bottom right corner
    traceStart = const localStart
    -- Next cell = max (d_i-1,j + gap, d_i,j-1 + gap, d_i-1,j-1 + s(i,j))
    dist :: forall s m m' . (Alignable m, Alignable m')
         => LocalAlignment AffineGap (IxValue m) (IxValue m')
         -> Matrix s m m'
         -> m
         -> m'
         -> (Index m, Index m', EditOp)
         -> ST s Int
    dist (LocalAlignment subC AffineGap{..}) mat s t (i, j, k) = result
      where
        sub :: Index m -> Index m' -> Int
        sub = substitute subC s t

        gapCost :: Int -> Int
        gapCost 0 = gapOpen
        gapCost _ = gapExtend

        (lowerS, _) = bounds s
        (lowerT, _) = bounds t

        replacement :: ST s Int
        replacement = readArray mat (pred i, pred j, Match) <&> (+ sub i j)

        insertions :: ST s Int
        insertions = readArray mat (pred i,      j, Insert)
        deletions :: ST s Int
        deletions = readArray mat (     i, pred j, Delete)

        insertion :: ST s Int
        insertion = do
            insertions' <- insertions
            readArray mat (pred i,      j, Match) <&> (+ gapCost insertions')
        deletion :: ST s Int
        deletion = do
            deletions' <- deletions
            readArray mat (     i, pred j, Match) <&> (+ gapCost deletions')

        maxIxValue :: ST s Int
        maxIxValue = maximum <$> sequence [replacement, insertion, deletion, pure 0]

        result :: ST s Int
        result = if | i == lowerS -> pure 0
                    | j == lowerT -> pure 0
                    | k == Insert -> do
                        maxIxValue' <- maxIxValue
                        insertion' <- insertion
                        if maxIxValue' == insertion'
                            then succ <$> insertions
                            else pure 0
                    | k == Delete -> do
                        maxIxValue' <- maxIxValue
                        deletion' <- deletion
                        if maxIxValue' == deletion'
                            then succ <$> deletions
                            else pure 0
                    | otherwise -> maxIxValue

instance SequenceAlignment (SemiglobalAlignment AffineGap) where
    -- The alignment uses affine gap penalty
    affine = const True
    -- The alignment is semiglobal, so we have to perform some additional operations
    semi = const True
    -- Conditions of traceback are described below
    cond (SemiglobalAlignment subC gap) = Conditions defStop (defDiag subC) (affVert gap) (affHoriz gap)
    -- Start from bottom right corner
    traceStart = const semiStart
    -- Next cell = max (d_i-1,j + gap, d_i,j-1 + gap, d_i-1,j-1 + s(i,j))
    dist :: forall s m m' . (Alignable m, Alignable m')
         => SemiglobalAlignment AffineGap (IxValue m) (IxValue m')
         -> Matrix s m m'
         -> m
         -> m'
         -> (Index m, Index m', EditOp)
         -> ST s Int
    dist (SemiglobalAlignment subC AffineGap{..}) mat s t (i, j, k) = result
      where
        sub :: Index m -> Index m' -> Int
        sub = substitute subC s t

        gapCost :: Int -> Int
        gapCost 0 = gapOpen
        gapCost _ = gapExtend

        (lowerS, _) = bounds s
        (lowerT, _) = bounds t

        replacement :: ST s Int
        replacement = readArray mat (pred i, pred j, Match) <&> (+ sub i j)

        insertions :: ST s Int
        insertions = readArray mat (pred i,      j, Insert)
        deletions :: ST s Int
        deletions = readArray mat (     i, pred j, Delete)

        insertion :: ST s Int
        insertion = do
            insertions' <- insertions
            readArray mat (pred i,      j, Match) <&> (+ gapCost insertions')
        deletion :: ST s Int
        deletion = do
            deletions' <- deletions
            readArray mat (     i, pred j, Match) <&> (+ gapCost deletions')

        maxIxValue :: ST s Int
        maxIxValue = maximum <$> sequence [replacement, insertion, deletion]

        result :: ST s Int
        result = if | i == lowerS -> pure 0
                    | j == lowerT -> pure 0
                    | k == Insert -> do
                        maxIxValue' <- maxIxValue
                        insertion' <- insertion
                        if maxIxValue' == insertion'
                            then succ <$> insertions
                            else pure 0
                    | k == Delete -> do
                        maxIxValue' <- maxIxValue
                        deletion' <- deletion
                        if maxIxValue' == deletion'
                            then succ <$> deletions
                            else pure 0
                    | otherwise -> maxIxValue
