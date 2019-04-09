{-# LANGUAGE InstanceSigs  #-}
{-# LANGUAGE TupleSections #-}
module Bio.Chain.Alignment.Algorithms where

import           Bio.Chain
import           Bio.Chain.Alignment.Type

import           Control.Lens             (Index, IxValue, ix, (^?!))
import           Control.Monad            (forM, forM_, when)
import           Control.Monad.ST         (ST)
import           Data.Array.Base          (unsafeRead)
import           Data.Array.ST            (Ix (..), getBounds, range)
import           Data.List                (maximumBy)
import           Data.Ord                 (comparing)
import           Data.STRef               (newSTRef, readSTRef, writeSTRef)

-- | Alignnment methods
--
newtype EditDistance e1 e2         = EditDistance        (e1 -> e2 -> Bool)
data GlobalAlignment gap e1 e2     = GlobalAlignment     (Scoring e1 e2) gap
data LocalAlignment gap e1 e2      = LocalAlignment      (Scoring e1 e2) gap
data SemiglobalAlignment gap e1 e2 = SemiglobalAlignment (Scoring e1 e2) gap

calculateStartPointGlobal
    :: (SequenceAlignment alignment, Alignable m, Alignable m')
    => alignment (IxValue m) (IxValue m')
    -> Matrix s m m' Int
    -> m
    -> m'
    -> ST s (Index m, Index m')
calculateStartPointGlobal _ scores _ _ = snd <$> getBounds scores

calculateStartPointLocal
    :: (Alignable m, Alignable m')
    => LocalAlignment gap (IxValue m) (IxValue m')
    -> Matrix s m m' Int
    -> m
    -> m'
    -> ST s (Index m, Index m')
calculateStartPointLocal _ scores _ _ = do
    bounds' <- getBounds scores
    maxValueRef <- newSTRef minBound
    maxIndexRef <- newSTRef undefined
    forM_ (range bounds') $ \ij -> do
        score <- unsafeRead scores (index bounds' ij)
        maxValue <- readSTRef maxValueRef
        when (maxValue <= score) $ do
            writeSTRef maxValueRef score
            writeSTRef maxIndexRef ij
    readSTRef maxIndexRef

calculateStartPointSemiglobal
    :: (Alignable m, Alignable m')
    => SemiglobalAlignment gap (IxValue m) (IxValue m')
    -> Matrix s m m' Int
    -> m
    -> m'
    -> ST s (Index m, Index m')
calculateStartPointSemiglobal _ scores s t = do
    let (lowerS, upperS) = bounds s
    let (lowerT, upperT) = bounds t
    let bounds' = ((pred lowerS, pred lowerT), (upperS, upperT))
    let range' = map (\i -> (i, upperT)) (range (pred lowerS, upperS))
              ++ map (\j -> (upperS, j)) (range (pred lowerT, upperT))
    positionsAndScores <-
        flip traverse range' $ \ij -> do
            score <- unsafeRead scores (index bounds' ij)
            pure (ij, score)
    pure $ fst (maximumBy (comparing snd) positionsAndScores)

postProcessOperationsDoNothing
    :: (Alignable m, Alignable m')
    => algorithm (IxValue m) (IxValue m')
    -> [Operation (Index m) (Index m')]
    -> (Index m, Index m')
    -> (Index m, Index m')
    -> m
    -> m'
    -> ST s ([Operation (Index m) (Index m')], (Index m, Index m'), (Index m, Index m'))
postProcessOperationsDoNothing _ operations startMatch endMatch _ _ =
    pure (operations, startMatch, endMatch)

postProcessOperationsSemiglobal
    :: (Alignable m, Alignable m')
    => algorithm (IxValue m) (IxValue m')
    -> [Operation (Index m) (Index m')]
    -> (Index m, Index m')
    -> (Index m, Index m')
    -> m
    -> m'
    -> ST s ([Operation (Index m) (Index m')], (Index m, Index m'), (Index m, Index m'))
postProcessOperationsSemiglobal _ operations startMatch endMatch s t = do
    let (startS, startT) = startMatch
    let (endS, endT) = endMatch
    let (lowerS, upperS) = bounds s
    let (lowerT, upperT) = bounds t
    -- List of operations might be incomplete on both ends so we add missing Deletions and
    -- Insertions.
    let match = map Delete [lowerS..startS]
                ++ map Insert [lowerT..startT]
                ++ operations
                ++ map Delete [succ endS..upperS]
                ++ map Insert [succ endT..upperT]
    pure (match, (pred lowerS, pred lowerT), (upperS, upperT))

instance SequenceAlignment EditDistance where
    isAffine _ _ _ = False
    calculateStartPoint = calculateStartPointGlobal
    postProcessOperations = postProcessOperationsDoNothing

    {-# SPECIALISE
        calculateMatrix
            :: EditDistance Char Char
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Chain Int Char
            -> Chain Int Char
            -> Int
            -> Int
            -> ST s (Int, Int, Int, Int) #-}
    {-# INLINE calculateMatrix #-}
    calculateMatrix
        :: (Alignable m, Alignable m')
        => EditDistance (IxValue m) (IxValue m')
        -> Matrix s m m' Int -- scores
        -> Matrix s m m' Int -- insertionCosts
        -> Matrix s m m' Int -- deletionCosts
        -> Matrix s m m' Int -- paths
        -> m                 -- s
        -> m'                -- t
        -> Index m           -- i
        -> Index m'          -- j
        -> ST s (Int, Int, Int, Int)
    calculateMatrix algo scores _ _ _ s t i j = do
        let EditDistance eq = algo
        let (lowerS, upperS) = bounds s
        let (lowerT, upperT) = bounds t
        let bounds' = ((pred lowerS, pred lowerT), (upperS, upperT))
        let makeResult score path = pure (score, undefined, undefined, path)

        if | j == pred lowerT -> makeResult (index (bounds s) i) 2
           | i == pred lowerS -> makeResult (index (bounds t) j) 3
           | otherwise -> do
                a <- unsafeRead scores (index bounds' (pred i, pred j))
                b <- unsafeRead scores (index bounds' (pred i, j))
                c <- unsafeRead scores (index bounds' (i, pred j))
                let a' = a + if eq (s ^?! ix i) (t ^?! ix j) then 1 else 0
                let b' = b - 1
                let c' = c - 1
                let m = maximum [a', b', c']
                if | m == a'   -> makeResult a' 1
                   | m == b'   -> makeResult b' 2
                   | otherwise -> makeResult c' 3

instance SequenceAlignment (GlobalAlignment SimpleGap) where
    isAffine _ _ _ = False
    calculateStartPoint = calculateStartPointGlobal
    postProcessOperations = postProcessOperationsDoNothing

    {-# SPECIALISE
        calculateMatrix
            :: GlobalAlignment SimpleGap Char Char
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Chain Int Char
            -> Chain Int Char
            -> Int
            -> Int
            -> ST s (Int, Int, Int, Int) #-}
    {-# INLINE calculateMatrix #-}
    calculateMatrix
        :: (Alignable m, Alignable m')
        => GlobalAlignment SimpleGap (IxValue m) (IxValue m')
        -> Matrix s m m' Int -- scores
        -> Matrix s m m' Int -- insertionCosts
        -> Matrix s m m' Int -- deletionCosts
        -> Matrix s m m' Int -- paths
        -> m                 -- s
        -> m'                -- t
        -> Index m           -- i
        -> Index m'          -- j
        -> ST s (Int, Int, Int, Int)
    calculateMatrix algo scores _ _ _ s t i j = do
        let GlobalAlignment substitutionCost gapPenalty = algo
        let (lowerS, upperS) = bounds s
        let (lowerT, upperT) = bounds t
        let bounds' = ((pred lowerS, pred lowerT), (upperS, upperT))
        let makeResult score path = pure (score, undefined, undefined, path)

        if | i == pred lowerS
           , j == pred lowerT -> makeResult 0 0

           | j == pred lowerT -> do
                a <- unsafeRead scores (index bounds' (pred i, j))
                makeResult (a + gapPenalty) 2

           | i == pred lowerS -> do
                a <- unsafeRead scores (index bounds' (i, pred j))
                makeResult (a + gapPenalty) 3

           | otherwise -> do
                a <- unsafeRead scores (index bounds' (pred i, pred j))
                b <- unsafeRead scores (index bounds' (pred i, j))
                c <- unsafeRead scores (index bounds' (i, pred j))
                let a' = a + substitutionCost (s ^?! ix i) (t ^?! ix j)
                let b' = b + gapPenalty
                let c' = c + gapPenalty
                let m = maximum [a', b', c']
                if | m == a'   -> makeResult a' 1
                   | m == b'   -> makeResult b' 2
                   | otherwise -> makeResult c' 3

instance SequenceAlignment (LocalAlignment SimpleGap) where
    isAffine _ _ _ = False
    calculateStartPoint = calculateStartPointLocal
    postProcessOperations = postProcessOperationsDoNothing

    {-# SPECIALISE
        calculateMatrix
            :: LocalAlignment SimpleGap Char Char
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Chain Int Char
            -> Chain Int Char
            -> Int
            -> Int
            -> ST s (Int, Int, Int, Int) #-}
    {-# INLINE calculateMatrix #-}
    calculateMatrix
        :: (Alignable m, Alignable m')
        => LocalAlignment SimpleGap (IxValue m) (IxValue m')
        -> Matrix s m m' Int -- scores
        -> Matrix s m m' Int -- insertionCosts
        -> Matrix s m m' Int -- deletionCosts
        -> Matrix s m m' Int -- paths
        -> m                 -- s
        -> m'                -- t
        -> Index m           -- i
        -> Index m'          -- j
        -> ST s (Int, Int, Int, Int)
    calculateMatrix algo scores _ _ _ s t i j = do
        let LocalAlignment substitutionCost gapPenalty = algo
        let (lowerS, upperS) = bounds s
        let (lowerT, upperT) = bounds t
        let bounds' = ((pred lowerS, pred lowerT), (upperS, upperT))
        let makeResult score path = pure (score, undefined, undefined, path)

        if i == pred lowerS || j == pred lowerT
            then makeResult 0 0
            else do
                a <- unsafeRead scores (index bounds' (pred i, pred j))
                b <- unsafeRead scores (index bounds' (pred i, j))
                c <- unsafeRead scores (index bounds' (i, pred j))
                let a' = a + substitutionCost (s ^?! ix i) (t ^?! ix j)
                let b' = b + gapPenalty
                let c' = c + gapPenalty
                let m = maximum [a', b', c', 0]
                if | m == 0    -> makeResult 0  0
                   | m == a'   -> makeResult a' 1
                   | m == b'   -> makeResult b' 2
                   | otherwise -> makeResult c' 3

instance SequenceAlignment (SemiglobalAlignment SimpleGap) where
    isAffine _ _ _ = False
    calculateStartPoint = calculateStartPointSemiglobal
    postProcessOperations = postProcessOperationsSemiglobal

    {-# SPECIALISE
        calculateMatrix
            :: SemiglobalAlignment SimpleGap Char Char
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Chain Int Char
            -> Chain Int Char
            -> Int
            -> Int
            -> ST s (Int, Int, Int, Int) #-}
    {-# INLINE calculateMatrix #-}
    calculateMatrix
        :: (Alignable m, Alignable m')
        => SemiglobalAlignment SimpleGap (IxValue m) (IxValue m')
        -> Matrix s m m' Int -- scores
        -> Matrix s m m' Int -- insertionCosts
        -> Matrix s m m' Int -- deletionCosts
        -> Matrix s m m' Int -- paths
        -> m
        -> m'
        -> Index m  -- i
        -> Index m' -- j
        -> ST s (Int, Int, Int, Int)
    calculateMatrix algo scores _ _ _ s t i j = do
        let SemiglobalAlignment substitutionCost gapPenalty = algo
        let (lowerS, upperS) = bounds s
        let (lowerT, upperT) = bounds t
        let bounds' = ((pred lowerS, pred lowerT), (upperS, upperT))
        let makeResult score path = pure (score, undefined, undefined, path)

        if i == pred lowerS || j == pred lowerT
            then makeResult 0 0
            else do
                a <- unsafeRead scores (index bounds' (pred i, pred j))
                b <- unsafeRead scores (index bounds' (pred i, j))
                c <- unsafeRead scores (index bounds' (i, pred j))
                let a' = a + substitutionCost (s ^?! ix i) (t ^?! ix j)
                let b' = b + gapPenalty
                let c' = c + gapPenalty
                let m = maximum [a', b', c']
                if | m == a'   -> makeResult a' 1
                   | m == b'   -> makeResult b' 2
                   | otherwise -> makeResult c' 3

instance SequenceAlignment (GlobalAlignment AffineGap) where
    isAffine _ _ _ = True
    calculateStartPoint = calculateStartPointGlobal
    postProcessOperations = postProcessOperationsDoNothing

    {-# SPECIALISE
        calculateMatrix
            :: GlobalAlignment AffineGap Char Char
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Chain Int Char
            -> Chain Int Char
            -> Int
            -> Int
            -> ST s (Int, Int, Int, Int) #-}
    {-# INLINE calculateMatrix #-}
    calculateMatrix
        :: (Alignable m, Alignable m')
        => GlobalAlignment AffineGap (IxValue m) (IxValue m')
        -> Matrix s m m' Int -- scores
        -> Matrix s m m' Int -- insertionCosts
        -> Matrix s m m' Int -- deletionCosts
        -> Matrix s m m' Int -- paths
        -> m
        -> m'
        -> Index m  -- i
        -> Index m' -- j
        -> ST s (Int, Int, Int, Int)
    calculateMatrix algo scores insertionCosts deletionCosts _ s t i j = do
        let GlobalAlignment substitutionCost gap = algo
        let AffineGap gapOpenPenalty gapExtendPenalty = gap
        let (lowerS, upperS) = bounds s
        let (lowerT, upperT) = bounds t
        let bounds' = ((pred lowerS, pred lowerT), (upperS, upperT))

        if | i == pred lowerS
           , j == pred lowerT -> pure (0, gapOpenPenalty, gapOpenPenalty, 0)

           | j == pred lowerT -> do
                a            <- unsafeRead scores (index bounds' (pred i, j))
                deletionCost <- unsafeRead deletionCosts (index bounds' (pred i, j))
                pure (a + deletionCost, gapOpenPenalty, gapExtendPenalty, 2)

           | i == pred lowerS -> do
                a             <- unsafeRead scores (index bounds' (i, pred j))
                insertionCost <- unsafeRead insertionCosts (index bounds' (i, pred j))
                pure (a + insertionCost, gapExtendPenalty, gapOpenPenalty, 3)

           | otherwise -> do
                a <- unsafeRead scores (index bounds' (pred i, pred j))
                b <- unsafeRead scores (index bounds' (pred i, j))
                c <- unsafeRead scores (index bounds' (i, pred j))
                deletionCost  <- unsafeRead deletionCosts (index bounds' (pred i, j))
                insertionCost <- unsafeRead insertionCosts (index bounds' (i, pred j))
                let a' = a + substitutionCost (s ^?! ix i) (t ^?! ix j)
                let b' = b + deletionCost
                let c' = c + insertionCost
                let m = maximum [a', b', c']
                if | m == a'   -> pure (a', gapOpenPenalty, gapOpenPenalty, 1)
                   | m == b'   -> pure (b', gapOpenPenalty, gapExtendPenalty, 2)
                   | otherwise -> pure (c', gapExtendPenalty, gapOpenPenalty, 3)

instance SequenceAlignment (LocalAlignment AffineGap) where
    isAffine _ _ _ = True
    calculateStartPoint = calculateStartPointLocal
    postProcessOperations = postProcessOperationsDoNothing

    {-# SPECIALISE
        calculateMatrix
            :: LocalAlignment AffineGap Char Char
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Chain Int Char
            -> Chain Int Char
            -> Int
            -> Int
            -> ST s (Int, Int, Int, Int) #-}
    {-# INLINE calculateMatrix #-}
    calculateMatrix
        :: (Alignable m, Alignable m')
        => LocalAlignment AffineGap (IxValue m) (IxValue m')
        -> Matrix s m m' Int -- scores
        -> Matrix s m m' Int -- insertionCosts
        -> Matrix s m m' Int -- deletionCosts
        -> Matrix s m m' Int -- paths
        -> m
        -> m'
        -> Index m  -- i
        -> Index m' -- j
        -> ST s (Int, Int, Int, Int)
    calculateMatrix algo scores insertionCosts deletionCosts _ s t i j = do
        let LocalAlignment substitutionCost gap = algo
        let AffineGap gapOpenPenalty gapExtendPenalty = gap
        let (lowerS, upperS) = bounds s
        let (lowerT, upperT) = bounds t
        let bounds' = ((pred lowerS, pred lowerT), (upperS, upperT))

        if i == pred lowerS || j == pred lowerT
            then pure (0, 0, 0, 0)
            else do
                deletionCost <- unsafeRead deletionCosts (index bounds' (pred i, j))
                insertionCost <- unsafeRead insertionCosts (index bounds' (i, pred j))
                a <- unsafeRead scores (index bounds' (pred i, pred j))
                b <- unsafeRead scores (index bounds' (pred i, j))
                c <- unsafeRead scores (index bounds' (i, pred j))
                let a' = a + substitutionCost (s ^?! ix i) (t ^?! ix j)
                let b' = b + deletionCost
                let c' = c + insertionCost
                let m = maximum [a', b', c', 0]
                if | m == 0    -> pure (0, 0, 0, 0)
                   | m == a'   -> pure (a', gapOpenPenalty, gapOpenPenalty, 1)
                   | m == b'   -> pure (b', gapOpenPenalty, gapExtendPenalty, 2)
                   | otherwise -> pure (c', gapExtendPenalty, gapOpenPenalty, 3)

instance SequenceAlignment (SemiglobalAlignment AffineGap) where
    isAffine _ _ _ = True
    calculateStartPoint = calculateStartPointSemiglobal
    postProcessOperations = postProcessOperationsSemiglobal

    {-# SPECIALISE
        calculateMatrix
            :: SemiglobalAlignment AffineGap Char Char
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Matrix s (Chain Int Char) (Chain Int Char) Int
            -> Chain Int Char
            -> Chain Int Char
            -> Int
            -> Int
            -> ST s (Int, Int, Int, Int) #-}
    {-# INLINE calculateMatrix #-}
    calculateMatrix
        :: (Alignable m, Alignable m')
        => SemiglobalAlignment AffineGap (IxValue m) (IxValue m')
        -> Matrix s m m' Int -- scores
        -> Matrix s m m' Int -- insertionCosts
        -> Matrix s m m' Int -- deletionCosts
        -> Matrix s m m' Int -- paths
        -> m
        -> m'
        -> Index m  -- i
        -> Index m' -- j
        -> ST s (Int, Int, Int, Int)
    calculateMatrix algo scores insertionCosts deletionCosts _ s t i j = do
        let SemiglobalAlignment substitutionCost gap = algo
        let AffineGap gapOpenPenalty gapExtendPenalty = gap
        let (lowerS, upperS) = bounds s
        let (lowerT, upperT) = bounds t
        let bounds' = ((pred lowerS, pred lowerT), (upperS, upperT))

        if i == pred lowerS || j == pred lowerT
            then pure (0, 0, 0, 0)
            else do
                deletionCost <- unsafeRead deletionCosts (index bounds' (pred i, j))
                insertionCost <- unsafeRead insertionCosts (index bounds' (i, pred j))
                a <- unsafeRead scores (index bounds' (pred i, pred j))
                b <- unsafeRead scores (index bounds' (pred i, j))
                c <- unsafeRead scores (index bounds' (i, pred j))
                let a' = a + substitutionCost (s ^?! ix i) (t ^?! ix j)
                let b' = b + deletionCost
                let c' = c + insertionCost
                let m = maximum [a', b', c']
                if | m == a'   -> pure (a', gapOpenPenalty, gapOpenPenalty, 1)
                   | m == b'   -> pure (b', gapOpenPenalty, gapExtendPenalty, 2)
                   | otherwise -> pure (c', gapExtendPenalty, gapOpenPenalty, 3)
