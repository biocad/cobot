module Bio.Chain.Alignment
  ( AlignmentResult (..), SimpleGap, AffineGap (..), Operation (..)
  , EditDistance (..)
  , GlobalAlignment (..), LocalAlignment (..), SemiglobalAlignment (..)
  , align
  , viewAlignment
--  , similarityGen
--  , differenceGen
--  , similarity
--  , difference
  ) where

import           Control.Lens                   (Index, IxValue, Ixed (..),
                                                 (^?!))
import           Data.Array.ST

import           Bio.Chain
import           Bio.Chain.Alignment.Algorithms
import           Bio.Chain.Alignment.Type
import           Bio.Utils.Geometry             (R)
import           Bio.Utils.Monomer              (Symbol (..))
import           Control.Monad.ST
import           Control.Monad

-- | Align chains using specifed algorithm
--
align :: forall algo m m'.(SequenceAlignment algo, Alignable m, Alignable m') => algo (IxValue m) (IxValue m') -> m -> m' -> AlignmentResult m m'
align algo s t = runST $ do
    -- Bounds of chains specify bounds of alignment matrix
    let (lowerS, upperS) = bounds s
    let (lowerT, upperT) = bounds t
    -- Number of matricies one Match in case of simple alignment
    -- and three (Insert, Delete, Match) in case of affine
    let (lowerOp, upperOp) = if affine algo then (Insert, Match) else (Match, Match)
    -- Bounds of alignment matrix
    -- bounds' :: ((Index m, Index m', EditOp), (Index m, Index m', EditOp))
    let bounds' = ((lowerS, lowerT, lowerOp), (succ upperS, succ upperT, upperOp))
    -- Fill the matrix
    mat <- newArray_ bounds'
    forM_ (range bounds') $ \ijk -> do
        value <- dist algo mat s t ijk
        writeArray mat ijk value
    -- Result coordinates
    coords@(x, y) <- traceStart algo mat s t
    -- Score of alignment
    alignmentScore <- readArray mat (x, y, Match)
    -- Traceback function
    let traceback i j ar = do
            c1 <- isStop  (cond algo) mat s t i j
            c2 <- isVert  (cond algo) mat s t i j
            c3 <- isHoriz (cond algo) mat s t i j
            c4 <- isDiag  (cond algo) mat s t i j
            if | c1 -> pure ar
               | c2 -> traceback (pred i)       j  (DELETE (pred i):ar)
               | c3 -> traceback       i  (pred j) (INSERT (pred j):ar)
               | c4 -> traceback (pred i) (pred j) (MATCH (pred i) (pred j):ar)
               | otherwise -> error "Alignment traceback: you cannot be here"
    -- Resulting alignment should contain additional deletions/insertions in case of semiglobal
    -- alignment
    alignmentResult <- do
        preResult <- uncurry traceback coords []
        -- Last index of FIRST chain affected by some operation in preResult or (lowerS - 1).
        let lastI = last . (pred lowerS :) . map getI $ filter (not . isInsert) preResult
        -- Last index of SECOND chain affected by some operation in preResult or (lowerS - 1).
        let lastJ = last . (pred lowerT :) . map getJ $ filter (not . isDelete) preResult
        -- Deletions and insertions of symbols after last operation in preResult
        let suffix = case last (MATCH (pred lowerS) (pred lowerT) : preResult) of
                        MATCH i j -> map DELETE [succ i .. upperS] ++ map INSERT [succ j .. upperT]
                        INSERT _ -> map DELETE [succ lastI .. upperS]
                        DELETE _ -> map INSERT [succ lastJ .. upperT]
        if semi algo
            then pure $ preResult ++ suffix
            else pure preResult

    pure $ AlignmentResult alignmentScore alignmentResult s t

---------------------------------------------------------------------------------------------------------
  --
  --                          Some TIPS for using the functions below
  --
  -- These are generic variants of similarity and difference functions alongside with their specialised variants.
  -- Generic versions take the alignment algorithm used for sequence alignment,
  -- an equality function on elements of both sequences to calculate hamming distance on aligned sequences,
  -- and the sequences themselves.
  --
  -- Sample usage of generic functions:
  --
  -- > similarityGen (GlobalAlignment (\x y -> if x == ord y then 1 else 0) (AffineGap (-11) (-1))) (\x y -> x == ord y) [ord 'R'.. ord 'z'] ['a'..'z']
  -- > 0.63414633
  --
  -- This one will calculate similarity between a list if `Int`s and a list of `Char`s.
  -- Generic scoring function used in alignment is `\x y -> if x == ord y then 1 else 0`
  -- Generic equality function used in hamming distance is `\x y -> x == ord y`
  --
  --
  -- Specialised versions do not take the equality function as the sequences are already constrained to have `Eq` elements.
  --
  -- Sample usage of specialised function is the same as before:
  --
  -- > seq1 :: String
  -- > seq1 = "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKVQLERYFDYWGQGTLVTVSS"
  -- >
  -- > seq2 :: String
  -- > seq2 = "EVQLLESGGGLVQPGGSLRLSAAASGFTFSTFSMNWVRQAPGKGLEWVSYISRTSKTIYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYVARGRFFDYWGQGTLVTVS"
  -- >
  -- > similarity (GlobalAlignment blosum62 (AffineGap (-11) (-1))) s1 s2
  -- > 0.8130081
  --
---------------------------------------------------------------------------------------------------------

--          -- | Calculate similarity and difference between two sequences, aligning them first using given algorithm.
--          --
--          similarityGen :: forall algo m m'.(SequenceAlignment algo, Alignable m, Alignable m')
--                        => algo (IxValue m) (IxValue m')
--                        -> (IxValue m -> IxValue m' -> Bool)
--                        -> m
--                        -> m'
--                        -> R
--          similarityGen algo genericEq s t = fromIntegral hamming / fromIntegral len
--            where
--              operations = alignment (align algo s t)
--              len        = length operations
--              hamming    = sum $ toScores <$> operations
--
--              toScores :: Operation (Index m) (Index m') -> Int
--              toScores (MATCH i j) = if (s ^?! ix i) `genericEq` (t ^?! ix j) then 1 else 0
--              toScores _           = 0
--
--          similarity :: forall algo m m'.(SequenceAlignment algo, Alignable m, Alignable m', IxValue m ~ IxValue m', Eq (IxValue m), Eq (IxValue m'))
--                     => algo (IxValue m) (IxValue m')
--                     -> m
--                     -> m'
--                     -> R
--          similarity algo = similarityGen algo (==)
--
--
--          differenceGen :: forall algo m m'.(SequenceAlignment algo, Alignable m, Alignable m')
--                        => algo (IxValue m) (IxValue m')
--                        -> (IxValue m -> IxValue m' -> Bool)
--                        -> m
--                        -> m'
--                        -> R
--          differenceGen algo genericEq s t = 1.0 - similarityGen algo genericEq s t
--
--
--          difference :: forall algo m m'.(SequenceAlignment algo, Alignable m, Alignable m', IxValue m ~ IxValue m', Eq (IxValue m), Eq (IxValue m'))
--                     => algo (IxValue m) (IxValue m')
--                     -> m
--                     -> m'
--                     -> R
--          difference algo = differenceGen algo (==)

-- | View alignment results as simple strings with gaps
--
viewAlignment :: forall m m'.(Alignable m, Alignable m', Symbol (IxValue m), Symbol (IxValue m')) => AlignmentResult m m' -> (String, String)
viewAlignment ar = unzip (toChars <$> alignment ar)
  where
    (s, t) = (sequence1 ar, sequence2 ar)

    toChars :: Operation (Index m) (Index m') -> (Char, Char)
    toChars (MATCH i j) = (symbol (s ^?! ix i), symbol (t ^?! ix j))
    toChars (DELETE i)  = (symbol (s ^?! ix i), '-')
    toChars (INSERT j)  = ('-', symbol (t ^?! ix j))

