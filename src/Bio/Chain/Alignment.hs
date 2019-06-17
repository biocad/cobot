module Bio.Chain.Alignment
  ( AlignmentResult (..), SimpleGap, SimpleGap2, AffineGap (..), AffineGap2, Operation (..)
  , EditDistance (..)
  , GlobalAlignment (..), LocalAlignment (..), SemiglobalAlignment (..)
  , IsGap (..)
  , align
  , viewAlignment
  , similarityGen
  , differenceGen
  , similarity
  , difference
  ) where

import           Control.Lens                   (Index, IxValue, Ixed (..),
                                                 (^?!))
import           Data.Array.Unboxed             ((!))

import           Bio.Chain                      hiding ((!))
import           Bio.Chain.Alignment.Algorithms
import           Bio.Chain.Alignment.Type
import           Bio.Utils.Geometry             (R)
import           Bio.Utils.Monomer              (Symbol (..))

-- | Align chains using specifed algorithm
--
{-# SPECIALISE align :: LocalAlignment SimpleGap Char Char -> Chain Int Char -> Chain Int Char -> AlignmentResult (Chain Int Char) (Chain Int Char) #-}
{-# SPECIALISE align :: LocalAlignment AffineGap Char Char -> Chain Int Char -> Chain Int Char -> AlignmentResult (Chain Int Char) (Chain Int Char) #-}
{-# SPECIALISE align :: SemiglobalAlignment SimpleGap Char Char -> Chain Int Char -> Chain Int Char -> AlignmentResult (Chain Int Char) (Chain Int Char) #-}
{-# SPECIALISE align :: SemiglobalAlignment AffineGap Char Char -> Chain Int Char -> Chain Int Char -> AlignmentResult (Chain Int Char) (Chain Int Char) #-}
{-# SPECIALISE align :: GlobalAlignment SimpleGap Char Char -> Chain Int Char -> Chain Int Char -> AlignmentResult (Chain Int Char) (Chain Int Char) #-}
{-# SPECIALISE align :: GlobalAlignment AffineGap Char Char -> Chain Int Char -> Chain Int Char -> AlignmentResult (Chain Int Char) (Chain Int Char) #-}
align :: forall algo m m'.(SequenceAlignment algo, Alignable m, Alignable m') => algo (IxValue m) (IxValue m') -> m -> m' -> AlignmentResult m m'
align algo s t = AlignmentResult alignmentScore alignmentResult s t
  where
    -- Bounds of chains specify bounds of alignment matrix
    (lowerS, upperS) = bounds s
    (lowerT, upperT) = bounds t
    -- Fill the matrix
    mat :: Matrix m m'
    mat = scoreMatrix algo s t
    -- Result coordinates
    coords :: (Index m, Index m')
    coords = traceStart algo mat s t
    -- Score of alignment
    alignmentScore :: Int
    alignmentScore = let (x, y) = coords in mat ! (x, y, Match)

    -- Resulting alignment should contain additional deletions/insertions in case of semiglobal
    -- alignment
    alignmentResult :: [Operation (Index m) (Index m')]
    alignmentResult
        | semi algo = preResult ++ suffix
        | otherwise = preResult
      where
        preResult = uncurry (traceback algo mat s t) coords
        -- Last index of FIRST chain affected by some operation in preResult or (lowerS - 1).
        lastI = last . (pred lowerS :) . map getI $ filter (not . isInsert) preResult
        -- Last index of SECOND chain affected by some operation in preResult or (lowerS - 1).
        lastJ = last . (pred lowerT :) . map getJ $ filter (not . isDelete) preResult
        -- Deletions and insertions of symbols after last operation in preResult
        suffix = case last (MATCH (pred lowerS) (pred lowerT) : preResult) of
                   MATCH i j -> map DELETE [succ i .. upperS] ++ map INSERT [succ j .. upperT]
                   INSERT _ -> map DELETE [succ lastI .. upperS]
                   DELETE _ -> map INSERT [succ lastJ .. upperT]

-- | Traceback function.
--
-- Builds traceback for alignment algorithm @algo@ in matrix @mat@, that is
-- result of alignment of sequences @s@ and @t@. Traceback will start from
-- position with coordinates (@i@, @j@) in matrix.
--
-- Traceback is represented as list of 'Operation's.
--
traceback :: (SequenceAlignment algo, Alignable m, Alignable m')
          => algo (IxValue m) (IxValue m')
          -> Matrix m m'
          -> m
          -> m'
          -> Index m
          -> Index m'
          -> [Operation (Index m) (Index m')]
traceback algo mat s t i' j' = helper i' j' []
  where
    helper i j ar | isStop  (cond algo) mat s t i j = ar
                  | isVert  (cond algo) mat s t i j = helper (pred i) j        (DELETE (pred i):ar)
                  | isHoriz (cond algo) mat s t i j = helper i        (pred j) (INSERT (pred j):ar)
                  | isDiag  (cond algo) mat s t i j = helper (pred i) (pred j) (MATCH (pred i) (pred j):ar)
                  | otherwise                       = error "Alignment traceback: you cannot be here"

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
  -- Sometimes for biological reasons gaps appearing in one of two sequences, that are being aligned,
  -- are not physical. For that reason we might want to use different gap penalties when aligning these sequences.
  --
  -- Example of usage of different gaps when aligning two sequences is presented below:
  --
  -- > seq1 :: String
  -- > seq1 = "AAAAAAGGGGGGGGGGGGTTTTTTTTT"
  -- >
  -- > seq2 :: String
  -- > seq2 = "AAAAAATTTTTTTTT"
  -- >
  -- > gapForSeq1 :: AffineGap
  -- > gapForSeq1 = AffineGap (-5) (-1)
  -- >
  -- > gapForSeq2 :: AffineGap
  -- > gapForSeq2 = AffineGap (-1000) (-1000) -- basically, we forbid gaps on @seq2@
  -- >
  -- > local = LocalAlignment nuc44 (gapForSeq1, gapForSeq2)
  -- >
  -- > viewAlignment (align local seq1 seq2) == ("TTTTTTTTT", "TTTTTTTTT")
  --
---------------------------------------------------------------------------------------------------------

-- | Calculate similarity and difference between two sequences, aligning them first using given algorithm.
--
similarityGen :: forall algo m m'.(SequenceAlignment algo, Alignable m, Alignable m')
              => algo (IxValue m) (IxValue m')
              -> (IxValue m -> IxValue m' -> Bool)
              -> m
              -> m'
              -> R
similarityGen algo genericEq s t = fromIntegral hamming / fromIntegral len
  where
    operations = alignment (align algo s t)
    len        = length operations
    hamming    = sum $ toScores <$> operations

    toScores :: Operation (Index m) (Index m') -> Int
    toScores (MATCH i j) = if (s ^?! ix i) `genericEq` (t ^?! ix j) then 1 else 0
    toScores _           = 0

similarity :: forall algo m m'.(SequenceAlignment algo, Alignable m, Alignable m', IxValue m ~ IxValue m', Eq (IxValue m), Eq (IxValue m'))
           => algo (IxValue m) (IxValue m')
           -> m
           -> m'
           -> R
similarity algo = similarityGen algo (==)


differenceGen :: forall algo m m'.(SequenceAlignment algo, Alignable m, Alignable m')
              => algo (IxValue m) (IxValue m')
              -> (IxValue m -> IxValue m' -> Bool)
              -> m
              -> m'
              -> R
differenceGen algo genericEq s t = 1.0 - similarityGen algo genericEq s t


difference :: forall algo m m'.(SequenceAlignment algo, Alignable m, Alignable m', IxValue m ~ IxValue m', Eq (IxValue m), Eq (IxValue m'))
           => algo (IxValue m) (IxValue m')
           -> m
           -> m'
           -> R
difference algo = differenceGen algo (==)

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

