module Bio.Chain.Alignment
  (
    -- * Alignment function
    align
    -- ** Alignment algorithms
  , EditDistance (..)
  , GlobalAlignment (..), LocalAlignment (..), SemiglobalAlignment (..)
  , SimpleGap, SimpleGap2, AffineGap (..), AffineGap2
  , IsGap (..)
    -- ** Alignment result types
  , AlignmentResult (..),  Operation (..)
    -- * Viewing alignment results
  , viewAlignment
  , prettyAlignmment
    -- * Similarity functions
    -- $similarity
  , similarityGen'
  , similarityGen
  , differenceGen'
  , differenceGen
  , similarity'
  , similarity
  , difference'
  , difference
  ) where

import           Control.Lens                   (Index, IxValue, Ixed (..), to,
                                                 (^?!))
import           Data.Array.Unboxed             ((!))
import           Data.List                      (intercalate)
import           Data.List.Split                (chunksOf)

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
traceback algo mat s t i' j' = helper i' j' Match []
  where
    helper i j prevOp ar
      | isStop  (cond algo) mat s t i j = ar
      | otherwise =
        let (nextOp, nextI, nextJ, op) = doMove (cond algo) mat s t i j prevOp
        in helper nextI nextJ nextOp $ op:ar

{- $similarity
These are generic variants of similarity and difference functions alongside with their specialised variants.
Generic versions take the alignment algorithm used for sequence alignment,
an equality function on elements of both sequences to calculate hamming distance on aligned sequences,
and the sequences themselves.

Sample usage of generic functions:

>>> similarityGen (GlobalAlignment (\x y -> if x == ord y then 1 else 0) (AffineGap (-11) (-1))) (\x y -> x == ord y) [ord 'R'.. ord 'z'] ['a'..'z']
0.63414633

This one will calculate similarity between a list of 'Int's and a list of 'Char's.
Generic scoring function used in alignment is @\\x y -> if x == ord y then 1 else 0@.
Generic equality function used in hamming distance is @\\x y -> x == ord y@.

Specialised versions do not take the equality function as the sequences are already constrained to have 'Eq' elements.

Sample usage of specialised function is the same as before:

>>> :{
seq1 :: String
seq1 = "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKVQLERYFDYWGQGTLVTVSS"
<BLANKLINE>
seq2 :: String
seq2 = "EVQLLESGGGLVQPGGSLRLSAAASGFTFSTFSMNWVRQAPGKGLEWVSYISRTSKTIYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYVARGRFFDYWGQGTLVTVS"
<BLANKLINE>
similarity (GlobalAlignment blosum62 (AffineGap (-11) (-1))) s1 s2
:}
0.8130081

Sometimes for biological reasons gaps appearing in one of two sequences, that are being aligned,
are not physical. For that reason we might want to use different gap penalties when aligning these sequences.

Example of usage of different gaps when aligning two sequences is presented below:

>>> :{
seq1 :: String
seq1 = "AAAAAAGGGGGGGGGGGGTTTTTTTTT"
<BLANKLINE>
seq2 :: String
seq2 = "AAAAAATTTTTTTTT"
<BLANKLINE>
gapForSeq1 :: AffineGap
gapForSeq1 = AffineGap (-5) (-1)
<BLANKLINE>
gapForSeq2 :: AffineGap
gapForSeq2 = AffineGap (-1000) (-1000) -- basically, we forbid gaps on @seq2@
<BLANKLINE>
local = LocalAlignment nuc44 (gapForSeq1, gapForSeq2)
<BLANKLINE>
viewAlignment (align local seq1 seq2)
:}
("TTTTTTTTT", "TTTTTTTTT")
-}

-- | Calculate similarity and difference between two sequences, aligning them first using given algorithm.
--
similarityGen :: forall algo m m'.(SequenceAlignment algo, Alignable m, Alignable m')
              => algo (IxValue m) (IxValue m')
              -> (IxValue m -> IxValue m' -> Bool)
              -> m
              -> m'
              -> R
similarityGen algo genericEq s t = similarityGen' (align algo s t) genericEq

-- | Calculate similarity by precomputed 'AlignmentResult'.
similarityGen' :: forall m m'. (Alignable m, Alignable m')
               => AlignmentResult m m'
               -> (IxValue m -> IxValue m' -> Bool)
               -> R
similarityGen' res genericEq = fromIntegral hamming / fromIntegral len
  where
    operations = alignment res
    s          = sequence1 res
    t          = sequence2 res
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

similarity' :: forall m m'.(Alignable m, Alignable m', IxValue m ~ IxValue m', Eq (IxValue m), Eq (IxValue m'))
            => AlignmentResult m m'
            -> R
similarity' res = similarityGen' res (==)

differenceGen :: forall algo m m'.(SequenceAlignment algo, Alignable m, Alignable m')
              => algo (IxValue m) (IxValue m')
              -> (IxValue m -> IxValue m' -> Bool)
              -> m
              -> m'
              -> R
differenceGen algo genericEq s t = 1.0 - similarityGen algo genericEq s t

differenceGen' :: forall m m'.(Alignable m, Alignable m')
               => AlignmentResult m m'
               -> (IxValue m -> IxValue m' -> Bool)
               -> R
differenceGen' res genericEq = 1.0 - similarityGen' res genericEq

difference :: forall algo m m'.(SequenceAlignment algo, Alignable m, Alignable m', IxValue m ~ IxValue m', Eq (IxValue m), Eq (IxValue m'))
           => algo (IxValue m) (IxValue m')
           -> m
           -> m'
           -> R
difference algo = differenceGen algo (==)

difference' :: forall m m'.(Alignable m, Alignable m', IxValue m ~ IxValue m', Eq (IxValue m), Eq (IxValue m'))
            => AlignmentResult m m'
            -> R
difference' res = differenceGen' res (==)

-- | View alignment results as simple strings with gaps.
--
viewAlignment :: forall m m'.(Alignable m, Alignable m', Symbol (IxValue m), Symbol (IxValue m')) => AlignmentResult m m' -> (String, String)
viewAlignment ar = unzip (toChars <$> alignment ar)
  where
    (s, t) = (sequence1 ar, sequence2 ar)

    toChars :: Operation (Index m) (Index m') -> (Char, Char)
    toChars (MATCH i j) = (symbol (s ^?! ix i), symbol (t ^?! ix j))
    toChars (DELETE i)  = (symbol (s ^?! ix i), '-')
    toChars (INSERT j)  = ('-', symbol (t ^?! ix j))

-- | Format alignment result as pretty columns of symbols.
--
-- Example with width equal to 20:
--
-- @
--  0 --------------------  0
--
--  0 TTTTTTTTTTTTTTTTTTTT 19
--
--  0 --GCCTGAATGGTGTGGTGT 17
--      || |||||| |||| |||
-- 20 TTGC-TGAATG-TGTG-TGT 36
--
-- 18 TCGGCGGAGGGACCCAGCTA 37
--     || |||||||||||||||
-- 37 -CG-CGGAGGGACCCAGCT- 53
--
-- 38 AAAAAAAAAA 47
--
-- 53 ---------- 53
-- @
prettyAlignmment
  :: forall m m'
  . (Alignable m, Alignable m', Symbol (IxValue m), Symbol (IxValue m'))
  => AlignmentResult m m' -- ^ Result of alignment to format
  -> Int                  -- ^ Desired width of one alignment row
  -> String
prettyAlignmment ar width =
  -- Due to construction 'resultRows' first element will be empty string, and we don't need it.
  intercalate "\n" $ tail resultRows
  where
    (s, t) = (sequence1 ar, sequence2 ar)
    rows = chunksOf width $ alignment ar

    chainLength :: forall c. ChainLike c => c -> Int
    chainLength ch = let (a, b) = bounds ch in fromEnum b - fromEnum a + 1

    -- Determine how many characters to leave for position numbers
    numWidth = length $ show $ max (chainLength s) (chainLength t)

    padLeft :: String -> String
    padLeft x = replicate (numWidth - length x) ' ' <> x

    -- Build one column of nice alignment like
    -- T
    -- |
    -- T
    toCharTriple :: Operation (Index m) (Index m') -> (Char, Char, Char)
    toCharTriple (MATCH i j) = (left, if left == right then '|' else ' ', right)
      where
        left  = s ^?! ix i . to symbol
        right = t ^?! ix j . to symbol
    toCharTriple (DELETE i) = (s ^?! ix i . to symbol, ' ', '-')
    toCharTriple (INSERT j) = ('-', ' ', t ^?! ix j . to symbol)

    -- Format one chunk of alignment, adding indices to start and end of strings
    -- (prevI, prevJ) must be 1-based indices of last printed characters in each string.
    formatRow :: (Int, Int) -> [Operation (Index m) (Index m')] -> ((Int, Int), [String])
    formatRow (prevI, prevJ) row = ((lastI, lastJ), [resLine1, resLine2, resLine3])
      where
        (line1, line2, line3) = unzip3 $ map toCharTriple row

        -- | Folding function to count lengths of both strings in alignment row
        countChars :: (Int, Int) -> Operation (Index m) (Index m') -> (Int, Int)
        countChars (li, lj) (MATCH _ _) = (li + 1, lj + 1)
        countChars (li, lj) (DELETE _)  = (li + 1, lj)
        countChars (li, lj) (INSERT _)  = (li, lj + 1)

        (lengthI, lengthJ) = foldl countChars (0, 0) row

        -- Indices of first printed non-gap characters in the current row.
        -- If the row contains any non-gap characters, this is equal to index
        -- of last printed character + 1
        (firstI, firstJ) =
          ( if lengthI > 0 then prevI + 1 else prevI
          , if lengthJ > 0 then prevJ + 1 else prevJ
          )
        (lastI, lastJ) = (prevI + lengthI, prevJ + lengthJ)

        -- It's easier to do everything in 1-based indices and convert before showing
        toZeroBased :: Int -> Int
        toZeroBased 0 = 0
        toZeroBased i = i - 1

        resLine1 = padLeft (show $ toZeroBased firstI) <> " " <> line1 <> " " <> padLeft (show $ toZeroBased lastI)
        resLine2 = padLeft ""                          <> " " <> line2
        resLine3 = padLeft (show $ toZeroBased firstJ) <> " " <> line3 <> " " <> padLeft (show $ toZeroBased lastJ)

    -- Go through all chunks of operations, accummulating current offsets in both strings
    (_, resultRows) =
      foldl
        (\(off, res) ops -> let (newOff, newRes) = formatRow off ops in (newOff, res <> [""] <> newRes))
        ((0, 0), [])
        rows
