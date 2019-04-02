module Bio.Chain.Alignment
  ( AlignmentResult (..), SimpleGap, AffineGap (..), Operation (..)
  , EditDistance (..)
  , GlobalAlignment (..), LocalAlignment (..), SemiglobalAlignment (..)
  , align
  , viewAlignment
  , similarityGen
  , differenceGen
  , similarity
  , difference
  ) where

import           Control.Lens                   (Index, IxValue, Ixed (..),
                                                 (^?!))
import           Data.Array                     (array, range)

import           Bio.Chain
import           Bio.Chain.Alignment.Algorithms
import           Bio.Chain.Alignment.Type
import           Bio.Utils.Geometry             (R)
import           Bio.Utils.Monomer              (Symbol (..))

-- | Align chains using specifed algorithm
--
align :: forall algo m m'.(SequenceAlignment algo, Alignable m, Alignable m') => algo (IxValue m) (IxValue m') -> m -> m' -> AlignmentResult m m'
align algo s t = AlignmentResult alignmentScore alignmentResult s t
  where
    -- Bounds of chains specify bounds of alignment matrix
    (lowerS, upperS) = bounds s
    (lowerT, upperT) = bounds t
    -- Number of matricies one Match in case of simple alignment
    -- and three (Insert, Delete, Match) in case of affine
    (lowerOp, upperOp) = if affine algo then (Insert, Match) else (Match, Match)
    -- Bounds of alignment matrix
    bounds' :: ((Index m, Index m', EditOp), (Index m, Index m', EditOp))
    bounds' = ((lowerS, lowerT, lowerOp), (succ upperS, succ upperT, upperOp))
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
    -- Complete list of edit operations which is traceback from traceStart appended with
    -- operations left
    --
    -- Function `traceback` walks over score matrix and builds up a list of operations.
    --
    -- For local alignment it starts in the global maximum of matrix, moves towards top left
    -- corner. But it might stop at any point without reaching it.
    --
    -- For semiglobal alignment the path is different. It starts in the maximum among cells in the
    -- bottom line or right column, moves towards top left corner and always reaches it.
    --
    -- For global alignment it starts in the bottom right corner, moves towards top left corner
    -- and always reaches it.
    --
    -- Each step represents some operation: insert, delete or match. So image the path doesn't
    -- start in bottom right cell or doesn't end in top left cell. For example in case of local
    -- alignment we get couple of matches in the middle of two chains: [MATCH 5 7, MATCH 6 8].
    -- This list doesn't tell us how to get the seconds chain from the first. Obviously before
    -- matching we need to delete first four symbols from chain one and insert six symbols from
    -- chain two. And of course there are symbols left after match in both chains. And with them
    -- we need to do the same. First delete left symbols from the first chain and then insert
    -- left symbols from the second chain.
    --
    -- Local                   Semiglobal              Global
    -- +------------------+    +-----------------+     +------------------+
    -- |..................|    |X               .|     |XXXX              |
    -- |..................|    |X               .|     |    XXXX          |
    -- |..................|    |X               .|     |        XX        |
    -- |...X..............|    | X              .|     |          XX      |
    -- |....X.............|    | X              .|     |            X     |
    -- |.....XXX..........|    |  X             .|     |             XX   |
    -- |........X.........|    |   X            .|     |               X  |
    -- |..................|    |    XX          .|     |                X |
    -- |..................|    |......X..........|     |                 .|
    -- +------------------+    +-----------------+     +------------------+
    alignmentResult :: [Operation (Index m) (Index m')]
    alignmentResult =
        let preResult = uncurry traceback coords []

            -- First index of FIRST chain affected by some operation in preResult or (upperS + 1).
            firstI = head . (++ [succ upperS]) . map getI $ filter (not . isInsert) preResult
            -- First index of SECOND chain affected by some operation in preResult or (upperT + 1).
            firstJ = head . (++ [succ upperT]) . map getJ $ filter (not . isDelete) preResult
            -- Deletions and insertions of symbols before first operation in preResult.
            prefix = case head (preResult ++ [MATCH (succ upperS) (succ upperT)]) of
                        MATCH i j -> map DELETE [lowerS .. pred i] ++ map INSERT [lowerT .. pred j]
                        INSERT _ -> map DELETE [lowerS .. pred firstI]
                        DELETE _ -> map INSERT [lowerT .. pred firstJ]

            -- Last index of FIRST chain affected by some operation in preResult or (lowerS - 1).
            lastI = last . (pred lowerS :) . map getI $ filter (not . isInsert) preResult
            -- Last index of SECOND chain affected by some operation in preResult or (lowerS - 1).
            lastJ = last . (pred lowerT :) . map getJ $ filter (not . isDelete) preResult
            -- Deletions and insertions of symbols after last operation in preResult
            suffix = case last (MATCH (pred lowerS) (pred lowerT) : preResult) of
                       MATCH i j -> map DELETE [succ i .. upperS] ++ map INSERT [succ j .. upperT]
                       INSERT _ -> map DELETE [succ lastI .. upperS]
                       DELETE _ -> map INSERT [succ lastJ .. upperT]

        in  -- if preResult has no matches at all
            if null (filter isMatch preResult)
                -- only make up suffix because prefix would be same
                then preResult ++ suffix
                -- otherwise make prefix too
                else prefix ++ preResult ++ suffix

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

