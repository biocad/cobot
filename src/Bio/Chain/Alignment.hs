module Bio.Chain.Alignment
  ( AlignmentResult (..)
  , AffineGap (..)
  , EditDistance (..)
  , GlobalAlignment (..)
  , LocalAlignment (..)
  , Operation (..)
  , SemiglobalAlignment (..)
  , SimpleGap
  , align
  , difference
  , differenceGen
  , hamming
  , isDelete
  , isInsert
  , isMatch
  , similarity
  , similarityGen
  , viewAlignment
  ) where

import           Bio.Chain
import           Bio.Chain.Alignment.Algorithms
import           Bio.Chain.Alignment.Type
import           Bio.Utils.Geometry             (R)
import           Bio.Utils.Monomer              (Symbol (..))

import           Control.Lens                   (Index, IxValue, Ixed (..), (^?!))
import           Control.Monad                  (when, forM_)
import           Control.Monad.ST               (runST)
import           Data.Array.Base                (unsafeRead, unsafeWrite, unsafeNewArray_)
import           Data.Array.ST                  (index, range, rangeSize)

hammingGeneric :: forall m m'.(Alignable m, Alignable m') => (IxValue m -> IxValue m' -> Bool) -> m -> m' -> Int
hammingGeneric genericEq chain1 chain2 = sum $ zipWith score (map snd $ assocs chain1) (map snd $ assocs chain2)
  where
    size1 = rangeSize (bounds chain1)
    size2 = rangeSize (bounds chain2)
    score a b | a `genericEq` b = 1
              | otherwise       = 0

hamming :: forall m.(Alignable m, Eq (IxValue m)) => m -> m -> Int
hamming = hammingGeneric (==)

{-# SPECIALISE align :: LocalAlignment SimpleGap Char Char -> Chain Int Char -> Chain Int Char -> AlignmentResult (Chain Int Char) (Chain Int Char) #-}
{-# SPECIALISE align :: LocalAlignment AffineGap Char Char -> Chain Int Char -> Chain Int Char -> AlignmentResult (Chain Int Char) (Chain Int Char) #-}
{-# SPECIALISE align :: SemiglobalAlignment SimpleGap Char Char -> Chain Int Char -> Chain Int Char -> AlignmentResult (Chain Int Char) (Chain Int Char) #-}
{-# SPECIALISE align :: SemiglobalAlignment AffineGap Char Char -> Chain Int Char -> Chain Int Char -> AlignmentResult (Chain Int Char) (Chain Int Char) #-}
{-# SPECIALISE align :: GlobalAlignment SimpleGap Char Char -> Chain Int Char -> Chain Int Char -> AlignmentResult (Chain Int Char) (Chain Int Char) #-}
{-# SPECIALISE align :: GlobalAlignment AffineGap Char Char -> Chain Int Char -> Chain Int Char -> AlignmentResult (Chain Int Char) (Chain Int Char) #-}
align :: forall algorithm m m'.(SequenceAlignment algorithm, Alignable m, Alignable m')
      => algorithm (IxValue m) (IxValue m') -> m -> m' -> AlignmentResult m m'
align algorithm s t = runST $ do
    let (lowerS, upperS) = bounds s
    let (lowerT, upperT) = bounds t
    let bounds' = ((pred lowerS, pred lowerT), (upperS, upperT))
    scores         <- unsafeNewArray_ bounds'
    insertionCosts <- if isAffine algorithm s t then unsafeNewArray_ bounds' else pure $ error "This shouldn't be used otherwise"
    deletionCosts  <- if isAffine algorithm s t then unsafeNewArray_ bounds' else pure $ error "This shouldn't be used otherwise"
    paths          <- unsafeNewArray_ bounds'
    forM_ (range bounds') $ \ij@(i, j) -> do
        (score, insertion, deletion, path) <-
            calculateMatrix algorithm scores insertionCosts deletionCosts paths s t i j
        unsafeWrite scores         (index bounds' ij) score
        when (isAffine algorithm s t) $ do
            unsafeWrite insertionCosts (index bounds' ij) insertion
            unsafeWrite deletionCosts  (index bounds' ij) deletion
        unsafeWrite paths          (index bounds' ij) path
    let traceback ij@(i, j) = do
            operation <- unsafeRead paths (index bounds' ij)
            case operation of
                0 -> pure (ij, [])
                1 -> ((Match i j :) <$>) <$> traceback (pred i, pred j)
                2 -> ((Delete i :) <$>) <$> traceback (pred i, j)
                3 -> ((Insert j :) <$>) <$> traceback (i, pred j)
                n -> error $ "0, 1, 2 or 3 are only allowed as operataion. But got: " <> show n
    startPoint <- calculateStartPoint algorithm scores s t
    (endPoint, operations) <- (reverse <$>) <$> traceback startPoint
    (operations', matchStart, matchEnd) <- postProcessOperations algorithm operations endPoint startPoint s t
    totalScore <- unsafeRead scores (index bounds' startPoint)
    pure $ AlignmentResult { arScore       = totalScore
                           , arOperations  = operations'
                           , arFirstChain  = s
                           , arSecondChain = t
                           , arMatchRange  = (matchStart, matchEnd)
                           }

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
similarityGen :: forall m m'.(Alignable m, Alignable m')
              => AlignmentResult m m'
              -> (IxValue m -> IxValue m' -> Bool)
              -> R
similarityGen ar genericEq = result
  where
    s          = arFirstChain ar
    t          = arSecondChain ar
    operations = arOperations ar
    hamming    = sum $ map toScore operations
    result     = fromIntegral hamming / fromIntegral (length operations)

    toScore :: Operation (Index m) (Index m') -> Int
    toScore (Match i j) = if (s ^?! ix i) `genericEq` (t ^?! ix j) then 1 else 0
    toScore _           = 0

similarity :: forall m m'.(Alignable m, Alignable m', IxValue m ~ IxValue m', Eq (IxValue m), Eq (IxValue m'))
           => AlignmentResult m m' -> R
similarity algo = similarityGen algo (==)

differenceGen :: forall m m'.(Alignable m, Alignable m')
              => AlignmentResult m m' -> (IxValue m -> IxValue m' -> Bool) -> R
differenceGen algo genericEq = 1.0 - similarityGen algo genericEq

difference :: forall m m'.(Alignable m, Alignable m', IxValue m ~ IxValue m', Eq (IxValue m), Eq (IxValue m'))
           => AlignmentResult m m' -> R
difference algo = differenceGen algo (==)

-- | View alignment results as simple strings with gaps
--
viewAlignment :: forall m m'.(Alignable m, Alignable m', Symbol (IxValue m), Symbol (IxValue m')) => AlignmentResult m m' -> (String, String)
viewAlignment ar = unzip $ map toChar (arOperations ar)
  where
    (s, t) = (arFirstChain ar, arSecondChain ar)

    toChar :: Operation (Index m) (Index m') -> (Char, Char)
    toChar (Match i j) = (symbol (s ^?! ix i), symbol (t ^?! ix j))
    toChar (Delete i)  = (symbol (s ^?! ix i), '-')
    toChar (Insert j)  = ('-', symbol (t ^?! ix j))
