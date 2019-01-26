module Bio.Chain.Alignment.Type where

import           Data.Array                     ( Array
                                                , Ix
                                                )
import           Control.Lens                   ( Index
                                                , IxValue
                                                )
import           Bio.Chain                      ( ChainLike )
import           Bio.Utils.Monomer              ( Symbol )

-- | Scoring function, returns substitution score for a couple of elements
--
type Scoring = Char -> Char -> Int

-- | Simple gap penalty
--
type SimpleGap = Int

-- | Affine gap penalty
--
data AffineGap = AffineGap { gapOpen   :: Int
                           , gapExtend :: Int
                           }

-- | Edit operation could be insertion, deletion or match/mismatch
--
data EditOp = Insert | Delete | Match
  deriving (Show, Eq, Ord, Bounded, Enum, Ix)

-- | Operation that was performed on current step of alignment
--
data Operation i j = INSERT j
                   | DELETE i
                   | MATCH i j
  deriving (Show, Eq, Ord)

-- | Alignment matrix type
--
type Matrix m m' = Array (Index m, Index m', EditOp) Int

-- | Traceback condition type
--
type Condition m m' = Matrix m m' -> m -> m' -> Index m -> Index m' -> Bool

-- | A set of traceback conditions
--
data Conditions m m' = Conditions { isStop  :: Condition m m' -- ^ Should we stop?
                                  , isDiag  :: Condition m m' -- ^ Should we go daigonally?
                                  , isVert  :: Condition m m' -- ^ Should we go vertically?
                                  , isHoriz :: Condition m m' -- ^ Should we go horizontally?
                                  }

-- | Sequence Alignment result
--
data AlignmentResult m m' = AlignmentResult { score     :: Int                              -- ^ Resulting score of alignment
                                            , alignment :: [Operation (Index m) (Index m')] -- ^ Alignment structure
                                            , sequence1 :: m                                -- ^ First chain
                                            , sequence2 :: m'                               -- ^ Second chain
                                            }

-- | Chain, that can be used for alignment
--
type Alignable m = (ChainLike m, Ix (Index m), Symbol (IxValue m))

-- |Method of sequence alignment
--
class SequenceAlignment a where
    -- | Defines wheater the alignment is affine or not
    --
    affine :: a -> Bool
    affine = const False
    -- | Defines wheater the alignment is semiglobal or not
    --
    semi :: a -> Bool
    semi = const False
    -- | Traceback conditions of alignment
    --
    cond :: (Alignable m, Alignable m') => a -> Conditions m m'
    -- | Starting position in matrix for traceback procedure
    --
    traceStart :: (Alignable m, Alignable m') => a -> Matrix m m' -> m -> m' -> (Index m, Index m')
    -- | Distance matrix element
    --
    dist :: (Alignable m, Alignable m') => a -> Matrix m m' -> m -> m' -> (Index m, Index m', EditOp) -> Int
  