{-# LANGUAGE DeriveAnyClass     #-}
{-# LANGUAGE StandaloneDeriving #-}

module Bio.Chain.Alignment.Type where

import           Bio.Chain          (ChainLike (..))
import           Control.DeepSeq    (NFData (..))
import           Control.Lens       (Index, IxValue)
import           Data.Array.Unboxed (Ix, UArray)
import           Data.Kind          (Type)
import           GHC.Generics       (Generic (..))

-- | Scoring function, returns substitution score for a couple of elements
--
-- type Scoring = Char -> Char -> Int

type Scoring a b = a -> b -> Int

-- | Simple gap penalty
--
type SimpleGap = Int

-- | Gap penalty with different 'SimpleGap' penalties for sequences.
--
-- First element of pair is penalty for first sequence passed to alignment
-- algorithm, second element — penalty for second passed sequence.
--
type SimpleGap2 = (SimpleGap, SimpleGap)

-- | Affine gap penalty
--
data AffineGap = AffineGap { gapOpen   :: Int
                           , gapExtend :: Int
                           }
  deriving (Show, Eq, Generic, NFData)

-- | Gap penalty with different 'AffineGap' penalties for sequences.
--
-- First element of pair is penalty for first sequence passed to alignment
-- algorithm, second element — penalty for second passed sequence.
--
type AffineGap2 = (AffineGap, AffineGap)

-- | Type class that describes possible gaps in alignments.
--
class IsGap a where
  -- | Insertions are gaps in the first argument of an alignment function.
  --
  insertCostOpen :: a -> Int
  insertCostExtend :: a -> Int

  -- | Deletions are gaps in the second argument of an alignment function.
  --
  deleteCostOpen :: a -> Int
  deleteCostExtend :: a -> Int

  isAffine :: a -> Bool
  isAffine x = insertCostOpen x /= insertCostExtend x || deleteCostOpen x /= deleteCostExtend x

instance IsGap SimpleGap where
  insertCostOpen = id
  insertCostExtend = id

  deleteCostOpen = id
  deleteCostExtend = id

instance IsGap SimpleGap2 where
  insertCostOpen = fst
  insertCostExtend = fst

  deleteCostOpen = snd
  deleteCostExtend = snd

instance IsGap AffineGap where
  insertCostOpen = gapOpen
  insertCostExtend = gapExtend

  deleteCostOpen = gapOpen
  deleteCostExtend = gapExtend

instance IsGap AffineGap2 where
  insertCostOpen = gapOpen . fst
  insertCostExtend = gapExtend . fst

  deleteCostOpen = gapOpen . snd
  deleteCostExtend = gapExtend . snd

-- | Edit operation could be insertion, deletion or match/mismatch
--
data EditOp = Insert | Delete | Match
  deriving (Show, Eq, Ord, Bounded, Enum, Ix, Generic, NFData)

-- | Operation that was performed on current step of alignment
--
data Operation i j = INSERT {            getJ :: j }
                   | DELETE { getI :: i            }
                   | MATCH  { getI :: i, getJ :: j }
  deriving (Show, Eq, Ord, Generic, NFData)

isInsert, isDelete, isMatch :: Operation i j -> Bool

isInsert INSERT{} = True
isInsert _        = False

isDelete DELETE{} = True
isDelete _        = False

isMatch MATCH{} = True
isMatch _       = False

-- | Alignment matrix type
--
type Matrix m m' = UArray (Index m, Index m', EditOp) Int

-- | Traceback stop condition type
--
type Stop m m' = Matrix m m' -> m -> m' -> Index m -> Index m' -> Bool

-- | Traceback next move generator type
--
type Move m m'
  =  Matrix m m'
  -> m -> m'
  -> Index m -> Index m'
  -> EditOp -- ^ Current matrix, used in affine alignment
  -> (EditOp, Index m, Index m', Operation (Index m) (Index m'))
     -- ^ Next matrix, next indices, new operation

-- | A set of traceback conditions
--
data Conditions m m' = Conditions { isStop :: Stop m m' -- ^ Should we stop?
                                  , doMove :: Move m m' -- ^ Where to go next?
                                  }

-- | Sequence Alignment result
--
data AlignmentResult m m' = AlignmentResult { score     :: Int                              -- ^ Resulting score of alignment
                                            , alignment :: [Operation (Index m) (Index m')] -- ^ Alignment structure
                                            , sequence1 :: m                                -- ^ First chain
                                            , sequence2 :: m'                               -- ^ Second chain
                                            }
  deriving (Generic)

instance (NFData a, NFData b) => NFData (UArray (a, b, EditOp) Int) where
  rnf a = seq a ()

deriving instance (NFData a, NFData b, NFData (Index a), NFData (Index b))
         => NFData (AlignmentResult a b)

-- | Chain, that can be used for alignment
--
type Alignable m = (ChainLike m, Ix (Index m))

-- |Method of sequence alignment
--
class SequenceAlignment (a :: Type -> Type -> Type) where
    -- | Defines wheater the alignment is semiglobal or not
    --
    semi :: a e1 e2 -> Bool
    {-# INLINABLE semi #-}
    semi = const False

    -- | Traceback conditions of alignment
    --
    cond :: (Alignable m, Alignable m') => a (IxValue m) (IxValue m') -> Conditions m m'

    -- | Starting position in matrix for traceback procedure
    --
    traceStart :: (Alignable m, Alignable m') => a (IxValue m) (IxValue m') -> Matrix m m' -> m -> m' -> (Index m, Index m')

    -- | Distance matrix element
    --
    scoreMatrix :: (Alignable m, Alignable m') => a (IxValue m) (IxValue m') -> m -> m' -> Matrix m m'
