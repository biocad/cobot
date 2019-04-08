{-# LANGUAGE StandaloneDeriving #-}

module Bio.Chain.Alignment.Type where

import           Data.Array.ST                  (STUArray, Ix)
import           Control.Lens                   (Index, IxValue)
import           Bio.Chain                      (ChainLike (..))
import           Control.Monad.ST               (ST)

-- | Scoring function, returns substitution score for a couple of elements
--
type Scoring a b = a -> b -> Int

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
data Operation i j = Insert {            getJ :: j }
                   | Delete { getI :: i            }
                   | Match  { getI :: i, getJ :: j }
  deriving (Show, Eq, Ord)

isInsert, isDelete, isMatch :: Operation i j -> Bool

isInsert Insert{} = True
isInsert _        = False

isDelete Delete{} = True
isDelete _        = False

isMatch Match{} = True
isMatch _       = False

-- | Alignment matrix type
--
type Matrix s m m' e = STUArray s (Index m, Index m') e

-- | Sequence Alignment result
--
data AlignmentResult m m'
    = AlignmentResult { arScore         :: Int      -- ^ Resulting score of alignment
                      , arOperations    :: [Operation (Index m) (Index m')] -- ^ Alignment structure
                      , arFirstChain    :: m        -- ^ First chain
                      , arSecondChain   :: m'       -- ^ Second chain
                      , arMatchRange    :: ((Index m, Index m'), (Index m, Index m')) -- ^ Indices which edit operations affect
                      }

deriving instance (Show (Index m), Show (Index m'), Show m, Show m') => Show (AlignmentResult m m')

-- | Chain, that can be used for alignment
--
type Alignable m = (ChainLike m, Ix (Index m))

class SequenceAlignment (a :: * -> * -> *) where
    isAffine
        :: a (IxValue m) (IxValue m')
        -> m  -- fictional parameter to make type inference work
        -> m' -- fictional parameter to make type inference work
        -> Bool
    calculateMatrix
        :: (Alignable m, Alignable m')
        => a (IxValue m) (IxValue m')
        -> Matrix s m m' Int -- scores
        -> Matrix s m m' Int -- insertionCosts
        -> Matrix s m m' Int -- deletionCosts
        -> Matrix s m m' Int -- paths
        -> m  -- s
        -> m' -- t
        -> Index m  -- i
        -> Index m' -- j
        -> ST s (Int, Int, Int, Int)
    calculateStartPoint
        :: (Alignable m, Alignable m')
        => a (IxValue m) (IxValue m')
        -> Matrix s m m' Int -- scores
        -> m   -- fictional parameter to make type inference work
        -> m'  -- fictional parameter to make type inference work
        -> ST s (Index m, Index m')
    postProcessOperations
        :: (Alignable m, Alignable m')
        => a (IxValue m) (IxValue m')
        -> [Operation (Index m) (Index m')]
        -> (Index m, Index m')
        -> (Index m, Index m')
        -> m
        -> m'
        -> ST s ([Operation (Index m) (Index m')], (Index m, Index m'), (Index m, Index m'))
