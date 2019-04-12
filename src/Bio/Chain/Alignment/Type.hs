{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE DeriveAnyClass     #-}

module Bio.Chain.Alignment.Type where

import           Bio.Chain                      (ChainLike (..))
import           Control.DeepSeq                (NFData (..))
import           Control.Lens                   (Index, IxValue)
import           Control.Monad.ST               (ST)
import           Data.Array.ST                  (STUArray, Ix)
import           GHC.Generics                   (Generic)

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
  deriving (Show, Eq, Ord, Generic, NFData)

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
    = AlignmentResult { score      :: Int      -- ^ Resulting score of alignment
                      , alignment  :: [Operation (Index m) (Index m')] -- ^ Alignment structure
                      , sequence1  :: m        -- ^ First chain
                      , sequence2  :: m'       -- ^ Second chain
                      , matchRange :: ((Index m, Index m'), (Index m, Index m')) -- ^ Range of indices which edit operations affect
                      }

deriving instance (Show (Index m), Show (Index m'), Show m, Show m') => Show (AlignmentResult m m')

instance (NFData a, NFData b, NFData (Index a), NFData (Index b))
         => NFData (AlignmentResult a b) where
    rnf (AlignmentResult score alignment s1 s2 range) = rnf (score, alignment, s1, s2, range)

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
