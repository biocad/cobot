{-# LANGUAGE DeriveAnyClass #-}

module Bio.NucleicAcid.Nucleotide.Type
  ( DNA (..)
  , RNA (..)
  , nucleoIso
  , toRNA
  , toDNA
  , Complementary (..)
  ) where

import Control.DeepSeq (NFData)
import Control.Lens    (Iso', iso)
import Data.Array      (Array, Ix, bounds, ixmap, listArray)
import Data.Foldable   (Foldable (..))
import GHC.Generics    (Generic)

data DNA
  = DA
  | DC
  | DG
  | DT
  deriving (Eq, Ord, Bounded, Enum, Generic, NFData)

instance Show DNA where
    show DA = "Adenine"
    show DC = "Cytosine"
    show DG = "Guanine"
    show DT = "Thymine"

data RNA
  = RA
  | RC
  | RG
  | RU
  deriving (Eq, Ord, Bounded, Enum, Generic, NFData)

instance Show RNA where
    show RA = "Adenine"
    show RC = "Cytosine"
    show RG = "Guanine"
    show RU = "Uracil"

-------------------------------------------------------------------------------
-- Transciption
-------------------------------------------------------------------------------

nucleoIso :: Iso' DNA RNA
nucleoIso = iso toRNA toDNA

{-# INLINE toRNA #-}
toRNA :: DNA -> RNA
toRNA DA = RA
toRNA DC = RC
toRNA DG = RG
toRNA DT = RU

{-# INLINE toDNA #-}
toDNA :: RNA -> DNA
toDNA RA = DA
toDNA RC = DC
toDNA RG = DG
toDNA RU = DT

------------------------------------------------------------------------------
-- Complementary and reverse complementary
-------------------------------------------------------------------------------

class Complementary a where
    -- | complement *NA (DNA or RNA)
    --
    cNA :: a -> a

    -- | reverce complement *NA (DNA or RNA)
    --
    rcNA :: a -> a

instance Complementary DNA where
    cNA DA = DT
    cNA DC = DG
    cNA DG = DC
    cNA DT = DA

    rcNA = cNA

instance Complementary RNA where
    cNA = toRNA . cNA . toDNA

    rcNA = cNA

instance Complementary Char where
    cNA 'A' = 'T'
    cNA 'C' = 'G'
    cNA 'G' = 'C'
    cNA 'T' = 'A'
    cNA 'a' = 't'
    cNA 'c' = 'g'
    cNA 'g' = 'c'
    cNA 't' = 'a'
    cNA c = c

    rcNA = cNA

instance Complementary a => Complementary [a] where
   cNA = fmap cNA

   rcNA = reverse . cNA

instance {-# OVERLAPPABLE #-} (Complementary a, Ix i) => Complementary (Array i a) where
   cNA = fmap cNA

   rcNA l = listArray (bounds l) rl
     where
       rl = rcNA . toList $ l

instance {-# OVERLAPPING  #-} (Complementary a) => Complementary (Array Int a) where
   cNA = fmap cNA

   rcNA l = cNA $ ixmap (lo, hi) (\i -> lo + (hi - i)) l
     where
       (lo, hi) = bounds l

