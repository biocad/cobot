{-# OPTIONS_GHC -fno-warn-orphans #-}

module Bio.NucleicAcid.Chain where

import           Bio.Chain                  as C
import           Bio.NucleicAcid.Nucleotide
import           Control.Lens
import           Data.Array                 (Ix (..))
import           Data.String                (IsString (..))

newtype NucleicAcidChain i a = NucleicAcidChain { getChain :: Chain i a }
  deriving (Show, Eq, Functor, Foldable, Traversable)

instance IsString (NucleicAcidChain Int DNA) where
  fromString = NucleicAcidChain . fromString

instance IsString (NucleicAcidChain Int RNA) where
  fromString = NucleicAcidChain . fromString

type instance Index (NucleicAcidChain i a) = i
type instance IxValue (NucleicAcidChain i a) = a

instance Ix i => Ixed (NucleicAcidChain i a) where
    ix i' = coerced . ix @(Chain i a) i'

instance (Enum i, Ix i) => ChainLike (NucleicAcidChain i a) where
    bounds           (NucleicAcidChain ar) = bounds ar
    assocs           (NucleicAcidChain ar) = assocs ar
    modify       i f (NucleicAcidChain ar) = NucleicAcidChain $ modify i f ar
    modifyBefore i f (NucleicAcidChain ar) = NucleicAcidChain $ modifyBefore i f ar
    modifyAfter  i f (NucleicAcidChain ar) = NucleicAcidChain $ modifyAfter i f ar
