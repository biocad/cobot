{-# OPTIONS_GHC -fno-warn-orphans #-}

module Bio.Protein.Chain
  ( module C
  , ProteinChain(..)
  ) where

import           Bio.Chain             as C
import           Bio.Protein.AminoAcid
import           Control.Lens
import           Data.Array            (Ix (..))
import           Data.String           (IsString (..))

newtype ProteinChain i a = ProteinChain { getChain :: Chain i a }
  deriving (Show, Eq, Functor, Foldable, Traversable, ChainLike)

type instance Index (ProteinChain i a) = i
type instance IxValue (ProteinChain i a) = a

instance Ix i => Ixed (ProteinChain i a) where
    ix i' = coerced . ix @(Chain i a) i'

instance IsString (ProteinChain Int AA) where
  fromString = ProteinChain . fromString

