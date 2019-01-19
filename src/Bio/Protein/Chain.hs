module Bio.Protein.Chain
  ( module C
  , ProteinChain(..)
  ) where

import           Control.Lens
import           Bio.Chain                     as C

newtype ProteinChain a = ProteinChain { getChain :: Chain a }
  deriving (Show, Eq)

type instance Index (ProteinChain a) = Int
type instance IxValue (ProteinChain a) = a

instance Ixed (ProteinChain a) where
    ix i = lens (\ar -> getChain ar ^? ix i) (\ar (Just x) -> ProteinChain (getChain ar & ix i .~ x)) . traverse

instance ChainLike (ProteinChain a) where
    modify       i f (ProteinChain ar) = ProteinChain $ modify i f ar
    modifyBefore i f (ProteinChain ar) = ProteinChain $ modifyBefore i f ar
    modifyAfter  i f (ProteinChain ar) = ProteinChain $ modifyAfter i f ar
