module Bio.NucleicAcid.Chain where

import           Bio.Chain                  as C
import           Control.Lens
import           Data.Array                 (Ix (..))

newtype NucleicAcidChain i a = NucleicAcidChain { getChain :: Chain i a }
  deriving (Show, Eq)

type instance Index (NucleicAcidChain i a) = i
type instance IxValue (NucleicAcidChain i a) = a

instance Ix i => Ixed (NucleicAcidChain i a) where
    ix i = lens (\ar -> getChain ar ^? ix i) (\ar (Just x) -> NucleicAcidChain (getChain ar & ix i .~ x)) . traverse

instance (Enum i, Ix i) => ChainLike (NucleicAcidChain i a) where
    bounds           (NucleicAcidChain ar) = bounds ar
    assocs           (NucleicAcidChain ar) = assocs ar
    modify       i f (NucleicAcidChain ar) = NucleicAcidChain $ modify i f ar
    modifyBefore i f (NucleicAcidChain ar) = NucleicAcidChain $ modifyBefore i f ar
    modifyAfter  i f (NucleicAcidChain ar) = NucleicAcidChain $ modifyAfter i f ar
