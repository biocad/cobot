{-# OPTIONS_GHC -fno-warn-orphans #-}
module Bio.Protein.Chain
  ( module C
  , ProteinChain(..)
  ) where

import           Data.String                   ( IsString (..) )
import           Data.Array                    ( Ix (..), listArray )
import           Control.Lens
import           Bio.Chain                     as C
import           Bio.Protein.AminoAcid

newtype ProteinChain i a = ProteinChain { getChain :: Chain i a }
  deriving (Show, Eq)

instance IsString [AA] where
  fromString = fmap fromS

instance IsString (ProteinChain Int AA) where
  fromString s = ProteinChain $ listArray (0, length s - 1) (fromS <$> s)

type instance Index (ProteinChain i a) = i
type instance IxValue (ProteinChain i a) = a

instance Ix i => Ixed (ProteinChain i a) where
    ix i = lens (\ar -> getChain ar ^? ix i) (\ar (Just x) -> ProteinChain (getChain ar & ix i .~ x)) . traverse

instance (Enum i, Ix i) => ChainLike (ProteinChain i a) where
    bounds           (ProteinChain ar) = bounds ar
    assocs           (ProteinChain ar) = assocs ar
    modify       i f (ProteinChain ar) = ProteinChain $ modify i f ar
    modifyBefore i f (ProteinChain ar) = ProteinChain $ modifyBefore i f ar
    modifyAfter  i f (ProteinChain ar) = ProteinChain $ modifyAfter i f ar

-- Helper functions

fromS :: Char -> AA
fromS 'A' = ALA
fromS 'C' = CYS
fromS 'D' = ASP
fromS 'E' = GLU
fromS 'F' = PHE
fromS 'G' = GLY
fromS 'H' = HIS
fromS 'I' = ILE
fromS 'K' = LYS
fromS 'L' = LEU
fromS 'M' = MET
fromS 'N' = ASN
fromS 'P' = PRO
fromS 'Q' = GLN
fromS 'R' = ARG
fromS 'S' = SER
fromS 'T' = THR
fromS 'V' = VAL
fromS 'W' = TRP
fromS 'Y' = TYR
fromS x   = error $ "'" ++ show x ++ "' is not an amino acid"