{-# OPTIONS_GHC -fno-warn-orphans  #-}

module Bio.NucleicAcid.Nucleotide.Instances () where

import           Bio.NucleicAcid.Nucleotide.Type
import           Bio.Utils.Monomer

instance Symbol DNA where
  symbol DA = 'A'
  symbol DC = 'C'
  symbol DG = 'G'
  symbol DT = 'T'

instance Symbol RNA where
  symbol RA = 'A'
  symbol RC = 'C'
  symbol RG = 'G'
  symbol RU = 'U'
