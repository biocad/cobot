module Bio.NucleicAcid.Nucleotide.Type where

import Control.Lens (Iso', iso)

data DNA = DA | DC | DG | DT
  deriving (Eq, Ord, Bounded, Enum)

data RNA = RA | RC | RG | RU
  deriving (Eq, Ord, Bounded, Enum)

nucleoIso :: Iso' DNA RNA
nucleoIso = iso toRNA toDNA

toRNA :: DNA -> RNA
toRNA DA = RA
toRNA DC = RC
toRNA DG = RG
toRNA DT = RU

toDNA :: RNA -> DNA
toDNA RA = DA
toDNA RC = DC
toDNA RG = DG
toDNA RU = DT
