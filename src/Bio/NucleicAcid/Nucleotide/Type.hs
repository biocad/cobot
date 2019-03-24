module Bio.NucleicAcid.Nucleotide.Type where

import Control.Lens (Iso', iso)

data DNA = DA | DC | DG | DT
  deriving (Eq, Ord, Bounded, Enum)

instance Show DNA where
    show DA = "Adenine"
    show DC = "Cytosine"
    show DG = "Guanine"
    show DT = "Thymine"

data RNA = RA | RC | RG | RU
  deriving (Eq, Ord, Bounded, Enum)

instance Show RNA where
    show RA = "Adenine"
    show RC = "Cytosine"
    show RG = "Guanine"
    show RU = "Uracil"

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
