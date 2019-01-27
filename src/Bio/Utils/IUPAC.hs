module Bio.Utils.IUPAC
  (
    AtomType (..)
  ) where

import           Control.DeepSeq (NFData (..))
import           GHC.Generics    (Generic)

-- | Atom types in IUPAC nomenclature.
--
data AtomType = N   | CA  | C    | O  | OXT
              | CB
              | CG  | CG1 | CG2
              | CD  | CD1 | CD2
              | CE  | CE1 | CE2  | CE3
              | CH3
              | CZ  | CZ2 | CZ3
              | CH2
              | SG
              | SD
              | OG  | OG1
              | OD1 | OD2
              | OE1 | OE2
              | OH
              | ND1 | ND2
              | NE  | NE1 | NE2
              | NZ
              | NH1 | NH2
              | H
              | HA  | HA2 | HA3
              | HB  | HB1 | HB2  | HB3
              | HG  | HG1 | HG2  | HG3  | HG11 | HG12 | HG13 | HG21 | HG22 | HG23
              | HD  | HD1 | HD2  | HD3  | HD11 | HD12 | HD13 | HD21 | HD22 | HD23
              | HE  | HE1 | HE2  | HE3  | HE21 | HE22
              | HH  | HH2 | HH11 | HH12 | HH21 | HH22
              | HZ  | HZ1 | HZ2  | HZ3
  deriving (Show, Read, Eq, Ord, Generic)

instance NFData AtomType
