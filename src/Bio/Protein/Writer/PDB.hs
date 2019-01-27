{-# OPTIONS_GHC -fno-warn-orphans #-}

module Bio.Protein.Writer.PDB
  (
    ToText (..)
  , ToFile (..)
  ) where

import           Bio.Protein.AminoAcid (AA (..), BBT, atom, c, ca, n,
                                        radicalType)
import           Bio.Utils.Geometry    (V3R)
import           Bio.Utils.IO.Class    (ToFile (..), ToText (..))
import           Bio.Utils.IO.PDB      (PDBLine (..), ToPDBLines (..))
import           Bio.Utils.Monomer     (ThreeSymbols (..))
import           Control.Lens          ((^.))
import           Data.Text             as T (Text, head, pack)
import           Linear.V3             (_x, _y, _z)

-- | TODO: change @N@, @CA@, @C@ to IUPAC names from module 'Bio.Utils.IUPAC'.
--
instance ToPDBLines (BBT V3R) where
  toPDBLines aa = [lineN, lineCA, lineC]
    where
      aaType = aa ^. radicalType
      lineN  = atomToPDBLine "N"  (aa ^. n  . atom) aaType
      lineCA = atomToPDBLine "CA" (aa ^. ca . atom) aaType
      lineC  = atomToPDBLine "C"  (aa ^. c  . atom) aaType

--------------------------------------------------------------------------------
-- INTERNAL
--------------------------------------------------------------------------------

-- | How to convert one atom into 'PDBLine'.
-- As you can see, we do not think about fields @atomSerial@, @atomChainID@, @atomResSeq@.
-- This fields will be updated in function from module 'Bio.Utils.IO.PDB'.
--
-- Moreover, as you can see, we use 'AA' to write @atomResName@.
-- This restriction is made consciously. We do not want to work with ACE and NMA and protonated forms of aminoacids.
--
atomToPDBLine :: Text -> V3R -> AA -> PDBLine
atomToPDBLine atomType coords aa =
    AtomLine { atomSerial     = 0
             , atomName       = atomType
             , atomAltLoc     = ' '
             , atomResName    = threeSymbols aa
             , atomChainID    = 'X'
             , atomResSeq     = 0
             , atomICode      = ' '
             , atomX          = coords ^. _x
             , atomY          = coords ^. _y
             , atomZ          = coords ^. _z
             , atomOccupancy  = 1.0
             , atomTempFactor = 0.0
             , atomElement    = pack [T.head atomType]
             , atomCharge     = "  "}
