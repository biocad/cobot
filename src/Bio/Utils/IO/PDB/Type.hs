module Bio.Utils.IO.PDB.Type
  ( -- * Types
    PDBChainID
  , PDBLine (..)
    -- * Exceptions
  , PDBParseException (..)
  ) where

import           Control.Exception (Exception)
import           Data.Text         (Text)
import           Data.Typeable     (Typeable)

-- | Alias for chain identifier.
--
type PDBChainID = Char

-- | Type for one line in PDB format.
--
data PDBLine
    -- | [ATOM](http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM)
    = AtomLine   { atomSerial     :: Int        -- ^ Atom  serial number.
                 , atomName       :: Text       -- ^ Atom name.
                 , atomAltLoc     :: Char       -- ^ Alternate location indicator.
                 , atomResName    :: Text       -- ^ Residue name.
                 , atomChainID    :: PDBChainID -- ^ Chain identifier.
                 , atomResSeq     :: Int        -- ^ Residue sequence number.
                 , atomICode      :: Char       -- ^ Code for insertion of residues.
                 , atomX          :: Float      -- ^ Orthogonal coordinates for X in Angstroms.
                 , atomY          :: Float      -- ^ Orthogonal coordinates for Y in Angstroms.
                 , atomZ          :: Float      -- ^ Orthogonal coordinates for Z in Angstroms.
                 , atomOccupancy  :: Float      -- ^ Occupancy.
                 , atomTempFactor :: Float      -- ^ Temperature  factor.
                 , atomElement    :: Text       -- ^ Element symbol, right-justified.
                 , atomCharge     :: Text       -- ^ Charge  on the atom.
                 }
    -- | [SEQRES](http://www.wwpdb.org/documentation/file-format-content/format33/sect3.html#SEQRES)
    | SeqResLine { atomSerNum   :: Int     -- ^ Serial number of the SEQRES record for  the current  chain. Starts at 1 and increments by one each line. Reset to 1 for each chain.
                 , atomChainID  :: Char    -- ^ Chain identifier. This may be any single legal  character, including a blank which is used if there is only one chain.
                 , atomNumRes   :: Int     -- ^ Number of residues in the chain. This  value is repeated on every record.
                 , atomResidues :: [Text]  -- ^ Residue names.
                 }
    -- | [SHEET](http://www.wwpdb.org/documentation/file-format-content/format33/sect5.html#SHEET)
    | SheetShortLine  { sheetStrand      :: Int  -- ^ Strand  number which starts at 1 for each strand within a sheet and increases by one.
                      , sheetID          :: Text -- ^ Sheet  identifier.
                      , sheetNumStrands  :: Int  -- ^ Number  of strands in sheet.
                      , sheetInitResName :: Text -- ^ Residue  name of initial residue.
                      , sheetInitChainID :: Char -- ^ Chain identifier of initial residue in strand.
                      , sheetInitSeqNum  :: Int  -- ^ Sequence number of initial residue in strand.
                      , sheetInitICode   :: Char -- ^ Insertion code of initial residue in strand.
                      , sheetEndResName  :: Text -- ^ Residue  name of terminal residue.
                      , sheetEndChainID  :: Char -- ^ Chain identifier of terminal residue in strand.
                      , sheetEndSeqNum   :: Int  -- ^ Sequence number of terminal residue in strand.
                      , sheetEndICode    :: Char -- ^ Insertion code of terminal residue in strand.
                      , sheetSense       :: Int  -- ^ Sense of strand with respect to previous strand in the sheet. 0 if first strand, 1 if  parallel,and -1 if anti-parallel.
                      }
    -- | [SHEET](http://www.wwpdb.org/documentation/file-format-content/format33/sect5.html#SHEET)
    | SheetLine  { sheetStrand      :: Int  -- ^ Strand  number which starts at 1 for each strand within a sheet and increases by one.
                 , sheetID          :: Text -- ^ Sheet  identifier.
                 , sheetNumStrands  :: Int  -- ^ Number  of strands in sheet.
                 , sheetInitResName :: Text -- ^ Residue  name of initial residue.
                 , sheetInitChainID :: Char -- ^ Chain identifier of initial residue in strand.
                 , sheetInitSeqNum  :: Int  -- ^ Sequence number of initial residue in strand.
                 , sheetInitICode   :: Char -- ^ Insertion code of initial residue in strand.
                 , sheetEndResName  :: Text -- ^ Residue  name of terminal residue.
                 , sheetEndChainID  :: Char -- ^ Chain identifier of terminal residue in strand.
                 , sheetEndSeqNum   :: Int  -- ^ Sequence number of terminal residue in strand.
                 , sheetEndICode    :: Char -- ^ Insertion code of terminal residue in strand.
                 , sheetSense       :: Int  -- ^ Sense of strand with respect to previous strand in the sheet. 0 if first strand, 1 if  parallel,and -1 if anti-parallel.
                 , sheetCurAtom     :: Text -- ^ Registration.  Atom name in current strand.
                 , sheetCurResName  :: Text -- ^ Registration.  Residue name in current strand
                 , sheetCurChainId  :: Char -- ^ Registration. Chain identifier in current strand.
                 , sheetCurResSeq   :: Int  -- ^ Registration.  Residue sequence number in current strand.
                 , sheetCurICode    :: Char -- ^ Registration. Insertion code in current strand.
                 , sheetPrevAtom    :: Text -- ^ Registration.  Atom name in previous strand.
                 , sheetPrevResName :: Text -- ^ Registration.  Residue name in previous strand
                 , sheetPrevChainId :: Char -- ^ Registration. Chain identifier in previous strand.
                 , sheetPrevResSeq  :: Int  -- ^ Registration.  Residue sequence number in previous strand.
                 , sheetPrevICode   :: Char -- ^ Registration. Insertion code in previous strand.
                 }
    -- | [HELIX](http://www.wwpdb.org/documentation/file-format-content/format33/sect5.html#HELIX)
    | HelixLine  { helixSerNum      :: Int  -- ^ Serial number of the helix. This starts at 1  and increases incrementally.
                 , helixID          :: Text -- ^ Helix  identifier. In addition to a serial number, each helix is given an alphanumeric character helix identifier.
                 , helixInitResName :: Text -- ^ Name of the initial residue.
                 , helixInitChainID :: Char -- ^ Chain identifier for the chain containing this  helix.
                 , helixInitSeqNum  :: Int  -- ^ Sequence number of the initial residue.
                 , helixInitICode   :: Char -- ^ Insertion code of the initial residue.
                 , helixEndResName  :: Text -- ^ Name of the terminal residue of the helix.
                 , helixEndChainID  :: Char -- ^ Chain identifier for the chain containing this  helix.
                 , helixEndSeqNum   :: Int  -- ^ Sequence number of the terminal residue.
                 , helixEndICode    :: Char -- ^ Insertion code of the terminal residue.
                 , helixClass       :: Int  -- ^ Helix class (see below).
                 , helixComment     :: Text -- ^ Comment about this helix.
                 , helixLength      :: Int  -- ^ Length of this helix.
                 }
    -- | Any other line
    | OtherLine  { otherHeader :: Text
                 , otherText   :: Text
                 }
  deriving (Show, Eq)

instance Ord PDBLine where
  SeqResLine { atomSerNum = n1, atomChainID = c1 } <= SeqResLine { atomSerNum = n2, atomChainID = c2 }
    = c1 < c2 || (c1 == c2 && n1 <= n2)
  SheetShortLine { sheetStrand = s1, sheetInitChainID = c1 } <= SheetShortLine { sheetStrand = s2, sheetInitChainID = c2 }
    = c1 < c2 || (c1 == c2 && s1 <= s2)
  SheetLine { sheetStrand = s1, sheetInitChainID = c1 } <= SheetLine { sheetStrand = s2, sheetInitChainID = c2 }
    = c1 < c2 || (c1 == c2 && s1 <= s2)
  HelixLine { helixSerNum = n1, helixInitChainID = c1 } <= HelixLine { helixSerNum = n2, helixInitChainID = c2 }
    = c1 < c2 || (c1 == c2 && n1 <= n2)
  AtomLine { atomSerial = a1, atomChainID = c1 } <= AtomLine { atomSerial = a2, atomChainID = c2 }
    = c1 < c2 || (c1 == c2 && a1 < a2)

  SeqResLine     {} <= _ = True
  SheetShortLine {} <= SheetLine {} = True
  SheetShortLine {} <= HelixLine {} = True
  SheetShortLine {} <= OtherLine {} = True
  SheetLine      {} <= HelixLine {} = True
  SheetLine      {} <= OtherLine {} = True
  HelixLine      {} <= OtherLine {} = True
  OtherLine      {} <= AtomLine  {} = True
  OtherLine      {} <= OtherLine {} = True
  _                 <= AtomLine  {} = True

  _            <= SeqResLine     {} = False
  SheetLine {} <= SheetShortLine {} = False
  HelixLine {} <= SheetShortLine {} = False
  OtherLine {} <= SheetShortLine {} = False
  HelixLine {} <= SheetLine      {} = False
  OtherLine {} <= SheetLine      {} = False
  OtherLine {} <= HelixLine      {} = False
  AtomLine  {} <= OtherLine      {} = False
  AtomLine  {} <= _                 = False

-- | Parsing exception with some meta-info.
--
data PDBParseException
    = AtomIsAbsent      { atomsToSearchIn :: [String] -- ^ List of atoms to look for absent atom.
                        , absentAtom      :: String   -- ^ Absent atom.
                        }
    | AminoAcidIsAbsent { absentAminoAcid :: String
                        }
    deriving (Show, Typeable)

instance Exception PDBParseException
