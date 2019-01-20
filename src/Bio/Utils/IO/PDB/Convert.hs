module Bio.Utils.IO.PDB.Convert
  ( pdbLineToText
  ) where

import           Bio.Utils.IO.PDB.Type (PDBLine (..))
import           Data.Char             (isDigit)
import           Data.List             (replicate)
import           Data.Monoid           ((<>))
import           Data.Text             (Text, append, concat, cons, head,
                                        length, pack, singleton, take)
import           Prelude               hiding (concat, head, length, take)
import           Text.Printf           (printf)

-- | Converts 'PDBLine' to 'Text'.
--
pdbLineToText :: PDBLine -> Text
pdbLineToText (AtomLine serial name altLoc resName chainID resSeq iCode x y z occupancy tempFactor element charge) =
  concat ["ATOM  ",
    fixedLine 5 (pack . show $ serial),
    " ",
    fixedAtomName name,
    singleton altLoc,
    fixedLine 3 resName,
    " ",
    singleton chainID,
    fixedLine 4 (pack . show $ resSeq),
    singleton iCode,
    "   ",
    fixedLine 8 (pack $ printf "%.3f" x),
    fixedLine 8 (pack $ printf "%.3f" y),
    fixedLine 8 (pack $ printf "%.3f" z),
    fixedLine 6 (pack $ printf "%f" occupancy),
    fixedLine 6 (pack $ printf "%f" tempFactor),
    "          ",
    fixedLine 2 element,
    fixedLine 2 charge
  ]
pdbLineToText (SeqResLine serNum chainID numRes residues) =
  fixedLineR 80 $ concat ["SEQRES ",
    fixedLine 3 (pack . show $ serNum),
    " ",
    singleton chainID,
    " ",
    fixedLine 4 (pack . show $ numRes),
    " ",
    concat $ fixedLine 4 <$> residues
  ]
pdbLineToText (SheetShortLine strand sheetId numStrands
                      initResName initChainID initSeqNum initICode
                      endResName endChainID endSeqNum endICode sense) =
  fixedLineR 80 $ concat ["SHEET  ",
    fixedLine 3 $ pack $ show strand,
    " ",
    fixedLine 3 sheetId,
    fixedLine 2 $ pack $ show numStrands,
    " ",
    fixedLine 3 initResName,
    " ",
    singleton initChainID,
    fixedLine 4 $ pack $ show initSeqNum,
    singleton initICode,
    " ",
    fixedLine 3 endResName,
    " ",
    singleton endChainID,
    fixedLine 4 $ pack $ show endSeqNum,
    singleton endICode,
    fixedLine 2 $ pack $ show sense
  ]
pdbLineToText (SheetLine strand sheetId numStrands
                      initResName initChainID initSeqNum initICode
                      endResName endChainID endSeqNum endICode
                      sense curAtom curResName curChainId curResSeq curICode
                      prevAtom prevResName prevChainId prevResSeq prevICode) =
  fixedLineR 80 $ concat ["SHEET  ",
    fixedLine 3 $ pack $ show strand,
    " ",
    fixedLine 3 sheetId,
    fixedLine 2 $ pack $ show numStrands,
    " ",
    fixedLine 3 initResName,
    " ",
    singleton initChainID,
    fixedLine 4 $ pack $ show initSeqNum,
    singleton initICode,
    " ",
    fixedLine 3 endResName,
    " ",
    singleton endChainID,
    fixedLine 4 $ pack $ show endSeqNum,
    singleton endICode,
    fixedLine 2 $ pack $ show sense,
    " ",
    fixedLine 4 curAtom,
    fixedLine 3 curResName,
    " ",
    singleton curChainId,
    fixedLine 4 $ pack $ show curResSeq,
    singleton curICode,
    " ",
    fixedLine 4 prevAtom,
    fixedLine 3 prevResName,
    " ",
    singleton prevChainId,
    fixedLine 4 $ pack $ show prevResSeq,
    singleton prevICode
  ]
pdbLineToText (HelixLine serNum hID
                      initResName initChainID initSeqNum initICode
                      endResName endChainID endSeqNum endICode
                      hClass comment len) =
  fixedLineR 80 $ concat ["HELIX  ",
    fixedLine 3 $ pack $ show serNum,
    " ",
    fixedLine 3 hID,
    " ",
    fixedLine 3 initResName,
    " ",
    singleton initChainID,
    " ",
    fixedLine 4 $ pack $ show initSeqNum,
    singleton initICode,
    " ",
    fixedLine 3 endResName,
    " ",
    singleton endChainID,
    " ",
    fixedLine 4 $ pack $ show endSeqNum,
    singleton endICode,
    fixedLine 2 $ pack $ show hClass,
    fixedLine 30 comment,
    " ",
    fixedLine 5 $ pack $ show len
  ]
pdbLineToText OtherLine {} = ""

fixedLine :: Int -> Text -> Text
fixedLine n s = if spacesCnt <= 0 then take n s
                else pack (replicate spacesCnt ' ') <> s
  where spacesCnt = n - length s

fixedLineR :: Int -> Text -> Text
fixedLineR n s = if spacesCnt <= 0 then s
                 else s <> pack (replicate spacesCnt ' ')
  where spacesCnt = n - length s

fixedAtomName :: Text -> Text
fixedAtomName s | lenS == 4 = s
                | isDigit f = s `append` pack (replicate (4 - lenS) ' ')
                | otherwise = ' ' `cons` s `append` pack (replicate (3 - lenS) ' ')
  where lenS = length s
        f = head s
