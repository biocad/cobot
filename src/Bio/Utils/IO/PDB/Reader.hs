{-# OPTIONS_GHC -fno-warn-orphans #-}

module Bio.Utils.IO.PDB.Reader
  ( -- * Reader
    -- $about
    FromPDBLines (..)
  ) where

import           Bio.Utils.IO.Class    (FromText (..))
import           Bio.Utils.IO.PDB.Type (PDBChainID, PDBLine (..))
import           Control.Arrow         (second, (&&&))
import           Control.Monad         (void)
import           Data.List             (groupBy, sortOn)
import           Data.Map              as M (Map, fromList)
import           Data.Text             as T (Text, pack, strip)
import           Text.Megaparsec       (ErrorFancy, Parsec, count, eof, many,
                                        parse, try, (<|>))
import           Text.Megaparsec.Char  (anyChar, char, noneOf, string)

{- $about

This module contains parser for 'PDBLine's.

Also, this module defines class 'FromPDBLines'.
If you have your own types, everything that your need is to define instance for 'FromPDBLines' for the type that corresponds residue.
Then you will have ability to read from any valid IO automatically.

TODO: more general instances for 'FromText'.

-}

-- | How to convert 'PDBLine's into residue in your own types.
--
class FromPDBLines a where
  fromPDBLines :: [PDBLine] -> a

instance FromText [PDBLine] where
  fromText = either (error . show) id . parse pdbLinesParser "[PDBLine] parser"

instance FromPDBLines a => FromText (M.Map PDBChainID [a]) where
  fromText = M.fromList . fmap splitToItems . splitToChains . fromText

splitToChains :: [PDBLine] -> [(PDBChainID, [PDBLine])]
splitToChains = fmap    extractChainID
              . groupBy sameChain
              . sortOn  atomChainID
              . filter  onlyAtomLine
  where
    onlyAtomLine :: PDBLine -> Bool
    onlyAtomLine AtomLine{} = True
    onlyAtomLine _          = False

    sameChain :: PDBLine -> PDBLine -> Bool
    sameChain a b = atomChainID a == atomChainID b

    extractChainID :: [PDBLine] -> (PDBChainID, [PDBLine])
    extractChainID = (atomChainID . head) &&& id

splitToItems :: FromPDBLines a => (PDBChainID, [PDBLine]) -> (PDBChainID, [a])
splitToItems = second (fmap fromPDBLines . groupBy sameResidue)
  where
    sameResidue :: PDBLine -> PDBLine -> Bool
    sameResidue a b = atomResSeq a == atomResSeq b &&
                      atomICode  a == atomICode  b


--------------------------------------------------------------------------------
-- PARSER
--------------------------------------------------------------------------------

type Parser = Parsec (ErrorFancy ()) Text

pdbLinesParser :: Parser [PDBLine]
pdbLinesParser = many (atomLineP <|> otherLineP <|> endFileP)

otherLineP :: Parser PDBLine
otherLineP = OtherLine . T.strip . T.pack <$> try (count 6 anyChar) <*> (pack <$> otherLineSymbols)

endFileP :: Parser PDBLine
endFileP = OtherLine <$> (string "END" <* (pack <$> otherLineSymbols)) <*> return ""

atomLineP :: Parser PDBLine
atomLineP = AtomLine <$> (try (string "ATOM ")                                    -- (1 -  6) # we increased atomSerial field for one symbol
                           *> (read <$> count 6 anyChar))                         -- (7 - 11)  atomSerial
                     <*> (T.strip . T.pack <$> (anyChar *> count 4 anyChar))      -- (13 - 16) atomName
                     <*> anyChar                                                  -- (17)      atomAltLoc
                     <*> (T.pack <$> count 3 anyChar)                             -- (18 - 20) atomResName
                     <*> (anyChar *> anyChar)                                     -- (22)      atomChainID
                     <*> (read <$> count 4 anyChar)                               -- (23 - 26) atomResSeq
                     <*> anyChar                                                  -- (27)      atomICode
                     <*> (read <$> (count 3 anyChar *> count 8 anyChar))          -- (31 - 38) atomX
                     <*> (read <$> count 8 anyChar)                               -- (39 - 46) atomY
                     <*> (read <$> count 8 anyChar)                               -- (47 - 54) atomZ
                     <*> (read <$> count 6 anyChar)                               -- (55 - 60) atomOccupancy
                     <*> (read <$> count 6 anyChar)                               -- (61 - 66) atomTempFactor
                     <*> (strip . pack <$> (count 10 anyChar *> count 2 anyChar)) -- (77 - 78) atomElement
                     <*> (strip . pack <$> count 2 anyChar                        -- (79 - 80) atomCharge
                            <* otherLineSymbols)

eol :: Parser ()
eol = void (char '\n') <|> eof

otherLineSymbols :: Parser String
otherLineSymbols = many (noneOf ['\n']) <* eol

-- atomLines = sortOn atomChainID $ filter onlyAtoms pdbLines

{-
moleculeToPDB :: forall g a . (Functor g, Foldable g, ToPDBLines a) => [g a] -> [PDBLine]
moleculeToPDB chains = renumerate . concat $ zipWith updateChainID chainCodes (convert <$> chains)
  where
    chainCodes :: [PDBChainID]
    chainCodes | length chains < 26 = take (length chains) ['A'..]
               | otherwise          = error "Bio.Protein.Writer.PDB.Molecule: impossible to save more than 26 chains"

    convert :: g a -> [[PDBLine]]
    convert = fmap toPDBLines . toList

renumerate :: [[PDBLine]] -> [PDBLine]
renumerate blocks = zipWith3 updateLine (concat blocks) [1..] resIndexes
  where
    lengths :: [Int]
    lengths = length <$> blocks

    resIndexes :: [Int]
    resIndexes = concat $ zipWith replicate lengths [1..]

    updateLine :: PDBLine -> Int -> Int -> PDBLine
    updateLine pdbLine atomIdx resIdx =
        pdbLine { atomSerial = atomIdx, atomResSeq = resIdx}

updateChainID :: PDBChainID -> [[PDBLine]] -> [[PDBLine]]
updateChainID chainID = fmap (fmap updateChainID')
  where
    updateChainID' :: PDBLine -> PDBLine
    updateChainID' pl = pl {atomChainID = chainID}
-}
