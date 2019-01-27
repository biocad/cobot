{-# OPTIONS_GHC -fno-warn-orphans #-}
{-# LANGUAGE AllowAmbiguousTypes   #-}
{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE ScopedTypeVariables   #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE TypeApplications #-}


module Bio.Utils.IO.PDB.Writer
  ( -- * Writer
    -- $about
    ToPDBLines (..)
  , chainsToPDB
  ) where

import           Bio.Utils.IO.Class       (ToText (..))
import           Bio.Utils.IO.PDB.Convert (pdbLineToText)
import           Bio.Utils.IO.PDB.Type    (PDBChainID,
                                           PDBLine (..))
import           Data.Foldable            (Foldable (..))
import           Data.List                (replicate)
import           Data.Text                as T (unlines)

{- $about 

This module contains convertion from 'PDBLine' to 'Text'.

Also, this module defines class 'ToPDBLines'.
If you have your own types, everything that your need is to define instance for 'ToPDBLines' for the type that corresponds residue.
Then you will have ability to write into any valid IO automatically.

-}

-- | How to convert residue in your own types into 'PDBLine's.
--
class ToPDBLines a where
  toPDBLines :: a -> [PDBLine]

instance (Foldable f, Functor g, Foldable g, ToPDBLines a) => ToText (f (g a)) where
  toText = T.unlines . fmap pdbLineToText . chainsToPDB . toList

-- | Converts list of chains @g@ with list of residues @a@ to list of 'PDBLine's.
-- TODO: write more about rules of the conversion.  
--
chainsToPDB :: forall g a . (Functor g, Foldable g, ToPDBLines a) => [g a] -> [PDBLine]
chainsToPDB chains = renumerate . concat $ zipWith updateChainID chainCodes (convert <$> chains)
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

