{- |

This description will be expanded later.

Minimum example (with writer example too):

> >>> import Data.Map as M
> >>> loadedChains <- fromFile "/tmp/1igt.pdb" :: IO (M.Map PDBChainID [BBT V3R])
> >>> keys loadedChains
> "ABCD"
> >>> loadedChains `toFile` "/tmp/myOwn.pdb"

-}

{-# OPTIONS_GHC -fno-warn-orphans #-}

module Bio.Protein.Reader.PDB
  ( FromText (..)
  , FromFile (..)
  ) where

import           Bio.Protein.AminoAcid (AA (..), BBT, Createable (..))
import           Bio.Utils.Geometry    (V3R)
import           Bio.Utils.IO.Class    (FromFile (..), FromText (..))
import           Bio.Utils.IO.PDB      (FromPDBLines (..), PDBLine (..))
import           Bio.Utils.IUPAC       (AtomType (..))
import           Bio.Utils.Monomer     (FromThreeSymbols (..))
import           Control.Arrow         ((&&&))
import           Data.Map.Strict       as M (Map, findWithDefault, fromList)
import           Data.Text             as T (Text, pack, unpack)
import           Linear.V3             (V3 (..))

instance FromPDBLines (BBT V3R) where
  fromPDBLines lines' = create @(BBT V3R) x y z aa'
    where
      atomTypes :: [AtomType]
      atomTypes = [N, CA, C]

      map' = extractCoords . linesToMap atomTypes $ lines'

      aa'  :: AA
      aa'  = fromThreeSymbols . atomResName . head $ lines'

      [x, y, z] = findAtoms map' atomTypes

linesToMap :: [AtomType] -> [PDBLine] -> Map AtomType PDBLine
linesToMap atomTypes = fromList
                     . fmap (read . unpack . atomName &&& id)
                     . filter byAtomType
  where
    atomTypesText :: [Text]
    atomTypesText = pack . show <$> atomTypes

    byAtomType :: PDBLine -> Bool
    byAtomType AtomLine{..} = atomName `elem` atomTypesText
    byAtomType _            = error "Bio.Protein.Reader.PDB.linesToMap: impossible"

extractCoords :: Map AtomType PDBLine -> Map AtomType V3R
extractCoords = fmap (\AtomLine{..} -> V3 atomX atomY atomZ)

findAtoms :: forall a f . Functor f => Map AtomType a -> f AtomType -> f a
findAtoms map' = fmap findOrCrash
  where
    findOrCrash :: AtomType -> a
    findOrCrash at = let err = error $ "could not find atom type " ++ show at ++ "in PDBLines"
                     in M.findWithDefault err at map'
