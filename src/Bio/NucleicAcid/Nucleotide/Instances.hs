{-# OPTIONS_GHC -fno-warn-orphans  #-}

module Bio.NucleicAcid.Nucleotide.Instances () where

import           Bio.NucleicAcid.Nucleotide.Type
import           Bio.Utils.Monomer               (FromSymbol (..), Symbol (..))
import           Data.Array                      (Array, listArray)
import           Data.String                     (IsString (..))

-------------------------------------------------------------------------------
-- Symbol and ThreeSymbols
-------------------------------------------------------------------------------

instance Symbol DNA where
    symbol DA = 'A'
    symbol DC = 'C'
    symbol DG = 'G'
    symbol DT = 'T'

instance FromSymbol DNA where
    fromSymbolE 'A' = Right DA
    fromSymbolE 'C' = Right DC
    fromSymbolE 'G' = Right DG
    fromSymbolE 'T' = Right DT
    fromSymbolE ch  = Left ch

instance Symbol RNA where
    symbol RA = 'A'
    symbol RC = 'C'
    symbol RG = 'G'
    symbol RU = 'U'

instance FromSymbol RNA where
    fromSymbolE 'A' = Right RA
    fromSymbolE 'C' = Right RC
    fromSymbolE 'G' = Right RG
    fromSymbolE 'U' = Right RU
    fromSymbolE ch  = Left ch

-------------------------------------------------------------------------------
-- IsString
-------------------------------------------------------------------------------

instance {-# OVERLAPPING #-} IsString [DNA] where
    fromString s =
        case traverse fromSymbolE s of
            Right l -> l
            Left e  -> error $ "Bio.NucleicAcid.Nucleotide.Instances: could not read nucleotide " <> [e]

instance IsString (Array Int DNA) where
    fromString s = listArray (0, length s - 1) $ fromString s

instance {-# OVERLAPPING #-} IsString [RNA] where
    fromString s =
        case traverse fromSymbolE s of
            Right l -> l
            Left e  -> error $ "Bio.NucleicAcid.Nucleotide.Instances: could not read nucleotide " <> [e]

instance IsString (Array Int RNA) where
    fromString s = listArray (0, length s - 1) $ fromString s
