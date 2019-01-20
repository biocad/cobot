module Bio.Utils.Monomer
  ( Symbol(..)
  , ThreeSymbols(..)
  , FromThreeSymbols (..)
  ) where

import           Data.Text                      ( Text )

class Symbol a where
    symbol :: a -> Char

class ThreeSymbols a where
    threeSymbols :: a -> Text

instance Symbol Char where
    symbol = id

class FromThreeSymbols a where
    fromThreeSymbols :: Text -> a
