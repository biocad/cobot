module Bio.Utils.Monomer
  ( Symbol(..)
  , FromSymbol(..)
  , ThreeSymbols(..)
  , FromThreeSymbols (..)
  ) where

import           Data.Text                      ( Text )

class Symbol a where
    symbol :: a -> Char

instance Symbol Char where
    symbol = id

class FromSymbol a where
    fromSymbol :: Char -> Maybe a
    fromSymbol = either (const Nothing) Just . fromSymbolE

    fromSymbolE :: Char -> Either Char a

instance FromSymbol Char where
    fromSymbolE = Right

class ThreeSymbols a where
    threeSymbols :: a -> Text

class FromThreeSymbols a where
    fromThreeSymbols :: Text -> Maybe a
