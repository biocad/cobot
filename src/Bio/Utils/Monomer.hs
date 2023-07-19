module Bio.Utils.Monomer
  ( Symbol(..)
  , FromSymbol(..)
  , ThreeSymbols(..)
  , FromThreeSymbols (..)
  , FullName(..)

  , ShowMonomer (..)
  ) where

import           Data.Text                      ( Text )
import           Data.Coerce                    ( coerce )

class Symbol a where
    symbol :: a -> Char

instance Symbol Char where
    symbol = id

class FromSymbol a where
    fromSymbol :: Char -> Maybe a
    fromSymbol = either (const Nothing) Just . fromSymbolE

    fromSymbolE :: Char -> Either Char a
    fromSymbolE c = maybe (Left c) Right $ fromSymbol c

instance FromSymbol Char where
    fromSymbolE = Right

class ThreeSymbols a where
    threeSymbols :: a -> Text

class FromThreeSymbols a where
    fromThreeSymbols :: Text -> Maybe a

-- | Full name of monomer, like "Adenosine" or "Lysine".
class FullName a where
    fullName :: a -> String

-- | DerivingVia wrapper to 'show' @[a]@ as a string of one-symbol encoding.
newtype ShowMonomer a = ShowMonomer a

instance Symbol a => Show (ShowMonomer a) where
    show (ShowMonomer nc) = [symbol nc]
    showList = foldMap (showChar . symbol . coerce @_ @a)
