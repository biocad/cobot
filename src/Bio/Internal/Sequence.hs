module Bio.Internal.Sequence where

import           Data.Text                      ( Text )

class Symbol a where
    symbol :: a -> Char

class ThreeSymbols a where
    threeSymbols :: a -> Text