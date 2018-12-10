{-# LANGUAGE RankNTypes        #-}

module Bio.Internal.Sequence where

import           Control.Lens
import           Data.Text                      ( Text )

class Symbol a where
    symbol :: a -> Char

class ThreeSymbols a where
    threeSymbols :: a -> Text

idx :: Ixed m => Index m -> Lens' m (IxValue m)
idx i = lens (^?! ix i) (flip $ set (ix i))
