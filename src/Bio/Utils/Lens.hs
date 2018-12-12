{-# LANGUAGE RankNTypes #-}
module Bio.Utils.Lens
    ( idx
    ) where

import           Control.Lens

idx :: Ixed m => Index m -> Lens' m (IxValue m)
idx i = lens (^?! ix i) (flip $ set (ix i))
