module Bio.Utils.Matrix
  ( eig3
  ) where

import Linear.Matrix
import Linear.V3

import Bio.Utils.Geometry ( R )

eig3 :: M33 R -> V3 R
eig3 m | isSym m = ei3sym m
       | otherwise = undefined
  where
    ei3sym :: M33 R -> V3 R
    ei3sym = undefined

isSym :: Eq a => M33 a -> Bool
isSym (V3 (V3 _   a21 a31)
          (V3 a12 _   a32)
          (V3 a13 a23 _  )) = a21 == a12 && a31 == a13 && a32 == a23