{-# OPTIONS_GHC -fno-warn-orphans #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE FlexibleInstances    #-}
{-# LANGUAGE TypeFamilies         #-}

module Bio.Protein.Builder where

import           Data.Array
import           Control.Lens
import           Linear.V3
import           Linear.Vector
import           Linear.Metric

import           Bio.Internal.Structure
import           Bio.Protein.AminoAcid

instance Measureable BackboneAtom where
    dist N CA = 1.460
    dist CA C = 1.509
    dist C N  = 1.290
    dist x y  = dist y x

    angle N CA C = pi * 110.990 / 180.0
    angle CA C N = pi * 118.995 / 180.0
    angle C N CA = angle CA C N
    angle x y z  = angle z y x

class Buildable a where
    type Monomer a :: *

    initB :: Monomer a -> a
    nextB :: Monomer a -> a -> a

build :: Buildable a => [Monomer a] -> Array Int a
build chain = result
  where
    result = array (0, length chain - 1) [(i, next i a) | (a, i) <- chain `zip` [0..]]
    next 0 a = initB a
    next k a = nextB a (result ! (k - 1))
    
instance Buildable (BB V3R) where
    type Monomer (BB V3R) = AA

    -- | Place first amino acid backbone in some chain
    -- The placement will be like this:
    --        y /|\
    --           |
    --           |
    --      N    | Ca
    -- ----*-----*------------->
    --           |     C        x
    --           |    *
    --           |
    --
    initB _ = let n_ = V3 n_x 0.0 0.0
                  a_ = V3 0.0 0.0 0.0
                  c_ = V3 c_x c_y 0.0
                  --
                  n_x = - dist N CA
                  c_x = dist CA C * cos (pi + angle N CA C)
                  c_y = dist CA C * sin (pi + angle N CA C)
              in  makeBB n_ a_ c_

    -- | Place next amino acid backbone in some chain
    -- The placement can be done by two cases.
    -- First:
    --               Ca_i      N_i+1     C_i+1
    --              *         *         *
    --                
    --         *         *         *
    --          N_i       C_i       Ca_i+1
    -- Second:
    --          N_i       C_i       Ca_i+1
    --         *         *         *
    --
    --              *         *         *
    --               Ca_i      N_i+1     C_i+1
    --
    -- Let us enumerate atoms: 1 for N_i, 2 for Ca_i, 3 for C_i, 4 for N_i+1, 5 for Ca_i+1, 6 for C_i+1.
    -- We have to find points 4, 5, 6 using 1, 2, 3. To find this points let us introduce vectors named
    -- like 'vij' from i to j, e.g. v12 is a vector from N_i to Ca_i. Our main idea will be to get a 
    -- direction vector from i+1 to i, rotate it and then upscale by specified bond length. One thing to
    -- look at is the direction of rotations. If we have the first case, then the first rotation should be
    -- conterclock-wise, otherwise — clock-wise. To detect it we have to understand whether 3 is on the left
    -- of 12 vector (first case) or on the right. We can understand it using v21 and v23:
    -- if (v21 `cross` v23) ^. _z < 0 then First else Second. First means that every angle should be negated.
    -- So, we can determine coordinate of 4. First we get the v32 and normalize it, then we will rotate it to
    -- CA-C-N angle (multiplied by -1 or not), next multiply this direction vector by typical C-N bond length
    -- and at last add the obtained vector to 3. The same idea is used to find point 5, but now we should
    -- make out rotation in the opposite direction. At last we will do the same with point 6.
    --
    nextB _ aa = let -- we will always rotate around Z
                     rot = rotate (unit _z)
                     -- determine the direction
                     v21 = aa ^. n - aa ^. ca
                     v23 = aa ^. c - aa ^. ca
                     cw  = if (v21 `cross` v23) ^. _z < 0 then -1.0 else 1.0 :: R
                     -- determine the coordinate of n (point 4)
                     v32 = negated v23
                     v34 = dist C N *^ rot (cw * angle CA C N) (normalize v32)
                     n_  = aa ^. c + v34
                     -- determine the coordinate of ca (point 5)
                     v43 = negated v34
                     v45 = dist N CA *^ rot (-cw * angle C N CA) (normalize v43)
                     ca_ = n_ + v45
                     -- determine the coordinate of ca (point 6)
                     v54 = negated v45
                     v56 = dist CA C *^ rot (cw * angle N CA C) (normalize v54)
                     c_  = ca_ + v56
                 in  makeBB n_ ca_ c_

instance Buildable (BBT V3R) where
    type Monomer (BBT V3R) = AA

    initB t = let aa = initB t :: BB V3R
              in  makeBBT (aa ^. n) (aa ^. ca) (aa ^. c) t

    nextB t aaT = let aa = makeBB (aaT ^. n) (aaT ^. ca) (aaT ^. c)
                      ab = nextB t aa :: BB V3R
                  in  makeBBT (ab ^. n) (ab ^. ca) (ab ^. c) t
