{-# LANGUAGE TypeSynonymInstances #-}

module Bio.Protein.Chain.Builder
    ( Buildable (..)
    , build
    ) where

import           Data.Ix                        ( Ix )
import           Control.Lens
import           Linear.V3                      ( V3 (..)
                                                , cross
                                                , _z
                                                )
import           Linear.Vector                  ( negated
                                                , unit
                                                , (*^)
                                                )

import           Bio.Utils.Geometry      hiding ( angle )
import           Bio.Protein.AminoAcid
import           Bio.Protein.Chain

class Buildable a where
    type Monomer a :: *
    initB :: Monomer a -> a
    nextB :: Monomer a -> a -> a

build :: forall a m.(Buildable a, ChainLike m, Ix (Index m), IxValue m ~ Monomer a) => m -> ProteinChain (Index m) a
build ch = ProteinChain result
  where
    result :: Chain (Index m) a
    result = chain (bounds ch) [ (i, next i x) | (i, x) <- assocs ch ]
    next :: Index m -> Monomer a -> a
    next k x | k == fst (bounds ch) = initB x
             | otherwise            = nextB x (result ! pred k)

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
              in  create @(BB V3R) n_ a_ c_

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
                     v21 = aa ^. n . atom - aa ^. ca . atom
                     v23 = aa ^. c . atom - aa ^. ca . atom
                     cw  = if (v21 `cross` v23) ^. _z < 0 then 1.0 else -1.0 :: R
                     -- determine the coordinate of n (point 4)
                     v32 = negated v23
                     v34 = dist C N *^ rot (cw * angle CA C N) (normalize v32)
                     n_  = aa ^. c . atom + v34
                     -- determine the coordinate of ca (point 5)
                     v43 = negated v34
                     v45 = dist N CA *^ rot (-cw * angle C N CA) (normalize v43)
                     ca_ = n_ + v45
                     -- determine the coordinate of ca (point 6)
                     v54 = negated v45
                     v56 = dist CA C *^ rot (cw * angle N CA C) (normalize v54)
                     c_  = ca_ + v56
                 in  create @(BB V3R) n_ ca_ c_

instance Buildable (BBT V3R) where
    type Monomer (BBT V3R) = AA

    initB t = let aa = initB t :: BB V3R
              in  create @(BBT V3R) (aa ^. n . atom) (aa ^. ca . atom) (aa ^. c . atom) t

    nextB t aaT = let aa = create @(BB V3R) (aaT ^. n . atom) (aaT ^. ca . atom) (aaT ^. c . atom)
                      ab = nextB t aa :: BB V3R
                  in  create @(BBT V3R) (ab ^. n . atom) (ab ^. ca . atom) (ab ^. c . atom) t

-- Helper types and functions

-- | Atoms of amino acid backbone
--
data BackboneAtom = N | CA | C
  deriving (Show, Eq, Ord, Bounded, Enum)

-- | Atoms of amino acid radicals (TODO: fill this)
--
-- data RadicalAtom

-- | Distance between two basic backbone atom types
dist :: BackboneAtom -> BackboneAtom -> R
dist N  CA = 1.460
dist CA C  = 1.509
dist C  N  = 1.290
dist x  y  = dist y x

-- | Angles between every triple of succesive atoms
angle :: BackboneAtom -> BackboneAtom -> BackboneAtom -> R
angle N  CA C  = pi * 110.990 / 180.0
angle CA C  N  = pi * 118.995 / 180.0
angle C  N  CA = angle CA C N
angle x  y  z  = angle z y x
