{-# LANGUAGE RankNTypes   #-}
{-# LANGUAGE TypeFamilies #-}
module Bio.Protein.Metric where

import           Data.Coerce
import           Control.Lens
import           Linear.V3
import           Linear.Metric                 ( dot, normalize )
import qualified Linear.Metric                 as L
                                                ( distance )
import           Bio.Internal.Structure
import           Bio.Protein.AminoAcid

distance :: (v3r ~ V3R) => Getting v3r a v3r -> Getting v3r a v3r -> Getting R a R
distance x y _ aa = Const $ L.distance (coerce $ aa ^. x) (aa ^. y)

angle :: (v3r ~ V3R) => Getting v3r a v3r -> Getting v3r a v3r -> Getting v3r a v3r -> Getting R a R
angle x y z _ aa = Const $ angleBetween (aa ^. x - aa ^. y) (aa ^. z - aa ^. y)

dihedral :: (v3r ~ V3R) => Getting v3r a v3r -> Getting v3r a v3r -> Getting v3r a v3r -> Getting v3r a v3r -> Getting R a R
dihedral x y z w _ aa = let b1 = (aa ^. y - aa ^. x)
                            b2 = (aa ^. z - aa ^. y)
                            b3 = (aa ^. w - aa ^. z)
                            n1 = normalize $ b1 `cross` b2
                            n2 = normalize $ b2 `cross` b3
                            m1 = n1 `cross` normalize b2
                        in  Const $ atan2 (m1 `dot` n2) (n1 `dot` n2)

phi :: (Ixed m, Functor f, Functor r, Functor g, v3r ~ V3R, IxValue m ~ AminoAcid f r g v3r) => m -> Lens' m R
phi _ = undefined

psi :: (Ixed m, Functor f, Functor r, Functor g, v3r ~ V3R, IxValue m ~ AminoAcid f r g v3r) => m -> Lens' m R
psi _ = undefined

omega :: (Ixed m, Functor f, Functor r, Functor g, v3r ~ V3R, IxValue m ~ AminoAcid f r g v3r) => m -> Lens' m R
omega _ = undefined