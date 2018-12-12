{-# LANGUAGE RankNTypes       #-}
{-# LANGUAGE TypeFamilies     #-}
{-# LANGUAGE FlexibleContexts #-}
module Bio.Protein.Metric
    ( distance
    , angle
    , dihedral
    ) where

import           Control.Lens
import           Bio.Utils.Geometry             ( V3R
                                                , R
                                                )
import qualified Bio.Utils.Geometry            as G

distance :: Getting V3R a V3R -> Getting V3R a V3R -> Getting R a R
distance x y _ aa = Const $ G.distance (aa ^. x) (aa ^. y)

angle :: Getting V3R a V3R -> Getting V3R a V3R -> Getting V3R a V3R -> Getting R a R
angle x y z _ aa = Const $ G.angle (aa ^. x - aa ^. y) (aa ^. z - aa ^. y)

dihedral :: Getting V3R a V3R -> Getting V3R a V3R -> Getting V3R a V3R -> Getting V3R a V3R -> Getting R a R
dihedral x y z w _ aa = Const $ G.dihedral (aa ^. x) (aa ^. y) (aa ^. z) (aa ^. w)