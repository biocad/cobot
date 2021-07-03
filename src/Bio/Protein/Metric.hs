module Bio.Protein.Metric
    ( Metricable (..)
    ) where

import           Data.Kind                      ( Type )
import           Data.Monoid                    ( First (..) )
import           Control.Lens
import           Bio.Utils.Geometry             ( V3R
                                                , R
                                                )
import qualified Bio.Utils.Geometry            as G

class Metricable m where
    type ReturnMetric m :: Type
    distance :: Getting m a V3R -> Getting m a V3R -> Getting (ReturnMetric m) a R
    angle    :: Getting m a V3R -> Getting m a V3R -> Getting m a V3R -> Getting (ReturnMetric m) a R
    dihedral :: Getting m a V3R -> Getting m a V3R -> Getting m a V3R -> Getting m a V3R -> Getting (ReturnMetric m) a R

instance Metricable (First V3R) where
    type ReturnMetric (First V3R) = First R
    distance x y     _ aa = Const . First $ G.distance <$> (aa ^? x) <*> (aa ^? y)
    angle    x y z   _ aa = Const . First $ G.angle <$> ((-) <$> aa ^? x <*> aa ^? y) <*> ((-) <$> aa ^? z <*> aa ^? y)
    dihedral x y z w _ aa = Const . First $ G.dihedral <$> (aa ^? x) <*> (aa ^? y) <*> (aa ^? z) <*> (aa ^? w)

instance Metricable V3R where
    type ReturnMetric V3R = R
    distance x y _ aa = Const $ G.distance (aa ^. x) (aa ^. y)
    angle x y z _ aa = Const $ G.angle (aa ^. x - aa ^. y) (aa ^. z - aa ^. y)
    dihedral x y z w _ aa = Const $ G.dihedral (aa ^. x) (aa ^. y) (aa ^. z) (aa ^. w)
