module Bio.Utils.Geometry
    ( R
    , V3R
    , AffineTransformable(..)
    , Epsilon (..)
    , cross, dot
    , norm , normalize
    , distance, angle, dihedral
    ) where

import           Linear.V3                      ( V3
                                                , cross
                                                )
import           Linear.Epsilon                 ( Epsilon (..) )
import           Linear.Metric                  ( dot
                                                , norm
                                                , normalize
                                                , distance
                                                )
import qualified Linear.Quaternion             as Q
                                                ( rotate
                                                , axisAngle
                                                )

-- | Default floating point type, switch here to move to Doubles
--
type R = Float

-- | Defalut type of 3D vectors
--
type V3R = V3 R

-- | Affine transformations for vectors and sets of vectors
--
class AffineTransformable a where
    -- | Rotate an object around the vector by some angle
    --
    rotate    :: V3R -> R -> a -> a
    -- | Translocate an object by some vectors
    --
    translate :: V3R -> a -> a

-- | We can apply affine transformations to vectors
--
instance AffineTransformable V3R where
    rotate v a = Q.rotate (Q.axisAngle v a)
    translate v = (v +)

-- | If we have any collection of vectors, than we can transform it too
--
instance Functor f => AffineTransformable (f V3R) where
    rotate v a = fmap (rotate v a)
    translate v = fmap (translate v)

-- | Measure angle between vectors
--
angle :: V3R -> V3R -> R
angle a b = atan2 (norm (a `cross` b)) (a `dot` b)

-- | Measure dihedral between four points
--
dihedral :: V3R -> V3R -> V3R -> V3R -> R
dihedral x y z w = let b1 = y - x
                       b2 = z - y
                       b3 = w - z
                       n1 = normalize $ b1 `cross` b2
                       n2 = normalize $ b2 `cross` b3
                       m1 = n1 `cross` normalize b2
                   in  atan2 (m1 `dot` n2) (n1 `dot` n2)