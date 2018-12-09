{-# LANGUAGE DeriveFunctor         #-}
{-# LANGUAGE TemplateHaskell       #-}
{-# LANGUAGE FlexibleInstances     #-}

module Bio.Internal.Structure where

import           Control.Lens
import           Linear.V3                      ( V3 )
import qualified Linear.Quaternion             as Q

-- |Atom environment, e.g. hydrogens or radicals
data Env r a = Env { _atom        :: a
                   , _environment :: r a
                   }
  deriving (Show, Functor)

makeLenses ''Env

-- |Hydrogens envrironment
type H a  = Env [] a

-- |Default floating point type, switch here to move to Doubles
type R = Float

-- |Defalut type of 3D vectors
type V3R = V3 R

-- |Affine transformations for vectors and sets of vectors
class AffineTransformable a where
  -- |Rotate an object around the vector by some angle
  rotate    :: V3R -> R -> a -> a
  -- |Translocate an object by some vectors
  translate :: V3R -> a -> a

-- |We can apply affine transformations to vectors
instance AffineTransformable V3R where
  rotate v a = Q.rotate (Q.axisAngle v a)
  translate v = (v +)

-- |If we have any collection of vectors, than we can transform it too
instance Functor f => AffineTransformable (f V3R) where
  rotate v a = fmap (rotate v a)
  translate v = fmap (translate v)

-- |Measurements between some objects in space
class Measureable a where
  dist :: a -> a -> R
  angle :: a -> a -> a -> R