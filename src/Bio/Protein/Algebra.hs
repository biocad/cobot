module Bio.Protein.Algebra
    ( phi
    , psi
    , omega
    , chi
    ) where

import           Data.Monoid                    ( First (..) )
import           Control.Lens
import           Bio.Utils.Geometry             ( V3R, R, Ray (..), normalize, rotateR )

import           Bio.Protein.AminoAcid
import           Bio.Protein.Metric
import           Bio.Protein.Chain

type Dihedral m f r g h = (ChainLike m, HasN f, HasCA r, HasC g, HasAtom h, IxValue m ~ AminoAcid f r g (h V3R))

-- | Measure and rotate Psi dihedral angle
--
psi :: forall m f r g h.Dihedral m f r g h => Index m -> Traversal' m R
psi i = rcd (\rot -> (& c %~ fmap rot)) (ix i . n . atom) (ix i . ca . atom) (ix i . c . atom) (ix (succ i) . n . atom) i

-- | Measure and rotate Phi dihedral angle
--
phi :: forall m f r g h.Dihedral m f r g h => Index m -> Traversal' m R
phi i = rcd (\rot -> (& ca %~ fmap rot) . (& c %~ fmap rot)) (ix (pred i) . c . atom) (ix i . n . atom) (ix i . ca . atom) (ix i . c . atom) i

-- | Measure and rotate Omega dihedral angle
--
omega :: forall m f r g h.Dihedral m f r g h => Index m -> Traversal' m R
omega i = rcd (fmap . fmap) (ix (pred i) . ca . atom) (ix (pred i) . c . atom) (ix i . n . atom) (ix i . ca . atom) i

-- | Measure and rotate Chi (1, 2, 3, 4) dihedral angles
--
chi :: (ChainLike m, IxValue m ~ AminoAcid f (Env Radical) g (h V3R)) => Int -> Index m -> Lens' m R
chi _ _ = undefined

-- Helper functions

type ModifyFunction m = (V3R -> V3R) -> IxValue m -> IxValue m

-- | Rotate cannonical dihedral in backbone
--
rcd :: forall m f r g h.Dihedral m f r g h => ModifyFunction m {- modify function -} ->
                                              Traversal' m V3R {- first  point    -} ->
                                              Traversal' m V3R {- second point    -} ->
                                              Traversal' m V3R {- third  point    -} ->
                                              Traversal' m V3R {- fourth point    -} ->
                                              Index m          {- dihedral index  -} ->
                                              Traversal' m R
rcd mf x1 x2 x3 x4 i = lens getRCD setRCD . traverse
  where
    getRCD :: m -> Maybe R
    getRCD = (^? dihedral @(First V3R) x1 x2 x3 x4)

    setRCD :: m -> Maybe R -> m
    setRCD ar Nothing  = ar
    setRCD ar (Just d) = case getRCD ar of
                           Nothing -> ar
                           Just cud ->
                              let ray = Ray (ar ^?! x2) (normalize $ ar ^?! x3 - ar ^?! x2)
                                  rot = rotateR ray (cud - d)
                                  mfy = modify i (mf rot) . modifyAfter i (fmap (fmap rot))
                              in  mfy ar