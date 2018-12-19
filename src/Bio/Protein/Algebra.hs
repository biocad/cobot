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
psi :: forall m f r g h.Dihedral m f r g h => Index m -> Traversal' m R
psi i = lens getPsi setPsi . traverse
  where
    getPsi :: m -> Maybe R
    getPsi = (^? dihedral @(First V3R) (ix i . n . atom) (ix i . ca . atom) (ix i . c . atom) (ix (succ i) . n . atom))

    setPsi :: m -> Maybe R -> m
    setPsi ar Nothing  = ar
    setPsi ar (Just d) = case getPsi ar of
                           Nothing -> ar
                           Just cud -> 
                             let ang = cud - d   -- rotation angle (use cud - d instead of right d - cud, as we use different way to measure dihedrals)
                                 ori = ar ^?! ix i . ca . atom
                                 dir = normalize $ ar ^?! ix i . c . atom - ar ^?! ix i . ca . atom
                                 rot = rotateR (Ray ori dir) ang
                                 mfy = modify i (& c %~ fmap rot) . modifyAfter i (fmap (fmap rot))
                             in  mfy ar

-- | Measure and rotate Phi dihedral angle
phi :: forall m f r g h.Dihedral m f r g h => Index m -> Traversal' m R
phi i = lens getPhi setPhi . traverse
  where
    getPhi :: m -> Maybe R
    getPhi = (^? dihedral @(First V3R) (ix (pred i) . c . atom) (ix i . n . atom) (ix i . ca . atom) (ix i . c . atom))

    setPhi :: m -> Maybe R -> m
    setPhi ar Nothing  = ar
    setPhi ar (Just d) = case getPhi ar of
                           Nothing -> ar
                           Just cud ->
                             let ang = cud - d   -- rotation angle (use cud - d instead of right d - cud, as we use different way to measure dihedrals)
                                 ori = ar ^?! ix i . n . atom
                                 dir = normalize $ ar ^?! ix i . ca . atom - ar ^?! ix i . n . atom
                                 rot = rotateR (Ray ori dir) ang
                                 mfy = modify i (& ca %~ fmap rot) . modify i (& c %~ fmap rot) . modifyAfter i (fmap (fmap rot))
                             in  mfy ar

-- | Measure and rotate Omega dihedral angle
omega :: forall m f r g h.Dihedral m f r g h => Index m -> Traversal' m R
omega i = lens getOmega setOmega . traverse
  where
    getOmega :: m -> Maybe R
    getOmega = (^? dihedral @(First V3R) (ix (pred i) . ca . atom) (ix (pred i) . c . atom) (ix i . n . atom) (ix i . ca . atom))

    setOmega :: m -> Maybe R -> m
    setOmega ar Nothing  = ar
    setOmega ar (Just d) = case getOmega ar of
                             Nothing -> ar
                             Just cud ->
                               let ang = cud - d     -- rotation angle (use cud - d instead of right d - cud, as we use different way to measure dihedrals)
                                   ori = ar ^?! ix (pred i) . c . atom
                                   dir = normalize $ ar ^?! ix i . n . atom - ar ^?! ix (pred i) . c . atom
                                   rot = rotateR (Ray ori dir) ang
                                   mfy = modifyAfter (pred i) (fmap (fmap rot))
                               in  mfy ar

chi :: (ChainLike m, IxValue m ~ AminoAcid f (Env Radical) g (h V3R)) => Int -> Index m -> Lens' m R
chi _ _ = undefined
