module Bio.Protein.Algebra
    ( phi
    , psi
    , omega
    , chi
    ) where

import           Control.Lens
import           Bio.Utils.Lens                 ( idx )
import           Bio.Utils.Geometry             ( V3R, R, normalize, rotate )

import           Bio.Protein.AminoAcid
import           Bio.Protein.Metric
import           Bio.Protein.Chain

type Dihedral m f r g h = (ChainLike m, HasN f, HasCA r, HasC g, HasAtom h, IxValue m ~ AminoAcid f r g (h V3R))

phi :: forall m f r g h.Dihedral m f r g h => Index m -> Lens' m R
phi i = lens getPhi setPhi
  where
    getPhi = (^. dihedral (idx i . n . atom) (idx i . ca . atom) (idx i . c . atom) (idx (succ i) . n . atom))

    setPhi ar d = let cud = getPhi ar -- current dihedral
                      ang = d - cud   -- rotation angle
                      axe = normalize $ ar ^. idx i . c . atom - ar ^. idx i . ca . atom
                      rot = rotate axe ang
                      mfy = modify i (& c %~ fmap rot) . modifyAfter i (fmap (fmap rot))
                  in  mfy ar

psi :: forall m f r g h.Dihedral m f r g h => Index m -> Lens' m R
psi i = lens getPsi setPsi
  where
    getPsi = (^. dihedral (idx (pred i) . c . atom) (idx i . n . atom) (idx i . ca . atom) (idx i . c . atom))

    setPsi ar d = let cud = getPsi ar -- current dihedral
                      ang = d - cud   -- rotation angle
                      axe = normalize $ ar ^. idx i . ca . atom - ar ^. idx i . n . atom
                      rot = rotate axe ang
                      mfy = modify i (& ca %~ fmap rot) . modify i (& c %~ fmap rot) . modifyAfter i (fmap (fmap rot))
                  in  mfy ar

omega :: forall m f r g h.Dihedral m f r g h => Index m -> Lens' m R
omega i = lens getOmega setOmega
  where
    getOmega = (^. dihedral (idx (pred i) . ca . atom) (idx (pred i) . c . atom) (idx i . n . atom) (idx i . ca . atom))

    setOmega ar d = let cud = getOmega ar -- current dihedral
                        ang = d - cud   -- rotation angle
                        axe = normalize $ ar ^. idx i . n . atom - ar ^. idx (pred i) . c . atom
                        rot = rotate axe ang
                        mfy = modifyAfter (pred i) (fmap (fmap rot))
                    in  mfy ar

chi :: (ChainLike m, IxValue m ~ AminoAcid f (Env Radical) g (h V3R)) => Int -> Index m -> Lens' m R
chi _ _ = undefined
