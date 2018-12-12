{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeFamilies     #-}
{-# LANGUAGE RankNTypes       #-}
module Bio.Protein.Algebra
    ( phi
    , psi
    , omega
    , chi
    ) where

import           Control.Lens
import           Bio.Utils.Lens                 ( idx )
import           Bio.Utils.Geometry             ( V3R, R )

import           Bio.Protein.AminoAcid
import           Bio.Protein.Metric

phi :: (Ixed m, Num (Index m), HasN f, HasCA r, HasC g, HasAtom h, IxValue m ~ AminoAcid f r g (h V3R)) => Index m -> Lens' m R
phi i = lens (^. dihedral (idx i . n . atom) (idx i . ca . atom) (idx i . c . atom) (idx (i + 1) . n . atom))
             (\ar _ -> ar)

psi :: (Ixed m, Num (Index m), HasN f, HasCA r, HasC g, HasAtom h, IxValue m ~ AminoAcid f r g (h V3R)) => Index m -> Lens' m R
psi i = lens (^. dihedral (idx (i - 1) . c . atom) (idx i . n . atom) (idx i . ca . atom) (idx i . c . atom))
             (\ar _ -> ar)

omega :: (Ixed m, Num (Index m), HasN f, HasCA r, HasC g, HasAtom h, IxValue m ~ AminoAcid f r g (h V3R)) => Index m -> Lens' m R
omega i = lens (^. dihedral (idx (i - 1) . ca . atom) (idx (i - 1) . c . atom) (idx i . n . atom) (idx i . ca . atom))
               (\ar _ -> ar)

chi :: (Ixed m, Num (Index m), IxValue m ~ AminoAcid f (Env Radical) g (h V3R)) => Int -> Index m -> Lens' m R
chi _ _ = undefined
