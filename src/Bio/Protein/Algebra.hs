module Bio.Protein.Algebra where
    -- ( phi
    -- , psi
    -- , omega
    -- , chi
    -- ) where

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

-- | Measure and rotate Chi (1, 2, 3, 4, 5) dihedral angles
--
chi :: forall nr cr h m.(HasN nr, Functor cr, HasAtom h, m ~ AminoAcid nr (Env Radical) cr (h V3R)) => Int -> Traversal' m R
chi i = lens getChi setChi . traverse
  where
    checkI :: Bool
    checkI = i > 0 && i < 6

    getChi :: m -> Maybe R
    getChi | checkI    = (^? dihedral @(First V3R) (chiP i) (chiP (i + 1)) (chiP (i + 2)) (chiP (i + 3)))
           | otherwise = const Nothing

    setChi :: m -> Maybe R -> m
    setChi m Nothing              = m
    setChi m (Just d) | checkI    = safeSetChi m d
                      | otherwise = m

    safeSetChi :: m -> R -> m
    safeSetChi m d = case getChi m of
                       Nothing  -> m
                       Just cud ->
                         let ray = Ray (m ^?! chiP (i + 1)) (normalize $ m ^?! chiP (i + 2) - m ^?! chiP (i + 1))
                             rot = rotateR ray (cud - d) :: V3R -> V3R
                         in  rotateRadical i rot m

    rotateRadical :: Int -> (V3R -> V3R) -> m -> m
    rotateRadical j rot m | j == 1 && m ^. radicalType /= PRO = rr $ m & radical %~ fmap (fmap rot)
                          -- Chi 2
                          | j == 2 && m ^. radicalType == ASP =      m & radical . cg  %~ fmap rot
                                                                       & radical . od1 %~ fmap rot
                                                                       & radical . od2 %~ fmap rot
                          | j == 2 && m ^. radicalType == PHE =      m & radical . cg  %~ fmap rot
                                                                       & radical . cd1 %~ fmap rot
                                                                       & radical . cd2 %~ fmap rot
                                                                       & radical . ce1 %~ fmap rot
                                                                       & radical . ce2 %~ fmap rot
                                                                       & radical . cz  %~ fmap rot
                          | j == 2 && m ^. radicalType == HIS =      m & radical . cg  %~ fmap rot
                                                                       & radical . nd1 %~ fmap rot
                                                                       & radical . cd2 %~ fmap rot
                                                                       & radical . ce1 %~ fmap rot
                                                                       & radical . ne2 %~ fmap rot
                          | j == 2 && m ^. radicalType == ILE =      m & radical . cg1 %~ fmap rot
                                                                       & radical . cg2 %~ fmap rot
                                                                       & radical . cd1 %~ fmap rot
                          | j == 2 && m ^. radicalType == LEU =      m & radical . cg  %~ fmap rot
                                                                       & radical . cd1 %~ fmap rot
                                                                       & radical . cd2 %~ fmap rot
                          | j == 2 && m ^. radicalType == ASN =      m & radical . cg  %~ fmap rot
                                                                       & radical . od1 %~ fmap rot
                                                                       & radical . nd2 %~ fmap rot
                          | j == 2 && m ^. radicalType == TRP =      m & radical . cg  %~ fmap rot
                                                                       & radical . cd1 %~ fmap rot
                                                                       & radical . cd2 %~ fmap rot
                                                                       & radical . ne1 %~ fmap rot
                                                                       & radical . ce2 %~ fmap rot
                                                                       & radical . ce3 %~ fmap rot
                                                                       & radical . cz2 %~ fmap rot
                                                                       & radical . cz3 %~ fmap rot
                                                                       & radical . ch2 %~ fmap rot
                          | j == 2 && m ^. radicalType == TYR =      m & radical . cg  %~ fmap rot
                                                                       & radical . cd1 %~ fmap rot
                                                                       & radical . cd2 %~ fmap rot
                                                                       & radical . ce1 %~ fmap rot
                                                                       & radical . ce2 %~ fmap rot
                                                                       & radical . cz  %~ fmap rot
                                                                       & radical . ch2 %~ fmap rot
                          | j == 2 && m ^. radicalType == GLU = rr $ m & radical . cg  %~ fmap rot
                          | j == 2 && m ^. radicalType == MET = rr $ m & radical . cg  %~ fmap rot
                          | j == 2 && m ^. radicalType == GLN = rr $ m & radical . cg  %~ fmap rot
                          | j == 2 && m ^. radicalType == LYS = rr $ m & radical . cg  %~ fmap rot
                          | j == 2 && m ^. radicalType == ARG = rr $ m & radical . cg  %~ fmap rot
                          -- Chi 3
                          | j == 3 && m ^. radicalType == GLU =      m & radical . cd  %~ fmap rot
                                                                       & radical . oe1 %~ fmap rot
                                                                       & radical . oe2 %~ fmap rot
                          | j == 3 && m ^. radicalType == MET =      m & radical . sd  %~ fmap rot
                                                                       & radical . ce  %~ fmap rot
                          | j == 3 && m ^. radicalType == GLN =      m & radical . cd  %~ fmap rot
                                                                       & radical . oe1 %~ fmap rot
                                                                       & radical . ne2 %~ fmap rot
                          | j == 3 && m ^. radicalType == LYS = rr $ m & radical . cd  %~ fmap rot
                          | j == 3 && m ^. radicalType == ARG = rr $ m & radical . cd  %~ fmap rot
                          -- Chi 4
                          | j == 4 && m ^. radicalType == LYS =      m & radical . ce  %~ fmap rot
                                                                       & radical . nz  %~ fmap rot
                          | j == 4 && m ^. radicalType == ARG = rr $ m & radical . ne  %~ fmap rot
                          -- Chi 5
                          | j == 5 && m ^. radicalType == ARG =      m & radical . cz  %~ fmap rot
                                                                       & radical . nh1 %~ fmap rot
                                                                       & radical . nh2 %~ fmap rot
                          | otherwise                         =      m
      where
        rr  = rotateRadical (j + 1) rot

-- Helper functions

-- | Chi angle point (one to eight)
--   Points 1, 2, 3 and 4 — Chi 1
--   Points 2, 3, 4 and 5 - Chi 2
--   Points 3, 4, 5 and 6 - Chi 3
--   Points 4, 5, 6 and 7 - Chi 4
--   Points 5, 6, 7 and 8 - Chi 5
--
chiP :: forall a nr cr h m.(HasN nr, Functor cr, HasAtom h, m ~ AminoAcid nr (Env Radical) cr (h a)) => Int -> Traversal' m a
chiP i = lens getChiP setChiP . traverse
  where
    checkI :: Bool
    checkI = i > 0 && i < 9

    chiPL :: Int -> AA -> Traversal' m a
    chiPL 1 _              = n . atom
    chiPL 2 _              = ca . atom
    chiPL 3 _              = radical . cb . atom
    chiPL 4 aa | aa == CYS = radical . sg  . atom
               | aa == ILE = radical . cg1 . atom
               | aa == SER = radical . og  . atom
               | aa == THR = radical . og1 . atom
               | otherwise = radical . cg  . atom
    chiPL 5 aa | aa == ASN = radical . od1 . atom
               | aa == ASP = radical . od1 . atom
               | aa == HIS = radical . nd1 . atom
               | aa == MET = radical . sd  . atom
               | aa == LEU = radical . cd1 . atom
               | aa == PHE = radical . cd1 . atom
               | aa == TRP = radical . cd1 . atom
               | aa == TYR = radical . cd1 . atom
               | otherwise = radical . cd  . atom
    chiPL 6 aa | aa == ARG = radical . ne  . atom
               | aa == GLN = radical . oe1 . atom
               | aa == GLU = radical . oe1 . atom
               | otherwise = radical . ce  . atom
    chiPL 7 aa | aa == LYS = radical . nz  . atom
               | otherwise = radical . cz  . atom
    chiPL 8 _              = radical . nh1 . atom
    chiPL _ _              = error "You cannot be here, as Chi dihedrals involves only 8 points"

    getChiP :: m -> Maybe a
    getChiP m | checkI    = m ^? chiPL i (m ^. radicalType)
              | otherwise = Nothing

    setChiP :: m -> Maybe a -> m
    setChiP m Nothing  = m
    setChiP m (Just v) | checkI    = over (chiPL i (m ^. radicalType)) (const v) m
                       | otherwise = m

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