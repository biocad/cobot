{-# LANGUAGE TypeFamilies          #-}
{-# LANGUAGE AllowAmbiguousTypes   #-}
{-# LANGUAGE TypeSynonymInstances  #-}
{-# LANGUAGE FlexibleInstances     #-}

module Bio.Protein.AminoAcid.Instances where

import           Data.Coerce
import           Control.Lens
import           Bio.Internal.Structure
import           Bio.Protein.AminoAcid.Type

-- | Single object can be created
--
class Createable a where
    type Create a :: *
    -- | Function to create single object
    --
    create :: Create a

instance Createable (BB a) where
    type Create (BB a) = a -> a -> a -> BB a
    create n_ ca_ c_ = coerce <$> AminoAcid (pure n_) (pure ca_) (pure c_)

instance Createable (BBCA a) where
    type Create (BBCA a) = a -> BBCA a
    create ca_ = coerce <$> AminoAcid (Const ()) (pure ca_) (Const ())

instance Createable (BBT a) where
    type Create (BBT a) = a -> a -> a -> AA -> BBT a
    create n_ ca_ c_ aa = coerce <$> AminoAcid (pure n_) (Env ca_ (Const aa)) (pure c_)

instance Createable (BBCAT a) where
    type Create (BBCAT a) = a -> AA -> BBCAT a
    create ca_ aa = coerce <$> AminoAcid (Const ()) (Env ca_ (Const aa)) (Const ())

instance Createable (BBCG a) where
    type Create (BBCG a) = a -> a -> a -> a -> AA -> BBCG a
    create n_ ca_ c_ cg_ aa = coerce <$> AminoAcid (pure n_) (Env ca_ (CG cg_ aa)) (pure c_)

instance Createable (BBO a) where
    type Create (BBO a) = a -> a -> a -> a -> BBO a
    create n_ ca_ c_ o_ = coerce <$> AminoAcid (pure n_) (pure ca_) (Env c_ (pure o_))

instance Createable (BBOT a) where
    type Create (BBOT a) = a -> a -> a -> a -> AA -> BBOT a
    create n_ ca_ c_ o_ aa = coerce <$> AminoAcid (pure n_) (Env ca_ (Const aa)) (Env c_ (pure o_))

instance Createable (BBOCG a) where
    type Create (BBOCG a) = a -> a -> a -> a -> a -> AA -> BBOCG a
    create n_ ca_ c_ o_ cg_ aa = coerce <$> AminoAcid (pure n_) (Env ca_ (CG cg_ aa)) (Env c_ (pure o_))

instance Createable (BBOR a) where
    type Create (BBOR a) = a -> a -> a -> a -> Radical a -> BBOR a
    create n_ ca_ c_ o_ r = coerce <$> AminoAcid (pure n_) (Env ca_ r) (Env c_ (pure o_))

instance Createable (BBOXTR a) where
    type Create (BBOXTR a) = a -> a -> a -> a -> a -> Radical a -> BBOXTR a
    create n_ ca_ c_ o_ oxt_ r = coerce <$> AminoAcid (pure n_) (Env ca_ r) (Env c_ (OXT o_ oxt_))

instance Createable (BBORH a) where
    type Create (BBORH a) = a -> a -> a -> a -> Radical a -> BBORH a
    create n_ ca_ c_ o_ r = flip Env [] <$> AminoAcid (pure n_) (Env ca_ r) (Env c_ (pure o_))

instance Createable (BBOXTRH a) where
    type Create (BBOXTRH a) = a -> a -> a -> a -> a -> Radical a -> BBOXTRH a
    create n_ ca_ c_ o_ oxt_ r = flip Env [] <$> AminoAcid (pure n_) (Env ca_ r) (Env c_ (OXT o_ oxt_))

-- | Has lens to observe, set and modify radicals
--
class Functor r => HasRadical r where
    type RadicalType r a :: *
    -- | Lens for radical atom or group
    --
    radical :: (Functor f, Functor g) => Lens' (AminoAcid f (Env r) g a) (RadicalType r a)

instance HasRadical (Const x) where
    type RadicalType (Const x) a = x
    radical = lens (getConst . (^. ca' . environment)) (\aa x -> set (ca' . environment) (Const x) aa)

instance HasRadical Radical where
    type RadicalType Radical a = Radical a
    radical = lens (^. ca' . environment) (\aa x -> set (ca' . environment) x aa)

instance HasRadical CG where
    type RadicalType CG a = a
    radical = lens (^. ca' . environment . cg') (\aa x -> set (ca' . environment . cg') x aa)

instance HasRadical Identity where
    type RadicalType Identity a = a
    radical = lens (runIdentity . (^. ca' . environment)) (\aa x -> set (ca' . environment) (Identity x) aa)

-- | Has lens to observe radical types
--
class Functor r => HasRadicalType r where
    -- | Getter for radical type
    --
    radicalType :: (Functor f, Functor g) => Getting AA (AminoAcid f (Env r) g a) AA

instance HasRadicalType (Const AA) where
    radicalType _ aa = Const $ getConst (aa ^. ca' . environment)

instance HasRadicalType CG where
    radicalType _ aa = Const (aa ^. ca' . environment . radical')

instance HasRadicalType Radical where
    radicalType _ aa = Const (rad2rad $ aa ^. ca' . environment)
      where
        rad2rad :: Radical a -> AA
        rad2rad Alanine{}       = ALA
        rad2rad Cysteine{}      = CYS
        rad2rad AsparticAcid{}  = ASP
        rad2rad GlutamicAcid{}  = GLU
        rad2rad Phenylalanine{} = PHE
        rad2rad Glycine         = GLY
        rad2rad Histidine{}     = HIS
        rad2rad Isoleucine{}    = ISO
        rad2rad Lysine{}        = LYS
        rad2rad Leucine{}       = LEU
        rad2rad Methionine{}    = MET
        rad2rad Asparagine{}    = ASN
        rad2rad Proline{}       = PRO
        rad2rad Glutamine{}     = GLN
        rad2rad Arginine{}      = ARG
        rad2rad Serine{}        = SER
        rad2rad Threonine{}     = THR
        rad2rad Valine{}        = VAL
        rad2rad Tryptophan{}    = TRP
        rad2rad Tyrosine{}      = TYR

-- | Has lens to observe, set and modify ca_ atom
--
class Functor r => HasCA r where
    -- | Lens for ca_ atom
    --
    ca :: (Functor f, Functor g) => Lens' (AminoAcid f r g a) a

instance HasCA Identity where
    ca = lens (runIdentity . (^. ca')) (\aa x -> set ca' (Identity x) aa)

instance Functor f => HasCA (Env f) where
    ca = lens (^. ca' . atom') (\aa x -> set (ca' . atom') x aa)

-- | Has lens to observe, set and modify c_ atom
--
class Functor r => HasC r where
    -- | Lens for c_ atom
    --
    c :: (Functor f, Functor g) => Lens' (AminoAcid f g r a) a

instance HasC Identity where
    c = lens (runIdentity . (^. c')) (\aa x -> set c' (Identity x) aa)

instance Functor f => HasC (Env f) where
    c = lens (^. c' . atom') (\aa x -> set (c' . atom') x aa)

-- | Has lens to observe, set and modify o_ atom
--
class Functor r => HasO r where
    -- | Lens for o_ atom
    --
    o :: (Functor f, Functor g) => Lens' (AminoAcid f g (Env r) a) a

instance HasO Identity where
    o = lens (runIdentity  . (^. c' . environment)) (\aa x -> set (c' . environment) (Identity x) aa)

instance HasO OXT where
    o = lens (^. c' . environment . o') (\aa x -> set (c' . environment . o') x aa)

-- | Has lens to observe, set and modify OXT atom
--
class Functor r => HasOXT r where
    -- | Lens for OXT atom
    --
    oxt :: (Functor f, Functor g) => Lens' (AminoAcid f g (Env r) a) a

instance HasOXT OXT where
    oxt = lens (^. c' . environment . oxt') (\aa x -> set (c' . environment . oxt') x aa)

-- | Has lens to observe, set and modify n_ atom
--
class Functor r => HasN r where
    -- | Lens for n_ atom
    --
    n :: (Functor f, Functor g) => Lens' (AminoAcid r f g a) a

instance HasN Identity where
    n = lens (runIdentity . (^. n')) (\aa x -> set n' (Identity x) aa)

instance Functor f => HasN (Env f) where
    n = lens (^. n' . atom') (\aa x -> set (n' . atom') x aa)

-- | Lens to get atom from some enviroment
--
class Functor f => HasAtom f where
    atom :: Lens' (f a) a

instance HasAtom Identity where
    atom = lens runIdentity (\_ x -> Identity x)

instance Functor r => HasAtom (Env r) where
    atom = lens (^. atom') (\env x -> set atom' x env)

-- | Lens to get hydrogens from hydrated atom
--
hydrogens :: Lens' (Env [] a) [a]
hydrogens = environment