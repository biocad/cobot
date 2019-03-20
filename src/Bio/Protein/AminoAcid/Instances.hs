{-# OPTIONS_GHC -fno-warn-orphans  #-}
{-# LANGUAGE TypeSynonymInstances  #-}

module Bio.Protein.AminoAcid.Instances where

import           Data.Coerce
import           Control.Lens
import           Bio.Utils.Monomer
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
    radical = ca' . environment . coerced

instance HasRadical Radical where
    type RadicalType Radical a = Radical a
    radical = ca' . environment

instance HasRadical CG where
    type RadicalType CG a = a
    radical = ca' . environment . cg'

instance HasRadical Identity where
    type RadicalType Identity a = a
    radical = ca' . environment . coerced

-- | Has lens to observe radical types
--
class Functor r => HasRadicalType r where
    -- | Getter for radical type
    --
    radicalType :: (Functor f, Functor g) => Getting AA (AminoAcid f (Env r) g a) AA

instance HasRadicalType (Const AA) where
    radicalType = ca' . environment . coerced

instance HasRadicalType CG where
    radicalType = ca' . environment . radical' . coerced

instance HasRadicalType Radical where
    radicalType = ca' . environment . to rad2rad

-- | Has lens to observe, set and modify ca_ atom
--
class Functor r => HasCA r where
    -- | Lens for ca_ atom
    --
    ca :: (Functor f, Functor g) => Lens' (AminoAcid f r g a) a

instance HasCA Identity where
    ca = ca' . coerced

instance Functor f => HasCA (Env f) where
    ca = ca' . atom'

-- | Has lens to observe, set and modify c_ atom
--
class Functor r => HasC r where
    -- | Lens for c_ atom
    --
    c :: (Functor f, Functor g) => Lens' (AminoAcid f g r a) a

instance HasC Identity where
    c = c' . coerced

instance Functor f => HasC (Env f) where
    c = c' . atom'

-- | Has lens to observe, set and modify o_ atom
--
class Functor r => HasO r where
    -- | Lens for o_ atom
    --
    o :: (Functor f, Functor g) => Lens' (AminoAcid f g (Env r) a) a

instance HasO Identity where
    o = c' . environment . coerced

instance HasO OXT where
    o = c' . environment . o'

-- | Has lens to observe, set and modify OXT atom
--
class Functor r => HasOXT r where
    -- | Lens for OXT atom
    --
    oxt :: (Functor f, Functor g) => Lens' (AminoAcid f g (Env r) a) a

instance HasOXT OXT where
    oxt = c' . environment . oxt'

-- | Has lens to observe, set and modify n_ atom
--
class Functor r => HasN r where
    -- | Lens for n_ atom
    --
    n :: (Functor f, Functor g) => Lens' (AminoAcid r f g a) a

instance HasN Identity where
    n = n' . coerced

instance Functor f => HasN (Env f) where
    n = n' . atom'

-- | Lens to get atom from some enviroment
--
class Functor f => HasAtom f where
    -- | Lens for exact atom get
    --
    atom :: Lens' (f a) a

instance HasAtom Identity where
    atom = coerced

instance Functor r => HasAtom (Env r) where
    atom = atom'

-- | Lens to get hydrogens from hydrated atom
--
hydrogens :: Lens' (Env [] a) [a]
hydrogens = environment

-- | Lens to get Symbol from every suitable amino acid
--
instance (Functor nr, HasRadicalType car, Functor cr) => Symbol (AminoAcid nr (Env car) cr a) where
    symbol = symbol . (^. radicalType)
