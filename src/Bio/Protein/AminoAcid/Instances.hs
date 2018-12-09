{-# LANGUAGE TypeFamilies          #-}

module Bio.Protein.AminoAcid.Instances where

import           Control.Lens
import           Bio.Internal.Structure
import           Bio.Protein.AminoAcid.Type

-- |Has lens to observe, set and modify radicals
class Functor r => HasRadical r where
    type RadicalType r a :: *
    -- |Lens for radical atom or group
    radical :: (Functor f, Functor g) => Lens' (AminoAcid f (Env r) g a) (RadicalType r a)

instance HasRadical (Const x) where
    type RadicalType (Const x) a = x
    radical = lens (getConst . (^. ca' . environment)) (\aa x -> set (ca' . environment) (Const x) aa)

instance HasRadical Radical where
    type RadicalType Radical a = Radical a
    radical = lens (^. ca' . environment) (\aa x -> set (ca' . environment) x aa)

instance HasRadical Identity where
    type RadicalType Identity a = a
    radical = lens (runIdentity . (^. ca' . environment)) (\aa x -> set (ca' . environment) (Identity x) aa)

-- |Has lens to observe, set and modify CA atom
class Functor r => HasCA r where
    -- |Lens for CA atom
    ca :: (Functor f, Functor g) => Lens' (AminoAcid f r g a) a

instance HasCA Identity where
    ca = lens (runIdentity . (^. ca')) (\aa x -> set ca' (Identity x) aa)

instance Functor f => HasCA (Env f) where
    ca = lens (^. ca' . atom) (\aa x -> set (ca' . atom) x aa)

-- |Has lens to observe, set and modify C atom
class Functor r => HasC r where
    -- |Lens for C atom
    c :: (Functor f, Functor g) => Lens' (AminoAcid f g r a) a

instance HasC Identity where
    c = lens (runIdentity . (^. c')) (\aa x -> set c' (Identity x) aa)

instance Functor f => HasC (Env f) where
    c = lens (^. c' . atom) (\aa x -> set (c' . atom) x aa)

-- |Has lens to observe, set and modify O atom
class Functor r => HasO r where
    -- |Lens for O atom
    o :: (Functor f, Functor g) => Lens' (AminoAcid f g (Env r) a) a

instance HasO Identity where
    o = lens (runIdentity  . (^. c' . environment)) (\aa x -> set (c' . environment) (Identity x) aa)

instance HasO OXT where
    o = lens (^. c' . environment . o') (\aa x -> set (c' . environment . o') x aa)

-- |Has lens to observe, set and modify OXT atom
class Functor r => HasOXT r where
    -- |Lens for OXT atom
    oxt :: (Functor f, Functor g) => Lens' (AminoAcid f g (Env r) a) a

instance HasOXT OXT where
    oxt = lens (^. c' . environment . oxt') (\aa x -> set (c' . environment . oxt') x aa)

-- |Has lens to observe, set and modify N atom
class Functor r => HasN r where
    -- |Lens for N atom
    n :: (Functor f, Functor g) => Lens' (AminoAcid r f g a) a

instance HasN Identity where
    n = lens (runIdentity . (^. n')) (\aa x -> set n' (Identity x) aa)

instance Functor f => HasN (Env f) where
    n = lens (^. n' . atom) (\aa x -> set (n' . atom) x aa)