{-# LANGUAGE DeriveFunctor   #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE DeriveAnyClass #-}

module Bio.Protein.AminoAcid.Type where

import           Control.DeepSeq        (NFData (..))
import           Control.Lens
import           Control.Monad.Identity (Identity)
import           GHC.Generics           (Generic (..))

-- | Proteinogenic amino acids
--
data AA = ALA -- A
        | CYS -- C
        | ASP -- D
        | GLU -- E
        | PHE -- F
        | GLY -- G
        | HIS -- H
        | ILE -- I
        | LYS -- K
        | LEU -- L
        | MET -- M
        | ASN -- N
        | PRO -- P
        | GLN -- Q
        | ARG -- R
        | SER -- S
        | THR -- T
        | VAL -- V
        | TRP -- W
        | TYR -- Y
  deriving (Eq, Ord, Bounded, Enum, Generic, NFData)

-- | Show full names of amino acids
--
instance Show AA where
    show ALA = "Alanine"
    show CYS = "Cysteine"
    show ASP = "AsparticAcid"
    show GLU = "GlutamicAcid"
    show PHE = "Phenylalanine"
    show GLY = "Glycine"
    show HIS = "Histidine"
    show ILE = "Isoleucine"
    show LYS = "Lysine"
    show LEU = "Leucine"
    show MET = "Methionine"
    show ASN = "Asparagine"
    show PRO = "Proline"
    show GLN = "Glutamine"
    show ARG = "Arginine"
    show SER = "Serine"
    show THR = "Threonine"
    show VAL = "Valine"
    show TRP = "Tryptophan"
    show TYR = "Tyrosine"

-- | Amino acid structure type
--
data AminoAcid nr car cr a = AminoAcid { _n'  :: nr a
                                       , _ca' :: car a
                                       , _c'  :: cr a
                                       }
  deriving (Show, Eq, Functor, Generic, NFData)

-- | Radical structure type
--
data Radical a = Alanine          --  no chi
                   { _cb  :: a    --  -CB
                   }              --
               | Cysteine         --
                   { _cb :: a    --  -CB-SG
                   , _sg :: a    --
                   }              --
               | AsparticAcid     --
                   { _cb  :: a    --  -CB-CG-OD1
                   , _cg  :: a    --       |
                   , _od1 :: a    --       OD2
                   , _od2 :: a    --
                   }              --
               | GlutamicAcid     --
                   { _cb  :: a    --  -CB-CG-CD-OE1
                   , _cg  :: a    --          |
                   , _cd  :: a    --          OE2
                   , _oe1 :: a    --
                   , _oe2 :: a    --
                   }              --
               | Phenylalanine    --
                   { _cb  :: a    --  -CB-CG-CD1-CE1
                   , _cg  :: a    --       |       |
                   , _cd1 :: a    --       CD2-CE2-CZ
                   , _cd2 :: a    --
                   , _ce1 :: a    --
                   , _ce2 :: a    --
                   , _cz  :: a    --
                   }              --
               | Glycine          --
               | Histidine        --
                   { _cb  :: a    --  -CB-CG-ND1-CE1
                   , _cg  :: a    --       |
                   , _nd1 :: a    --       CD2-NE2
                   , _cd2 :: a    --
                   , _ce1 :: a    --
                   , _ne2 :: a    --
                   }              --
               | Isoleucine       --
                   { _cb  :: a    --  -CB-CG1-CD1
                   , _cg1 :: a    --    |
                   , _cg2 :: a    --    CG2
                   , _cd1 :: a    --
                   }              --
               | Lysine           --
                   { _cb :: a    --  -CB-CG-CD-CE-NZ
                   , _cg :: a    --
                   , _cd :: a    --
                   , _ce :: a    --
                   , _nz :: a    --
                   }              --
               | Leucine          --
                   { _cb  :: a    --  -CB-CG-CD1
                   , _cg  :: a    --       |
                   , _cd1 :: a    --       CD2
                   , _cd2 :: a    --
                   }              --
               | Methionine       --
                   { _cb :: a    --  -CB-CG-SD-CE
                   , _cg :: a    --
                   , _sd :: a    --
                   , _ce :: a    --
                   }              --
               | Asparagine       --
                   { _cb  :: a    --  -CB-CG-OD1
                   , _cg  :: a    --       |
                   , _od1 :: a    --       ND2
                   , _nd2 :: a    --
                   }              --
               | Proline          --
                   { _cb :: a    --  -CB-CG-CD(-N)
                   , _cg :: a    --
                   , _cd :: a    --
                   }              --
               | Glutamine        --
                   { _cb  :: a    --  -CB-CG-CD-OE1
                   , _cg  :: a    --          |
                   , _cd  :: a    --          NE2
                   , _oe1 :: a    --
                   , _ne2 :: a    --
                   }              --
               | Arginine         --
                   { _cb  :: a    --  -CB-CG-CD-NE-CZ-NH1
                   , _cg  :: a    --                |
                   , _cd  :: a    --                NH2
                   , _ne  :: a    --
                   , _cz  :: a    --
                   , _nh1 :: a    --
                   , _nh2 :: a    --
                   }              --
               | Serine           --
                   { _cb :: a    --  -CB-OG
                   , _og :: a    --
                   }              --
               | Threonine        --
                   { _cb  :: a    --  -CB-OG1
                   , _og1 :: a    --    |
                   , _cg2 :: a    --    CG2
                   }              --
               | Valine           --
                   { _cb  :: a    --  -CB-CG1
                   , _cg1 :: a    --    |
                   , _cg2 :: a    --    CG2
                   }              --
               | Tryptophan       --
                   { _cb  :: a    --  -CB-CG-CD1-NE1
                   , _cg  :: a    --       |     |
                   , _cd1 :: a    --       CD2-CE2-CZ2
                   , _cd2 :: a    --       |        |
                   , _ne1 :: a    --       CE3-CZ3-CH2
                   , _ce2 :: a    --
                   , _ce3 :: a    --
                   , _cz2 :: a    --
                   , _cz3 :: a    --
                   , _ch2 :: a    --
                   }              --
               | Tyrosine         --
                   { _cb  :: a    --  -CB-CG-CD1-CE1
                   , _cg  :: a    --       |       |
                   , _cd1 :: a    --       CD2-CE2-CZ-OH
                   , _cd2 :: a    --
                   , _ce1 :: a    --
                   , _ce2 :: a    --
                   , _cz  :: a    --
                   , _oh  :: a    --
                   }              --
                     deriving (Show, Eq, Functor, Generic, NFData)

-- | Atom environment, e.g. hydrogens or radicals
--
data Env r a = Env { _atom'       :: a
                   , _environment :: r a
                   }
  deriving (Show, Eq, Functor, Generic, NFData)

-- | Hydrogens envrironment
--
type H a = Env [] a

-- | Oxigen and hydroxi group, connected to C-terminal of amino acid
--
data OXT a = OXT { _o'   :: a
                 , _oxt' :: a
                 }
  deriving (Show, Eq, Functor, Generic, NFData)

-- | CG atom with radical type
--
data CG a = CG { _cg'      :: a
               , _radical' :: AA
               }
  deriving (Show, Eq, Functor, Generic, NFData)

makeLenses ''AminoAcid
makeLenses ''Radical
makeLenses ''Env
makeLenses ''OXT
makeLenses ''CG

-- | BackBone
--
type BB a      = AminoAcid Identity   Identity                  Identity       (Identity a)

-- | BackBone CA-only
--
type BBCA a    = AminoAcid (Const ()) Identity                  (Const ())     (Identity a)

-- | BackBone with radical Type
--
type BBT a     = AminoAcid Identity   (Env (Const AA))          Identity       (Identity a)

-- | BackBone CA-only with radical Type
--
type BBCAT a   = AminoAcid (Const ()) (Env (Const AA))          (Const ())     (Identity a)

-- | BackBone with CG-radical
--
type BBCG a    = AminoAcid Identity   (Env CG)                  Identity       (Identity a)

-- | BackBone with Oxigen
--
type BBO a     = AminoAcid Identity    Identity                 (Env Identity) (Identity a)

-- | BackBone with Oxigen and radical Type
--
type BBOT a    = AminoAcid Identity    (Env (Const AA))         (Env Identity) (Identity a)

-- | BackBone with Oxigen and CG-radical
--
type BBOCG a   = AminoAcid Identity    (Env CG)                 (Env Identity) (Identity a)

-- | BackBone with Oxigen and Radical
--
type BBOR a    = AminoAcid Identity    (Env Radical)            (Env Identity) (Identity a)

-- | BackBone with Oxigen, oXigen Two and Radical
--
type BBOXTR a  = AminoAcid Identity    (Env Radical)            (Env OXT)      (Identity a)

-- | BackBone with Oxigen, Radical and Hydrogens
--
type BBORH a   = AminoAcid Identity    (Env Radical)            (Env Identity) (H a)

-- | BackBone with Oxigen, oXigen Two, Radical and Hydrogens
--
type BBOXTRH a = AminoAcid Identity    (Env Radical)            (Env OXT)      (H a)
--
-- | Convert radical to radical name
--
rad2rad :: Radical a -> AA
rad2rad Alanine{}       = ALA
rad2rad Cysteine{}      = CYS
rad2rad AsparticAcid{}  = ASP
rad2rad GlutamicAcid{}  = GLU
rad2rad Phenylalanine{} = PHE
rad2rad Glycine         = GLY
rad2rad Histidine{}     = HIS
rad2rad Isoleucine{}    = ILE
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
