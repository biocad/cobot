{-# LANGUAGE DeriveFunctor         #-}
{-# LANGUAGE TemplateHaskell       #-}

module Bio.Protein.AminoAcid.Type where

import           Control.Lens
import           Control.Monad.Identity         ( Identity )
import           Bio.Utils.Monomer

-- | Proteinogenic amino acids
--
data AA = ALA -- A
        | CYS -- C
        | ASP -- D
        | GLU -- E
        | PHE -- F
        | GLY -- G
        | HIS -- H
        | ISO -- I
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
  deriving (Eq, Ord, Bounded, Enum)

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
    show ISO = "Isoleucine"
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

-- | Show one symbol encoding
--
instance Symbol AA where
    symbol ALA = 'A'
    symbol CYS = 'C'
    symbol ASP = 'D'
    symbol GLU = 'E'
    symbol PHE = 'F'
    symbol GLY = 'G'
    symbol HIS = 'H'
    symbol ISO = 'I'
    symbol LYS = 'K'
    symbol LEU = 'L'
    symbol MET = 'M'
    symbol ASN = 'N'
    symbol PRO = 'P'
    symbol GLN = 'Q'
    symbol ARG = 'R'
    symbol SER = 'S'
    symbol THR = 'T'
    symbol VAL = 'V'
    symbol TRP = 'W'
    symbol TYR = 'Y'

-- | Show three symbols encoding
--
instance ThreeSymbols AA where
    threeSymbols ALA = "ALA"
    threeSymbols CYS = "CYS"
    threeSymbols ASP = "ASP"
    threeSymbols GLU = "GLU"
    threeSymbols PHE = "PHE"
    threeSymbols GLY = "GLY"
    threeSymbols HIS = "HIS"
    threeSymbols ISO = "ISO"
    threeSymbols LYS = "LYS"
    threeSymbols LEU = "LEU"
    threeSymbols MET = "MET"
    threeSymbols ASN = "ASN"
    threeSymbols PRO = "PRO"
    threeSymbols GLN = "GLN"
    threeSymbols ARG = "ARG"
    threeSymbols SER = "SER"
    threeSymbols THR = "THR"
    threeSymbols VAL = "VAL"
    threeSymbols TRP = "TRP"
    threeSymbols TYR = "TYR"

-- | Amino acid structure type
--
data AminoAcid nr car cr a = AminoAcid { _n'  :: nr a
                                       , _ca' :: car a
                                       , _c'  :: cr a
                                       }
  deriving (Show, Eq, Functor)

-- | Radical structure type
--
data Radical a = Alanine
                   { _cb  :: a
                   }
               | Cysteine
                   { _cb  :: a
                   , _sg  :: a
                   }
               | AsparticAcid
                   { _cb  :: a
                   , _cg  :: a
                   , _od1 :: a
                   , _od2 :: a
                   }
               | GlutamicAcid
                   { _cb  :: a
                   , _cg  :: a
                   , _cd  :: a
                   , _oe1 :: a
                   , _oe2 :: a
                   }
               | Phenylalanine
                   { _cb  :: a
                   , _cg  :: a
                   , _cd1 :: a
                   , _cd2 :: a
                   , _ce1 :: a
                   , _ce2 :: a
                   , _cz  :: a
                   }
               | Glycine
               | Histidine
                   { _cb  :: a
                   , _cg  :: a
                   , _nd1 :: a
                   , _cd2 :: a
                   , _ce1 :: a
                   , _ne2 :: a
                   }
               | Isoleucine
                   { _cb  :: a
                   , _cg1 :: a
                   , _cg2 :: a
                   , _cd1 :: a
                   }
               | Lysine
                   { _cb  :: a
                   , _cg  :: a
                   , _cd  :: a
                   , _ce  :: a
                   , _nz  :: a
                   }
               | Leucine
                   { _cb  :: a
                   , _cg  :: a
                   , _cd1 :: a
                   , _cd2 :: a
                   }
               | Methionine
                   { _cb  :: a
                   , _cg  :: a
                   , _sd  :: a
                   , _ce  :: a
                   }
               | Asparagine
                   { _cb  :: a
                   , _cg  :: a
                   , _od1 :: a
                   , _nd2 :: a
                   }
               | Proline
                   { _cb  :: a
                   , _cg  :: a
                   , _cd  :: a
                   }
               | Glutamine
                   { _cb  :: a
                   , _cg  :: a
                   , _cd  :: a
                   , _oe1 :: a
                   , _ne2 :: a
                   }
               | Arginine
                   { _cb  :: a
                   , _cg  :: a
                   , _cd  :: a
                   , _ne  :: a
                   , _cz  :: a
                   , _nh1 :: a
                   , _nh2 :: a
                   }
               | Serine
                   { _cb  :: a
                   , _og  :: a
                   }
               | Threonine
                   { _cb  :: a
                   , _og1 :: a
                   , _cg2 :: a
                   }
               | Valine
                   { _cb  :: a
                   , _cg1 :: a
                   , _cg2 :: a
                   }
               | Tryptophan
                   { _cb  :: a
                   , _cg  :: a
                   , _cd1 :: a
                   , _cd2 :: a
                   , _ne1 :: a
                   , _ce2 :: a
                   , _ce3 :: a
                   , _cz2 :: a
                   , _cz3 :: a
                   , _ch2 :: a
                   }
               | Tyrosine
                   { _cb  :: a
                   , _cg  :: a
                   , _cd1 :: a
                   , _cd2 :: a
                   , _ce1 :: a
                   , _ce2 :: a
                   , _cz  :: a
                   , _oh  :: a
                   }
  deriving (Show, Eq, Functor)

-- | Atom environment, e.g. hydrogens or radicals
--
data Env r a = Env { _atom'       :: a
                   , _environment :: r a
                   }
  deriving (Show, Eq, Functor)

-- | Hydrogens envrironment
--
type H a = Env [] a

-- | Oxigen and hydroxi group, connected to C-terminal of amino acid
--
data OXT a = OXT { _o'   :: a
                 , _oxt' :: a
                 }
  deriving (Show, Eq, Functor)

-- | CG atom with radical type
--
data CG a = CG { _cg'      :: a
               , _radical' :: AA
               }
  deriving (Show, Eq, Functor)

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