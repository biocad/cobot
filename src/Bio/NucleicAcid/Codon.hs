module Bio.NucleicAcid.Codon where

import Bio.NucleicAcid.Nucleotide (DNA (..), RNA (..))

type Codon a = (a, a, a)

type DNACodon = Codon DNA
type RNACodon = Codon RNA

class Eq a => CodonLike a where
    startCodon :: Codon a
    amberCodon :: Codon a
    ochreCodon :: Codon a
    opalCodon  :: Codon a

    isStart :: Codon a -> Bool
    isStart = (==) startCodon

    isStop :: Codon a -> Bool
    isStop = flip elem [amberCodon, ochreCodon, opalCodon]

instance CodonLike DNA where
    startCodon = (DA, DT, DG)
    amberCodon = (DT, DA, DG)
    ochreCodon = (DT, DA, DA)
    opalCodon  = (DT, DG, DA)

instance CodonLike RNA where
    startCodon = (RA, RU, RG)
    amberCodon = (RU, RA, RG)
    ochreCodon = (RU, RA, RA)
    opalCodon  = (RU, RG, RA)
