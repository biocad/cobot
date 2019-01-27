{- |

This module contains realization of PDB format, which is one of the main formats for the [Protein Data Bank](https://www.rcsb.org/).
Read the docs about this format [here](http://www.wwpdb.org/documentation/file-format).

-}

module Bio.Utils.IO.PDB
  ( -- * About PDB
    -- $aboutPDB
    -- * Type
    module Bio.Utils.IO.PDB.Type
  , module Bio.Utils.IO.PDB.Writer
  , module Bio.Utils.IO.PDB.Reader
  , module Bio.Utils.IO.Class
  ) where

import           Bio.Utils.IO.PDB.Type
import           Bio.Utils.IO.PDB.Writer
import           Bio.Utils.IO.PDB.Reader
import           Bio.Utils.IO.Class

{- $aboutPDB

In two words, PDB format consist from ATOM lines (that contains information about atom) and other lines, that are used for the meta-information.
Our main interest will be in area of atom lines. The following description about them.

PDB format has the following tree:

> - PDB file
>   - Chains
>     - Residues
>       - Atoms

In other words, PDB file consists from chains, every chain consists from residues, every residue consists from atoms.
Every atom is one ATOM line in file.

Actually, this is one more layer between PDB file and chains - MODEL. 
But models are out from ATOM lines, they are not nesessary by specification, and there are a lot of cases when they are not written by authors of PDB file. 
So we would not take them.

-}

{- $writer

Module writer consist from two parts.

(1) First part is how to convert list of 'PDBLine's into 'Text'.

(1) Second part is how to convert your own types into list of 'PDBLine's (then you can write them into any IO automatically).
Actually, due to PDB format, we can reduce code a little bit.
Main idea is that you have to write only converter for residue (how convert your type that corresponds to PDB residue into '[PDBLine]').

-}

{- $reader

-}
