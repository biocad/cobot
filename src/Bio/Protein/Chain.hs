module Bio.Protein.Chain
    ( Chain
    , fromList
    , (!)
    ) where

import           Data.Array                     ( Array
                                                , listArray
                                                , (!)
                                                )

type Chain a = Array Int a

fromList :: [a] -> Chain a
fromList lst = listArray (0, length lst - 1) lst
