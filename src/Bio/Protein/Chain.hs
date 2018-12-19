module Bio.Protein.Chain
    ( ChainLike (..)
    , Chain
    , fromList
    , (!), (//)
    ) where

import           Control.Lens
import           Data.Array                     ( Array
                                                , Ix
                                                , bounds
                                                , listArray
                                                , (!)
                                                , (//)
                                                )

type Chain a = Array Int a

fromList :: [a] -> Chain a
fromList lst = listArray (0, length lst - 1) lst

-- | Chain-like sequence, by default it is an array or a list
--
class (Ixed m, Enum (Index m)) => ChainLike m where
    modify       :: Index m -> (IxValue m -> IxValue m) -> m -> m
    modifyBefore :: Index m -> (IxValue m -> IxValue m) -> m -> m
    modifyAfter  :: Index m -> (IxValue m -> IxValue m) -> m -> m

instance ChainLike [a] where
    modify       _ _ []      = []
    modify       0 f (x:xs)  = f x:xs
    modify       i f (x:xs)  = x:modify (i - 1) f xs

    modifyBefore i f lst = (f <$> take i lst) ++ drop i lst
    modifyAfter  i f lst = take (i + 1) lst ++ (f <$> drop (i + 1) lst)

instance (Ix i, Enum i) => ChainLike (Array i a) where
    modify       i f ar = ar // [(i, f (ar ! i))]

    modifyBefore i f ar = let (mi, _) = bounds ar
                          in ar // [(j, f (ar ! j)) | j <- [mi .. pred i]]
    modifyAfter  i f ar = let (_, ma) = bounds ar
                          in ar // [(j, f (ar ! j)) | j <- [succ i .. ma]]