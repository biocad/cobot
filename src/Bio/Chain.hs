{-# LANGUAGE TupleSections #-}
module Bio.Chain
    ( ChainLike (..)
    , Chain
    , chain, fromList
    , (!), (//)
    ) where

import           Control.Lens
import qualified Data.Array                as A ( bounds
                                                , assocs
                                                )
import           Data.Array                     ( Array
                                                , Ix
                                                , array
                                                , listArray
                                                , (!)
                                                , (//)
                                                )

type Chain i a = Array i a

-- | Construct new chain from list
--
chain :: Ix i => (i, i) -> [(i, a)] -> Chain i a
chain = array

-- | Construct new int-labeled chain from list
--
fromList :: [a] -> Chain Int a
fromList lst = listArray (0, length lst - 1) lst

-- | Chain-like sequence, by default it is an array or a list
--
class (Ixed m, Enum (Index m)) => ChainLike m where
    bounds       :: m -> (Index m, Index m)
    assocs       :: m -> [(Index m, IxValue m)]
    modify       :: Index m -> (IxValue m -> IxValue m) -> m -> m
    modifyBefore :: Index m -> (IxValue m -> IxValue m) -> m -> m
    modifyAfter  :: Index m -> (IxValue m -> IxValue m) -> m -> m

instance ChainLike [a] where
    bounds = (0,) . pred . length
    
    assocs  = zip [0..]

    modify       _ _ []      = []
    modify       0 f (x:xs)  = f x:xs
    modify       i f (x:xs)  = x:modify (i - 1) f xs

    modifyBefore i f lst = (f <$> take i lst) ++ drop i lst
    modifyAfter  i f lst = take (i + 1) lst ++ (f <$> drop (i + 1) lst)

instance (Ix i, Enum i) => ChainLike (Array i a) where
    bounds = A.bounds

    assocs = A.assocs

    modify       i f ar = ar // [(i, f (ar ! i))]

    modifyBefore i f ar = let (mi, _) = bounds ar
                          in ar // [(j, f (ar ! j)) | j <- [mi .. pred i]]
    modifyAfter  i f ar = let (_, ma) = bounds ar
                          in ar // [(j, f (ar ! j)) | j <- [succ i .. ma]]