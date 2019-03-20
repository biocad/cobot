module Bio.Molecule
  ( Molecule(..)
  , MoleculeLike(..)
  , singleton
  ) where

import           Control.Lens                   ( (^?)
                                                , Index
                                                , IxValue
                                                , Ixed (..)
                                                , lens
                                                , (&)
                                                , (.~)
                                                )

newtype Molecule t c = Molecule { getChains :: [(t, c)] }
  deriving (Show, Eq)

type instance Index (Molecule t c) = t
type instance IxValue (Molecule t c) = c

class (Eq (Index m), Ixed m) => MoleculeLike m where
    -- | Create empty molecule without chains
    --
    empty :: m
    -- | Delete chain with specified index (returns error if chain doesn't present)
    --
    deleteAt :: m -> Index m -> m
    -- | Create chain with specified index (returns error if chain is already present)
    --
    create :: m -> Index m -> IxValue m -> m
    -- | Set new chain with speficied index (creates new if does not present)
    --
    set :: m -> Index m -> IxValue m -> m

-- | Create molecule with single chain
--
singleton :: MoleculeLike m => Index m -> IxValue m -> m
singleton = create empty

instance Eq t => Ixed (Molecule t c) where
    ix idx = lens (lookup idx . getChains) (\(Molecule m) my -> Molecule $ setL my m) . traverse
      where
        setL :: Maybe c -> [(t, c)] -> [(t, c)]
        setL Nothing  xs = xs
        setL (Just _) [] = error "Chain should be present"
        setL y@(Just a) ((x', y') : xs) | x' == idx = (idx, a) : xs
                                        | otherwise = (x', y') : setL y xs

instance Eq t => MoleculeLike (Molecule t c) where
    empty = Molecule []

    deleteAt (Molecule xs) idx = Molecule $ deleteFromList xs
      where
        deleteFromList :: [(t, c)] -> [(t, c)]
        deleteFromList [] = error "Chain is not present"
        deleteFromList (a@(x', _) : ys) | x' == idx = ys
                                        | otherwise = a : deleteFromList ys

    create (Molecule xs) idx c = Molecule $ createInList xs
      where
        createInList :: [(t, c)] -> [(t, c)]
        createInList [] = [(idx, c)]
        createInList (a@(x', _) : ys)
            | x' == idx = error "Chain should not be present at molecule"
            | otherwise = a : createInList ys

    set m idx c = case m ^? ix idx of
        Nothing -> create m idx c
        Just _  -> m & ix idx .~ c
