-- |
-- This module contains basic classes with read and write functions.
-- For this moment we suppose that basic IO representation is 'Text' and other IO respresentations can be obtained from it.
-- Thus for every format you should know how to read it from 'Text' and write to 'Text'.

module Bio.Utils.IO.Class
  (
  -- * Readers
    FromText (..)
  , FromFile (..)
  -- * Writers
  , ToText (..)
  , ToFile (..)
  ) where

import           Data.Text    (Text)
import           Data.Text.IO (appendFile, readFile, writeFile)
import           Prelude      hiding (appendFile, readFile, writeFile)

-- | Ability to read from `Text`.
-- 
class FromText a where
  fromText :: Text -> a

-- | Ability to read from file.
--
class FromFile a where
  fromFile :: FilePath -> IO a

instance FromText a => FromFile a where
  fromFile f = fromText <$> readFile f

-- | Ability to write to `Text`.
--
class ToText a where
  toText :: a -> Text

-- | Ability to write to file.
--
class ToFile a where
  toFile       :: a -> FilePath -> IO ()
  appendToFile :: a -> FilePath -> IO ()

instance ToText a => ToFile a where
  x `toFile`       f = writeFile  f $ toText x
  x `appendToFile` f = appendFile f $ toText x
