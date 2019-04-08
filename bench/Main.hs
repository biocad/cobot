{-# LANGUAGE StandaloneDeriving #-}
{-# OPTIONS_GHC -Wno-orphans #-}

module Main where

import Bio.Chain.Alignment          (SimpleGap, AffineGap (..), LocalAlignment (..),
                                     GlobalAlignment (..), SemiglobalAlignment (..),
                                     AlignmentResult (..), Operation (..), align)
import Bio.Chain                    (Chain, fromList)
import Bio.Chain.Alignment.Scoring  (nuc44)
import Bio.Protein.AminoAcid        (AA (..))
import Control.DeepSeq              (NFData (..), force, deepseq)
import Control.Lens                 (Index)
import Control.Monad.State          (State, evalState, state)
import GHC.Generics                 (Generic)
import System.Random                (RandomGen, randomR, getStdGen)
import Control.Monad                (replicateM)
import System.Clock                 (getTime, Clock ( Monotonic ), diffTimeSpec, sec, nsec)
import Control.Parallel.Strategies  (rdeepseq, rpar, dot, parListChunk, withStrategy)
import GHC.Conc                     (numCapabilities)
import Criterion
import Criterion.Main

deriving instance Generic AA

instance NFData AA

instance (NFData a, NFData b) => NFData (Operation a b) where
    rnf (DELETE i) = rnf i
    rnf (INSERT j) = rnf j
    rnf (MATCH i j) = rnf (i, j)

instance (NFData a, NFData b, NFData (Index a), NFData (Index b))
         => NFData (AlignmentResult a b) where
    rnf (AlignmentResult score alignment s1 s2) = rnf (score', alignment', s1', s2')
      where
        score' = force score
        alignment' = force alignment
        s1' = force s1
        s2' = force s2

makeRandomChain :: RandomGen g => Int -> State g String
makeRandomChain 0 = pure ""
makeRandomChain len = do
    c <- ("ATGC" !!) <$> state (randomR (0, 3))
    cs <- makeRandomChain (len - 1)
    pure (c : cs)

makeRandomChainIO :: Int -> IO (Chain Int Char)
makeRandomChainIO len = do
    list <- evalState (makeRandomChain len) <$> getStdGen
    pure (fromList list)

measureTime :: NFData a => String -> a -> IO ()
measureTime label value = do
    t1 <- getTime Monotonic
    t2 <- value `deepseq` getTime Monotonic
    let dt = diffTimeSpec t2 t1
    let timeInSeconds = fromIntegral (sec dt) + fromIntegral (nsec dt) / 1000000000 :: Double
    let padding = " " <> replicate (max 0 (50 - length label)) '.' <> " "
    putStrLn $ label <> padding <> show timeInSeconds <> "s"

parMap' :: NFData b => Int -> (a -> b) -> [a] -> [b]
parMap' chunkSize f = withStrategy (parListChunk chunkSize (rdeepseq `dot` rpar)) . map f

setupEnv :: IO (Chain Int Char, [Chain Int Char], Int)
setupEnv = do
    a <- makeRandomChainIO 100
    bs <- replicateM 100 $ makeRandomChainIO 100
    let chunkSize = length bs `div` numCapabilities
    pure (a, bs, chunkSize)

main :: IO ()
main = defaultMain [
        env setupEnv $ \ ~(a, bs, chunkSize) -> bgroup "main" [
            bench "Local alignment" $
                let align' = align (LocalAlignment nuc44 (-10 :: SimpleGap)) a
                in  nfIO . pure $ parMap' chunkSize align' bs,
            bench "Global alignment" $
                let align' = align (GlobalAlignment nuc44 (-10 :: SimpleGap)) a
                in  nfIO . pure $ parMap' chunkSize align' bs,
            bench "Semiglobal alignment" $
                let align' = align (SemiglobalAlignment nuc44 (-10 :: SimpleGap)) a
                in  nfIO . pure $ parMap' chunkSize align' bs,
            bench "Local alignment with affine gap" $
                let align' = align (LocalAlignment nuc44 (AffineGap (-10) (-1))) a
                in  nfIO . pure $ parMap' chunkSize align' bs,
            bench "Global alignment with affine gap" $
                let align' = align (GlobalAlignment nuc44 (AffineGap (-10) (-1))) a
                in  nfIO . pure $ parMap' chunkSize align' bs,
            bench "Semiglobal alignment with affine gap" $
                let align' = align (SemiglobalAlignment nuc44 (AffineGap (-10) (-1))) a
                in  nfIO . pure $ parMap' chunkSize align' bs
        ]
    ]
