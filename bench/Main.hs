module Main where

import           Bio.Chain                   (Chain, fromList)
import           Bio.Chain.Alignment         (AffineGap (..),
                                              GlobalAlignment (..),
                                              LocalAlignment (..),
                                              SemiglobalAlignment (..),
                                              SimpleGap, align)
import           Bio.Chain.Alignment.Scoring (nuc44)
import           Control.DeepSeq             (NFData (..), deepseq)
import           Control.Monad               (replicateM)
import           Control.Monad.State         (State, evalState, state)
import           Control.Parallel.Strategies (dot, parListChunk, rdeepseq, rpar,
                                              withStrategy)
import           Criterion                   (bench, bgroup, env, nfIO)
import           Criterion.Main              (defaultMain)
import           GHC.Conc                    (numCapabilities)
import           System.Clock                (Clock (Monotonic), diffTimeSpec,
                                              getTime, nsec, sec)
import           System.Random               (RandomGen, getStdGen, randomR)

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
    a <- makeRandomChainIO 600
    bs <- replicateM 20 $ makeRandomChainIO 4500
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
