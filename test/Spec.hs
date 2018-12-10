{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TypeApplications  #-}

import           Test.Hspec

import Data.Text
import Data.Array
import Control.Lens
import Linear.Metric
import Linear.Epsilon
import Bio.Internal.Structure
import Bio.Protein.AminoAcid
import Bio.Protein.Builder

cosBV :: V3R -> V3R -> R
cosBV a b = dot a b / (norm a * norm b)

buildChainSpec :: Spec
buildChainSpec =
  describe "Chain builder (BBT)" $ do
    let chain = build [ALA, CYS, ASP] :: Array Int (BBT V3R)
    let aa1 = chain ! 0
        aa2 = chain ! 1
        aa3 = chain ! 2
    let nac = cos $ angle N CA C
        acn = cos $ angle CA C N
        cna = cos $ angle C N CA
        na  = dist N CA
        ac  = dist CA C
        cn  = dist C N
    let a_  = aa1 ^. ca
        n_  = aa1 ^. n
        c_  = aa1 ^. c
    let a2_ = aa2 ^. ca
        n2_ = aa2 ^. n
        c2_ = aa2 ^. c
    let a3_ = aa3 ^. ca
        n3_ = aa3 ^. n
        c3_ = aa3 ^. c
    it "builds single amino acid" $ do
        distance n_ a_ - na `shouldSatisfy` nearZero
        distance a_ c_ - ac `shouldSatisfy` nearZero
        cosBV (a_ - n_) (a_ - c_) - nac `shouldSatisfy` nearZero
    it "builds even amino acids" $ do
        distance n2_ c_  - cn `shouldSatisfy` nearZero
        cosBV (c_ - a_) (c_ - n2_) - acn `shouldSatisfy` nearZero
        cosBV (n2_ - c_) (n2_ - a2_) - cna `shouldSatisfy` nearZero
        distance n2_ a2_ - na `shouldSatisfy` nearZero
        distance a2_ c2_ - ac `shouldSatisfy` nearZero
        cosBV (a2_ - n2_) (a2_ - c2_) - nac `shouldSatisfy` nearZero
    it "builds odd amino acids" $ do
        distance n3_ c2_  - cn `shouldSatisfy` nearZero
        cosBV (c2_ - a2_) (c2_ - n3_) - acn `shouldSatisfy` nearZero
        cosBV (n3_ - c2_) (n3_ - a3_) - cna `shouldSatisfy` nearZero
        distance n3_ c2_ - cn `shouldSatisfy` nearZero
        distance n3_ a3_ - na `shouldSatisfy` nearZero
        distance a3_ c3_ - ac `shouldSatisfy` nearZero
        cosBV (a3_ - n3_) (a3_ - c3_) - nac `shouldSatisfy` nearZero

lensesSpec :: Spec
lensesSpec =
  describe "Amino acid lenses" $ do
    it "works on BB" $ do
        let aa = create @(BB Text) "N" "CA" "C"
        aa ^. n  `shouldBe` "N"
        aa ^. ca `shouldBe` "CA"
        aa ^. c  `shouldBe` "C"
    it "works on BBT" $ do
        let aa = create @(BBT Text) "N" "CA" "C" ALA
        aa ^. n       `shouldBe` "N"
        aa ^. ca      `shouldBe` "CA"
        aa ^. c       `shouldBe` "C"
        aa ^. radical `shouldBe` ALA
    it "works on BBO" $ do
        let aa = create @(BBO Text) "N" "CA" "C" "O"
        aa ^. n  `shouldBe` "N"
        aa ^. ca `shouldBe` "CA"
        aa ^. c  `shouldBe` "C"
        aa ^. o  `shouldBe` "O"

main :: IO ()
main = hspec $ do
    lensesSpec
    buildChainSpec