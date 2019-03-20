{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TypeApplications  #-}

import           Test.Hspec

import           Data.Text
import           Control.Lens

import           Bio.Utils.Geometry
import           Bio.Protein.AminoAcid
import           Bio.Protein.Chain
import           Bio.Protein.Chain.Builder

buildChainSpec :: Spec
buildChainSpec = describe "Chain builder (BBT)" $ do
    let chain = build [ALA, CYS, ASP] :: ProteinChain Int (BBT V3R)
    let aa1   = chain ^?! ix 0
        aa2   = chain ^?! ix 1
        aa3   = chain ^?! ix 2
    let nac   = pi * 110.990 / 180.0
        acn   = pi * 118.995 / 180.0
        cna   = pi * 118.995 / 180.0
        na    = 1.460
        ac    = 1.509
        cn    = 1.290
    let a_    = aa1 ^. ca . atom
        n_    = aa1 ^. n . atom
        c_    = aa1 ^. c . atom
    let a2_   = aa2 ^. ca . atom
        n2_   = aa2 ^. n . atom
        c2_   = aa2 ^. c . atom
    let a3_   = aa3 ^. ca . atom
        n3_   = aa3 ^. n . atom
        c3_   = aa3 ^. c . atom
    it "builds single amino acid" $ do
        distance n_ a_ - na `shouldSatisfy` nearZero
        distance a_ c_ - ac `shouldSatisfy` nearZero
        angle (a_ - n_) (a_ - c_) - nac `shouldSatisfy` nearZero
    it "builds even amino acids" $ do
        distance n2_ c_ - cn `shouldSatisfy` nearZero
        angle (c_ - a_) (c_ - n2_) - acn `shouldSatisfy` nearZero
        angle (n2_ - c_) (n2_ - a2_) - cna `shouldSatisfy` nearZero
        distance n2_ a2_ - na `shouldSatisfy` nearZero
        distance a2_ c2_ - ac `shouldSatisfy` nearZero
        angle (a2_ - n2_) (a2_ - c2_) - nac `shouldSatisfy` nearZero
    it "builds odd amino acids" $ do
        distance n3_ c2_ - cn `shouldSatisfy` nearZero
        angle (c2_ - a2_) (c2_ - n3_) - acn `shouldSatisfy` nearZero
        angle (n3_ - c2_) (n3_ - a3_) - cna `shouldSatisfy` nearZero
        distance n3_ c2_ - cn `shouldSatisfy` nearZero
        distance n3_ a3_ - na `shouldSatisfy` nearZero
        distance a3_ c3_ - ac `shouldSatisfy` nearZero
        angle (a3_ - n3_) (a3_ - c3_) - nac `shouldSatisfy` nearZero

lensesSpec :: Spec
lensesSpec = describe "Amino acid lenses" $ do
    it "works on BB" $ do
        let aa = create @(BB Text) "N" "CA" "C"
        aa ^. n . atom `shouldBe` "N"
        aa ^. ca . atom `shouldBe` "CA"
        aa ^. c . atom `shouldBe` "C"
    it "works on BBT" $ do
        let aa = create @(BBT Text) "N" "CA" "C" ALA
        aa ^. n . atom `shouldBe` "N"
        aa ^. ca . atom `shouldBe` "CA"
        aa ^. c . atom `shouldBe` "C"
        aa ^. radical `shouldBe` ALA
    it "works on BBO" $ do
        let aa = create @(BBO Text) "N" "CA" "C" "O"
        aa ^. n . atom `shouldBe` "N"
        aa ^. ca . atom `shouldBe` "CA"
        aa ^. c . atom `shouldBe` "C"
        aa ^. o . atom `shouldBe` "O"

main :: IO ()
main = hspec $ do
    lensesSpec
    buildChainSpec
