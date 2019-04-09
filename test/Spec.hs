{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TypeApplications  #-}

import           Test.Hspec

import           Data.Text
import           Control.Lens

import           Bio.Chain.Alignment
import           Bio.Utils.Geometry
import           Bio.Protein.AminoAcid
import           Bio.Protein.Chain
import           Bio.Protein.Chain.Builder
import           Bio.Chain.Alignment.Type

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

alignLocal :: String -> String -> AlignmentResult String String
alignLocal = align (LocalAlignment nuc44 (-10 :: SimpleGap))

alignGlobal :: String -> String -> AlignmentResult String String
alignGlobal = align (GlobalAlignment nuc44 (-10 :: SimpleGap))

alignSemiglobal :: String -> String -> AlignmentResult String String
alignSemiglobal = align (SemiglobalAlignment nuc44 (-10 :: SimpleGap))

alignLocalAffine:: String -> String -> AlignmentResult String String
alignLocalAffine = align (LocalAlignment nuc44 (AffineGap (-10) (-1)))

alignGlobalAffine :: String -> String -> AlignmentResult String String
alignGlobalAffine = align (GlobalAlignment nuc44 (AffineGap (-10) (-1)))

alignSemiglobalAffine :: String -> String -> AlignmentResult String String
alignSemiglobalAffine = align (SemiglobalAlignment nuc44 (AffineGap (-10) (-1)))

alignAllSpec :: String -> String -> [(String, String)] -> Spec
alignAllSpec a b expected = do
    alignLocalSpec a b expected
    alignGlobalSpec a b expected
    alignSemiglobalSpec a b expected

alignGlobalSpec :: String -> String -> [(String, String)] -> Spec
alignGlobalSpec a b expected = do
    it "Global" $ do
        expected `shouldContain` [viewAlignment (alignGlobal a b)]
    it "Global with Affine gap " $ do
        expected `shouldContain` [viewAlignment (alignGlobalAffine a b)]

alignLocalSpec :: String -> String -> [(String, String)] -> Spec
alignLocalSpec a b expected = do
    it "Local" $ do
        expected `shouldContain` [viewAlignment (alignLocal a b)]
    it "Local with Affine gap" $ do
        expected `shouldContain` [viewAlignment (alignLocalAffine a b)]

alignSemiglobalSpec :: String -> String -> [(String, String)] -> Spec
alignSemiglobalSpec a b expected = do
    it "Semiglobal" $ do
        expected `shouldContain` [viewAlignment (alignSemiglobal a b)]
    it "Semiglobal with Affine gap" $ do
        expected `shouldContain` [viewAlignment (alignSemiglobalAffine a b)]

smallAlignmentSpec :: Spec
smallAlignmentSpec = describe "Should work for all algorithms" $ do
    describe "Align single letter over itself" $ do
        alignAllSpec "A" "A" [("A", "A")]
    describe "Align several letters over itself" $ do
        alignAllSpec "ATGC" "ATGC" [("ATGC", "ATGC")]
    describe "Align single letter over mismatched one" $ do
        let a = "A"
        let b = "G"
        alignLocalSpec a b [("", "")]
        alignSemiglobalSpec a b [("-A", "G-"), ("A-", "-G")]
        alignGlobalSpec a b [("A", "G"), ("-A", "G-"), ("A-", "-G")]
    describe "Align multiple letters without match" $ do
        let a = "ATATA"
        let b = "GCGCG"
        let expected = [("-----ATATA", "GCGCG-----"), ("ATATA-----", "-----GCGCG")]
        alignLocalSpec a b [("", "")]
        alignSemiglobalSpec a b expected
        alignGlobalSpec a b (("ATATA", "GCGCG") : expected)
    describe "Align internal similarity" $ do
        let a = "ATGCATGCATGC"
        let b = "CATGCA"
        let expected = [("ATGCATGCATGC", "---CATGCA---")]
        alignLocalSpec a b [("CATGCA", "CATGCA")]
        alignSemiglobalSpec a b expected
        alignGlobalSpec a b expected
    describe "Align prefix similarity" $ do
        let a = "ATGCATGCATGC"
        let b = "ATGCATGCA"
        let expected = [("ATGCATGCATGC", "ATGCATGCA---")]
        alignLocalSpec a b [("ATGCATGCA", "ATGCATGCA")]
        alignSemiglobalSpec a b expected
        alignGlobalSpec a b expected
    describe "Align suffix similarity" $ do
        let a = "ATGCATGCATGC"
        let b = "CATGCATGC"
        let expected = [("ATGCATGCATGC", "---CATGCATGC")]
        alignLocalSpec a b [("CATGCATGC", "CATGCATGC")]
        alignSemiglobalSpec a b expected
        alignGlobalSpec a b expected

semiglobalSpec :: Spec
semiglobalSpec = describe "Semiglobal alignment" $ do
    it "may not end with MATCH (single letter)" $ do
        let result = alignSemiglobal "A" "T"
        viewAlignment result `shouldBe` ("A-", "-T")
    it "may not end with MATCH (many letters)" $ do
        let result = alignSemiglobal "ATAT" "GCGC"
        viewAlignment result `shouldBe` ("ATAT----", "----GCGC")
    it "may not end with MATCH (many letters and have match in the middle)" $ do
        let result = alignSemiglobal "AAATT" "TTGGG"
        viewAlignment result `shouldBe` ("AAATT---", "---TTGGG")
    it "may not end with MATCH (many letters and have match in the middle and gaps)" $ do
        let result = alignSemiglobal "CCCCCCCATGAGATAGATA" "ATGACCGATAGATAGGGGGG"
        let a = "CCCCCCCATGA--GATAGATA------"
        let b = "-------ATGACCGATAGATAGGGGGG"
        viewAlignment result `shouldBe` (a, b)

affineSpec :: Spec
affineSpec = describe "Affine gap and simple gap work differently" $ do
    it "local alignment" $ do
        let a = "CCCCCCATGATGACCCCCC"
        let b = "TTTTTTATGAGGGGTGAGGGGGG"
        viewAlignment (alignLocal a b) `shouldBe` ("TGATGA", "TTATGA")
        viewAlignment (alignLocalAffine a b) `shouldBe` ("ATGA----TGA", "ATGAGGGGTGA")
    it "global alignment" $ do
        let a = "AAAAAAAAATTTTTTTTTTATTTTTTTTTAATTTTTTTTTTTTTTTTTTT"
        let b = "TTTTTTTTTTATAAAAAAAAAAATTTTTTTTTTTTATTTTTTTTTTTTT"
        let a' = "AAAAAAAAATTTTTTTTTTATTTTTTTTTAATTTTTTTTTTTTTTTTTTT"
        let b' = "-TTTTTTTTTTATAAAAAAAAAAATTTTTTTTTTTTATTTTTTTTTTTTT"
        let a'' = "AAAAAAAAATTTTTTTTTTAT-----------TTTTTTTTAATTTTTTTTTTTTTTTTTTT"
        let b'' = "---------TTTTTTTTTTATAAAAAAAAAAATTTTTTTT---TTTTATTTTTTTTTTTTT"
        viewAlignment (alignGlobal a b) `shouldBe` (a', b')
        viewAlignment (alignGlobalAffine a b) `shouldBe` (a'', b'')
    it "semiglobal alignment" $ do
        let a = "ATAAAAAAAAAAATTTTTTTTTTTATTTTTTTTTTTTATTTTTTTTTTTTTT"
        let b = "AAATTTTTTTTTTTTTTTTATTTTTTTTTTTTAATTTTTTTTTTTT"
        let a' = "ATAAAAAAAAAAATTTTTTTTTTTATTTTTTTTTTTTATTTTTTTTTTTTTT"
        let b' = "-----AAATTTTTTTTTTTTTTTTATTTTTTTTTTTTAATTTTTTTTTTTT-"
        let a'' = "ATAAAAAAAAAAATTTTTTTTTTTATTTTTTTTTTTTATTTTTTTTTTTTTT"
        let b'' = "-----AAATTTTTTTTTTTTTTTTATTTTTTTTTTTTAATTTTTTTTTTTT-"
        viewAlignment (alignSemiglobal a b) `shouldBe` (a', b')
        viewAlignment (alignSemiglobalAffine a b) `shouldBe` (a'', b'')

alignmentSanitySpec  :: Spec
alignmentSanitySpec = describe "Sanity tests" $ do
    describe "Local alignment works" $ do
        let a = "TTTTTTTTAAAAAATTTTTTTTTTTT"
        let b = "CCCCCCCCCCCCCAAAAAACCCCCCCCCCC"
        let a' = "AAAAAA"
        let b' = "AAAAAA"
        alignLocalSpec a b [(a', b')]
    describe "Global alignment works" $ do
        let a = "AAAAAATAAAAAACAAAAAGAAA"
        let b = "AAAAGAAAAGAAAAAAACAAAAA"
        let a' = "AAAAAATAAAAAACAAAAAGAAA"
        let b' = "AAAAGAAAAGAAAAAAACAAAAA"
        alignGlobalSpec a b [(a', b')]
    describe "Semiglobal alignment works" $ do
        let a = "AAAAAAAAAATTTTTTTTTTT"
        let b = "TTTTTTTTTCCCCCCCCCCC"
        let a' = "AAAAAAAAAATTTTTTTTTTT-----------"
        let b' = "------------TTTTTTTTTCCCCCCCCCCC"
        alignSemiglobalSpec a b [(a', b')]

alignmentSpec :: Spec
alignmentSpec = describe "Alignment" $ do
    affineSpec
    semiglobalSpec
    smallAlignmentSpec
    alignmentSanitySpec

main :: IO ()
main = hspec $ do
    lensesSpec
    buildChainSpec
    alignmentSpec
