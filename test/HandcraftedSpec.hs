module HandcraftedSpec where

import           Bio.Chain.Alignment
import           Bio.Chain.Alignment.Scoring
import           Test.Hspec

handcraftedTests :: Spec
handcraftedTests = do
  handcraftedLocal

handcraftedLocal :: Spec
handcraftedLocal = do
  describe "LocalAlignment. AffineGap2" $ do
    let testAlign = align (LocalAlignment nuc44 (AffineGap (-5) (-1), AffineGap (-1000) (-1000)))
    let a = "AAAAAAGGGGGGGGGGGGTTTTTTTTT"
        b = "AAAAAATTTTTTTTT"
        aAns = "TTTTTTTTT"
        bAns = "TTTTTTTTT"
        res = testAlign (a :: String) (b :: String)
        (a', b') = viewAlignment res
        score1 = score res
    it "first sequence"  $ a'     `shouldBe` aAns
    it "second sequence" $ b'     `shouldBe` bAns
    it "score"           $ score1 `shouldBe` 45

    let a2 = "AAAAAATTTTTTTTT"
        b2 = "AAAAAAGGGGGGGGGGGGTTTTTTTTT"
        a2Ans = "AAAAAA------------TTTTTTTTT"
        b2Ans = "AAAAAAGGGGGGGGGGGGTTTTTTTTT"
        res2 = testAlign (a2 :: String) (b2 :: String)
        (a2', b2') = viewAlignment res2
        score2 = score res2
    it "first sequence"  $ a2'    `shouldBe` a2Ans
    it "second sequence" $ b2'    `shouldBe` b2Ans
    it "score"           $ score2 `shouldBe` 59

    let a3 = "AACCTTCCAACCGGCCAACCTTCCAACCGGCCTT"
        b3 = "AATTAAGGAATTAAGGTT"
        a3Ans = "GGCCTT"
        b3Ans = "GGAATT"
        res3 = testAlign (a3 :: String) (b3 :: String)
        (a3', b3') = viewAlignment res3
        score3 = score res3
    it "first sequence"  $ a3'    `shouldBe` a3Ans
    it "second sequence" $ b3'    `shouldBe` b3Ans
    it "score"           $ score3 `shouldBe` 12

    let a4 = "AATTAAGGAATTAAGGTT"
        b4 = "AACCTTCCAACCGGCCAACCTTCCAACCGGCCTT"
        a4Ans = "AA--TT--AA--GG--AA--TT--AA--GG--TT"
        b4Ans = "AACCTTCCAACCGGCCAACCTTCCAACCGGCCTT"
        res4 = testAlign (a4 :: String) (b4 :: String)
        (a4', b4') = viewAlignment res4
        score4 = score res4
    it "first sequence"  $ a4'    `shouldBe` a4Ans
    it "second sequence" $ b4'    `shouldBe` b4Ans
    it "score"           $ score4 `shouldBe` 42
