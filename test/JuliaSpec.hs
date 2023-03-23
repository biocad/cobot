module JuliaSpec where

import           Bio.Chain.Alignment
import           Bio.Chain.Alignment.Scoring
import           Test.Hspec

juliaBasedTests :: Spec
juliaBasedTests = do
  juliaGlobal
  juliaSemiglobal
  juliaLocal

juliaGlobal :: Spec
juliaGlobal = do
  describe "Global Alignment, Simple Gap" $ do
    let testAlign = align (GlobalAlignment nuc44 (-9 :: Int))

    let resBig1   = testAlign (bigA :: String) (bigB :: String)
        resBig2   = testAlign (bigB :: String) (bigA :: String)
        resMed1   = testAlign (bigA :: String) (bigB :: String)
        resMed2   = testAlign (bigB :: String) (bigA :: String)
        resSmall1 = testAlign (bigA :: String) (bigB :: String)
        resSmall2 = testAlign (bigB :: String) (bigA :: String)

    it "Is symmetric (big): " $   score resBig1   `shouldBe` score resBig2
    it "Is symmetric (med): " $   score resMed1   `shouldBe` score resMed2
    it "Is symmetric (small): " $ score resSmall1 `shouldBe` score resSmall2
    -- scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=0, gap_extend=-9);
    -- pairalign(GlobalAlignment(), a, b, scoremodel)
    let aAns = "CTCTGCAAGCATAGAGATATTTCACCGGCAATATTTTTCGTTGAAGTGTATTTGTCCCTATTATCACACC\
                \AGTCATTCCGATGTCGTTAGGGCCGCCTGTTTTACTGGAGTTCCCGTTCGATGTACCCTCCTCCATAGA\
                \GATTTAGGGGTATAGACCGCGACGGACGGGGCTTGCGAAAGCTGCCCGATAATGGACTTTGACAATAGA\
                \CGTAATCCGTGAACGGCTGGCGTTTCACACTCAACCGCATTGAACAGCCATGCCATCCGCAGGATGCGC\
                \CAGGGAAAAGCCACCCCCAAGAACGACTCCAGGGGCCCAAGTATGATGATAGCGTTAGGCTTGCTCGAA\
                \GGCAAACAACTCCGAGGATGCTTAGCAGGCCACAACGTGCTCGATAGCACGAGCGTAAGATGCGAACGC\
                \AGTAGTACCTAATGCGGGCGAGACATTTCGGGCCGTTGAAGCCCTGTCTGTTTTTTCGTACAGAGGGTGGTTTGAAGTATCGCCC"
        bAns = "-T--G-AAGCA--GC-AT-TCG-ACT--CACTCATGT-CGCCTA-GAG-A---GAA--T-TTAT-A----\
                \-GTCAA----ATGAC-TT-GTG--GC--GTTCT--TGGACTT----TTAGA-G--CCGAAATGCAT-G-\
                \GA----GG---A-A-AC----A----------TT---AA-G--G----ATG-T--A-TTTG----T---\
                \--TA--C--TGA--G-C---------A---T-----G-ATT-AA--GC-A----AT--------TGC-C\
                \-AGGTAAA--C----------------T---G---CC-A-G---G-T--T----T--GG--------A-\
                \--C---CAA-TC------T--T--G--GG--A----G--CT---TAGTAC---C-TAA-A-GA-AA-G-\
                \A--AGT-C--A-TGCA--C-AGTCA---CG--CCGC--A--CC----------------A-A-----T------A--TA-CAC--"
        res = testAlign (bigA :: String) (bigB :: String)
        (a', b') = viewAlignment res
        score1 = score res
    it "first sequence"  $ a'     `shouldBe` aAns
    it "second sequence" $ b'     `shouldBe` bAns
    it "score"           $ score1 `shouldBe` (-1605)

    let a2Ans = "-GT-TTGGG-CATATTC-AAG-ATCA-GAC-C-A----ATCGGTC-GAT-G--TGAAACAAGT--TA-\
                 \---ATAATGCTA-CACAGTGTTC--GCTG-T-TTT-A-CTTCGG-GT-CCTCCTGCCCCTTGAGG-A\
                 \--TATGAT--A--CAC----AACCG--CT--CTCTAGACA-GGAAGAAA-ACTGCACCGGATA--"
        b2Ans = "AGCATTAGCTCACATAATAAGGATTGTGGGGCTACGGCATTTGTATGATAGCCTGACCGAGGTGGTAG\
                 \GGCATTAAGGTATCGTCGTCCTCCAGCTGGTGTTTTATCTGCCGCGCACCTCAAGTCGACTACGGTA\
                 \ACTATGGTTTATGCACGTGCAAACGTACTAACACTGGAGGCGGAATAATGATTAGAGTGGCCCGT"
        res2 = testAlign (medA :: String) (medB :: String)
        (a2', b2') = viewAlignment res2
        score2 = score res2

    it "first sequence"  $ a2'    `shouldBe` a2Ans
    it "second sequence" $ b2'    `shouldBe` b2Ans
    it "score"           $ score2 `shouldBe` (-168)

    let a3Ans = "ATCATCTCGGACTATGGTTGGGGATCTAAGACGTCCTCTTGGCTCTCGCC"
        b3Ans = "AAGAGATCGT--TATCGC-GCGGAT-TAAT---TCCGAGTCG-TC-CAC-"
        res3 = testAlign (smallA :: String) (smallB :: String)
        (a3', b3') = viewAlignment res3
        score3 = score res3

    it "first sequence"  $ a3'    `shouldBe` a3Ans
    it "second sequence" $ b3'    `shouldBe` b3Ans
    it "score"           $ score3 `shouldBe` (-16)

  describe "GlobalAlignment. AffineGap" $ do
    let testAlign = align (GlobalAlignment nuc44 (AffineGap (-5) (-1)))
        res1 = testAlign "AGGT" "CAAT"
        res2 = testAlign "AGGTACC" "CAATG"

    -- Честно посчитано руками в тетради
    it "score" $ score res1 `shouldBe` (-2)
    it "score" $ score res2 `shouldBe` (-8)
    it "alignment" $ viewAlignment res1 `shouldBe` ("--AGGT", "CAA--T")
    it "alignment" $ viewAlignment res2 `shouldBe` ("--AGGTACC","CAATG----")

  describe "GlobalAlignment. Sanity test" $ do
    let resBig1   = align (GlobalAlignment nuc44 (AffineGap (-5) (-5))) (bigA :: String) (bigB :: String)
        resBig2   = align (GlobalAlignment nuc44 (-5 :: Int)) (bigB :: String) (bigA :: String)

        resMed1   = align (GlobalAlignment nuc44 (AffineGap (-5) (-5))) (medA :: String) (medB :: String)
        resMed2   = align (GlobalAlignment nuc44 (-5 :: Int)) (medB :: String) (medA :: String)

        resSmall1 = align (GlobalAlignment nuc44 (AffineGap (-5) (-5))) (smallA :: String) (smallB :: String)
        resSmall2 = align (GlobalAlignment nuc44 (-5 :: Int)) (smallB :: String) (smallA :: String)
    it "SimpleGap is a special case of AffineGap (big)"   $ score resBig1 `shouldBe` score resBig2
    it "SimpleGap is a special case of AffineGap (med)"   $ score resMed1 `shouldBe` score resMed2
    it "SimpleGap is a special case of AffineGap (small)" $ score resSmall1 `shouldBe` score resSmall2


juliaSemiglobal :: Spec
juliaSemiglobal = do
  describe "Semiglobal Alignment, Simple Gap" $ do
    let testAlign = align (SemiglobalAlignment nuc44 (-3 :: Int))
    let resBig1   = testAlign (bigA :: String) (bigB :: String)
        resBig2   = testAlign (bigB :: String) (bigA :: String)
        resMed1   = testAlign (bigA :: String) (bigB :: String)
        resMed2   = testAlign (bigB :: String) (bigA :: String)
        resSmall1 = testAlign (bigA :: String) (bigB :: String)
        resSmall2 = testAlign (bigB :: String) (bigA :: String)

    it "Is symmetric (big): " $   score resBig1   `shouldBe` score resBig2
    it "Is symmetric (med): " $   score resMed1   `shouldBe` score resMed2
    it "Is symmetric (small): " $ score resSmall1 `shouldBe` score resSmall2
    -- scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=0, gap_extend=-3);
    -- pairalign(OverlapAlignment(), a, b, scoremodel)
    let aAns = "CTCTGCAAGCATAGAGATATTTCACCGGCAATATTTTTCGTTGAAGTGTATTTGTCCCTATTATCACACC\
                \AGTCATTCCGATGTCGTTAGGGCCGCCTGTTTTACTG--G-AG--TTC--C-CGTTCGATGTACCCTCC\
                \TCCATAGAGA-TTTAGGGGTA-T-AGACCGCGACGGACGGGGCTTGCGAAAGCTGCCCGATAATGGACT\
                \TTGACAATAGACGTAATCCGTGAACGGC-TGGCGTTTCACACTCAACCGCATTGAACAGCCATGCCATC\
                \CGC-AGGATGCGCCAGGGAA-AAGCCACCCCCAAG-AACGACT-CCAGG---GGCCCAAGTATGATGAT\
                \AGCGTTAGGCTTGCTCGAAGGCAAACAACTCCGAGGATGCTTAG-CAGGCCACAACGTGC--TCGATAG\
                \CACGAGCGTAAGATGCGAACGCAGTAGTACCTAATGCGGGCGAGACATTTCGGGCCGTTGAAGCCCTGT\
                \CTGTTTTTTCGTACAGAGGGTGGTTTGAAGTATCGCCC"
        bAns = "----------------------------------------------------------------------\
                \-----------------------------------TGAAGCAGCATTCGACTCACTC-ATGT-CGC-C-\
                \T--AGAGAGAATTTA----TAGTCA-A-----ATG-AC-----TTGTG---GC-GTTC--T--TGGACT\
                \TT-----TAGA-G----CCGA-AATG-CATGGAGG---A-A---A-C---ATT-AAG-G--ATGT-ATT\
                \TGTTAC--TGAGC-ATG-ATTAAGCAATTGCCAGGTAA--ACTGCCAGGTTTGGACCAA-TCT--TGGG\
                \AGC-TTAG--TACCTA-AAG--AAAGAAGTC--A---TGCACAGTCACGCCGCA-C---CAAT--ATA-\
                \CAC------------------------------------------------------------------\
                \--------------------------------------"
        res = testAlign (bigA :: String) (bigB :: String)
        (a', b') = viewAlignment res
        score1 = score res
    it "first sequence"  $ a'     `shouldBe` aAns
    it "second sequence" $ b'     `shouldBe` bAns
    it "score"           $ score1 `shouldBe` 358

    let a2Ans = "---------------------G-TT-TGGG-C-AT---ATTCA-A-GATCAGACCAATCGGTCGATGT\
                 \GAAA---CAAGTTAA--TAAT-G-C-TACAC-AG-TGTTCGCTGTTTTACTTCGGGTCCTCCTGCCC\
                 \CTTGAG--GA-TATGATACACAACCGCTCTCTA-G-ACAG-G-AA--G-A--AA-ACTGCAC-CGGA\
                 \-TA------------------"
        b2Ans = "AGCATTAGCTCACATAATAAGGATTGTGGGGCTACGGCATTTGTATGAT-AG-CC--T-GACCGAGGT\
                 \GGTAGGGCA--TTAAGGTA-TCGTCGTCCTCCAGCTG---G-TGTTTTA-T-CTG--CCGC--GCAC\
                 \CTCAAGTCGACTACGGTA-ACTATGG-T-T-TATGCAC-GTGCAAACGTACTAACACTGGAGGCGGA\
                 \ATAATGATTAGAGTGGCCCGT"
        res2 = testAlign (medA :: String) (medB :: String)
        (a2', b2') = viewAlignment res2
        score2 = score res2

    it "first sequence"  $ a2'    `shouldBe` a2Ans
    it "second sequence" $ b2'    `shouldBe` b2Ans
    it "score"           $ score2 `shouldBe` 266

    let a3Ans = "ATCATCTCG--GACTATGGTT---G-G-GGATCTAA----GA--CGTCCTCTTGGCTCTCGCC"
        b3Ans = "---------AAGAG-ATCGTTATCGCGCGGAT-TAATTCCGAGTCGTCCAC------------"
        res3 = testAlign (smallA :: String) (smallB :: String)
        (a3', b3') = viewAlignment res3
        score3 = score res3

    it "first sequence"  $ a3'    `shouldBe` a3Ans
    it "second sequence" $ b3'    `shouldBe` b3Ans
    it "score"           $ score3 `shouldBe` 63
  describe "SemiglobalAlignment. Sanity test" $ do
    let resBig1   = align (SemiglobalAlignment nuc44 (AffineGap (-3) (-3))) (bigA :: String) (bigB :: String)
        resBig2   = align (SemiglobalAlignment nuc44 (-3 :: Int)) (bigB :: String) (bigA :: String)

        resMed1   = align (SemiglobalAlignment nuc44 (AffineGap (-3) (-3))) (medA :: String) (medB :: String)
        resMed2   = align (SemiglobalAlignment nuc44 (-3 :: Int)) (medB :: String) (medA :: String)

        resSmall1 = align (SemiglobalAlignment nuc44 (AffineGap (-3) (-3))) (smallA :: String) (smallB :: String)
        resSmall2 = align (SemiglobalAlignment nuc44 (-3 :: Int)) (smallB :: String) (smallA :: String)
    it "SimpleGap is a special case of AffineGap (big)"   $ score resBig1 `shouldBe` score resBig2
    it "SimpleGap is a special case of AffineGap (med)"   $ score resMed1 `shouldBe` score resMed2
    it "SimpleGap is a special case of AffineGap (small)" $ score resSmall1 `shouldBe` score resSmall2

juliaLocal :: Spec
juliaLocal = do
  describe "Local Alignment, Simple Gap" $ do
    let testAlign = align (LocalAlignment nuc44 (-7 :: Int))
    let resBig1   = testAlign (bigA :: String) (bigB :: String)
        resBig2   = testAlign (bigB :: String) (bigA :: String)
        resMed1   = testAlign (bigA :: String) (bigB :: String)
        resMed2   = testAlign (bigB :: String) (bigA :: String)
        resSmall1 = testAlign (bigA :: String) (bigB :: String)
        resSmall2 = testAlign (bigB :: String) (bigA :: String)

    it "Is symmetric (big): " $   score resBig1   `shouldBe` score resBig2
    it "Is symmetric (med): " $   score resMed1   `shouldBe` score resMed2
    it "Is symmetric (small): " $ score resSmall1 `shouldBe` score resSmall2
    -- scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=0, gap_extend=-7);
    -- pairalign(LocalAlignment(), a, b, scoremodel)
    let aAns = "TGAA-CAGCCATGCCA-TC-CGCAGGATGCGCC-AGGGAAAAGCCACCCCCAAGAACGACTCCA\
                \GGGGCCCAAGTATGATGATAGC-GTTAGGCTTGCTCGAAGGCAAACAACTCCGAGG-ATGC-T\
                \TAGCAGGCCACAACGTGCTCGATAGCACGAGCGTAAGATGCGAACGCAGTAGTACCTAATGCGGGCGAGACATT"
        bAns = "TGAAGCAGC-ATTCGACTCACTCATG-T-CGCCTAGAGAGAATTTATAGTCAA-AT-GACTTGT\
                \GGCGTTCTTGGACTTTTAGAGCCGAAATGCATG---GA-GG-AAACATTAAGGATGTATTTGT\
                \TA-CTGAGCATGAT-TAAGCAATTGC-C-AG-GTAAACTGCCAG-GTT-TGG-ACC-AATCTTGG-GAG-C-TT"
        res = testAlign (bigA :: String) (bigB :: String)
        (a', b') = viewAlignment res
        score1 = score res
    it "first sequence"  $ a'     `shouldBe` aAns
    it "second sequence" $ b'     `shouldBe` bAns
    it "score"           $ score1 `shouldBe` 108

    let a2Ans = "GGCATATTCAAGATCAGACCAATCG-GTCGAT-GTGAAACAAGTTAATAATGCTACACAGTGT\
                 \TCGCTGTTTTA-CTTCGG-GT-CCTCCTGCCCCTTGAGG-A--TATGATACACAACCGCTCT\
                 \CTAGACAGGAAGAAAACTGCAC-CGGA-TA"
        b2Ans = "GGCATTTGTATGAT-AGCCTGACCGAGGTGGTAGGGCATTAAGGTA-TCGTCGTCCTCCA-GC\
                 \TGG-TGTTTTATCTGCCGCGCACCTCAAGTCGACTACGGTAACTATGGTTTATGCACG-TG-\
                 \CAA-AC-GTACTAACACTGGAGGCGGAATA"
        res2 = testAlign (medA :: String) (medB :: String)
        (a2', b2') = viewAlignment res2
        score2 = score res2

    it "first sequence"  $ a2'    `shouldBe` a2Ans
    it "second sequence" $ b2'    `shouldBe` b2Ans
    it "score"           $ score2 `shouldBe` 91

    let a3Ans = "TCTAAGACGTCCTC"
        b3Ans = "TCCGAGTCGTCCAC"
        res3 = testAlign (smallA :: String) (smallB :: String)
        (a3', b3') = viewAlignment res3
        score3 = score res3

    it "first sequence"  $ a3'    `shouldBe` a3Ans
    it "second sequence" $ b3'    `shouldBe` b3Ans
    it "score"           $ score3 `shouldBe` 34

  describe "LocalAlignment. Sanity test" $ do
    let resBig1   = align (LocalAlignment nuc44 (AffineGap (-7) (-7))) (bigA :: String) (bigB :: String)
        resBig2   = align (LocalAlignment nuc44 (-7 :: Int)) (bigB :: String) (bigA :: String)

        resMed1   = align (LocalAlignment nuc44 (AffineGap (-7) (-7))) (medA :: String) (medB :: String)
        resMed2   = align (LocalAlignment nuc44 (-7 :: Int)) (medB :: String) (medA :: String)

        resSmall1 = align (LocalAlignment nuc44 (AffineGap (-7) (-7))) (smallA :: String) (smallB :: String)
        resSmall2 = align (LocalAlignment nuc44 (-7 :: Int)) (smallB :: String) (smallA :: String)
    it "SimpleGap is a special case of AffineGap (big)"   $ score resBig1 `shouldBe` score resBig2
    it "SimpleGap is a special case of AffineGap (med)"   $ score resMed1 `shouldBe` score resMed2
    it "SimpleGap is a special case of AffineGap (small)" $ score resSmall1 `shouldBe` score resSmall2

bigA :: String
bigA = "CTCTGCAAGCATAGAGATATTTCACCGGCAATATTTTTCGTTGAAGTGTATTTGTCCCTATTATCACACCAG\
        \TCATTCCGATGTCGTTAGGGCCGCCTGTTTTACTGGAGTTCCCGTTCGATGTACCCTCCTCCATAGAGATT\
        \TAGGGGTATAGACCGCGACGGACGGGGCTTGCGAAAGCTGCCCGATAATGGACTTTGACAATAGACGTAAT\
        \CCGTGAACGGCTGGCGTTTCACACTCAACCGCATTGAACAGCCATGCCATCCGCAGGATGCGCCAGGGAAA\
        \AGCCACCCCCAAGAACGACTCCAGGGGCCCAAGTATGATGATAGCGTTAGGCTTGCTCGAAGGCAAACAAC\
        \TCCGAGGATGCTTAGCAGGCCACAACGTGCTCGATAGCACGAGCGTAAGATGCGAACGCAGTAGTACCTAA\
        \TGCGGGCGAGACATTTCGGGCCGTTGAAGCCCTGTCTGTTTTTTCGTACAGAGGGTGGTTTGAAGTATCGCCC"

bigB :: String
bigB = "TGAAGCAGCATTCGACTCACTCATGTCGCCTAGAGAGAATTTATAGTCAAATGACTTGTGGCGTTCTTGGACT\
        \TTTAGAGCCGAAATGCATGGAGGAAACATTAAGGATGTATTTGTTACTGAGCATGATTAAGCAATTGCCAGG\
        \TAAACTGCCAGGTTTGGACCAATCTTGGGAGCTTAGTACCTAAAGAAAGAAGTCATGCACAGTCACGCCGCACCAATATACAC"

medA :: String
medA = "GTTTGGGCATATTCAAGATCAGACCAATCGGTCGATGTGAAACAAGTTAATAATGCTACACAGTGTTCGCTG\
        \TTTTACTTCGGGTCCTCCTGCCCCTTGAGGATATGATACACAACCGCTCTCTAGACAGGAAGAAAACTGCACCGGATA"

medB :: String
medB = "AGCATTAGCTCACATAATAAGGATTGTGGGGCTACGGCATTTGTATGATAGCCTGACCGAGGTGGTAGGGCA\
        \TTAAGGTATCGTCGTCCTCCAGCTGGTGTTTTATCTGCCGCGCACCTCAAGTCGACTACGGTAACTATGGT\
        \TTATGCACGTGCAAACGTACTAACACTGGAGGCGGAATAATGATTAGAGTGGCCCGT"

smallA :: String
smallA = "ATCATCTCGGACTATGGTTGGGGATCTAAGACGTCCTCTTGGCTCTCGCC"

smallB :: String
smallB = "AAGAGATCGTTATCGCGCGGATTAATTCCGAGTCGTCCAC"
