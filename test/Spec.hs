import           Test.Hspec

main :: IO ()
main = hspec $
    describe "Future tests" $
      it "should pass" $
        True `shouldBe` True
