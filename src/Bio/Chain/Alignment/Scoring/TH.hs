{-# LANGUAGE CPP             #-}
{-# LANGUAGE TemplateHaskell #-}
module Bio.Chain.Alignment.Scoring.TH where

import Data.Char                 (toLower, toUpper)
import Language.Haskell.TH
import Language.Haskell.TH.Quote

import Bio.Chain.Alignment.Scoring.Loader

type Substitution a = a -> a -> Int

matrix :: QuasiQuoter
matrix = QuasiQuoter { quotePat = undefined
                     , quoteType = undefined
                     , quoteExp = undefined
                     , quoteDec = matrixDec
                     }


matrixDec :: String -> Q [Dec]
matrixDec s = do let slines = lines s
                 let txt = (unlines . tail) slines
                 let name = head slines
                 (typeName, dataDecl) <- typeDec name
                 (funcName, funDecls) <- functionDec name txt
                 smDecl <- instDec typeName funcName
                 return $ dataDecl : smDecl : funDecls

instDec :: Name -> Name -> Q Dec
instDec typeN funN = return $ decl [body]
  where decl = InstanceD Nothing [] (AppT (ConT ''ScoringMatrix) (ConT typeN))
        body = FunD 'scoring [Clause [WildP] (NormalB (VarE funN)) []]

typeDec :: String -> Q (Name, Dec)
typeDec name = do typeN <- newName name
                  let dataN = mkName (nameBase typeN)
#if !MIN_VERSION_template_haskell(2,12,0)
                  let dervs = [ConT ''Show, ConT ''Eq]
#else
                  let dervs = [DerivClause Nothing [ConT ''Show, ConT ''Eq]]
#endif
                  return (typeN, DataD [] typeN [] Nothing [NormalC dataN []] dervs)

functionDec :: String -> String -> Q (Name, [Dec])
functionDec name txt = do let subM = loadMatrix txt
                          funName <- newName (toLower <$> name)
                          let funSign = SigD funName (AppT (ConT ''Substitution) (ConT ''Char))
                          let clauses = concatMap mkClause subM
                          let funDecl = FunD funName clauses
                          return (funName, [funSign, funDecl])

mkClause :: ((Char, Char), Int) -> [Clause]
mkClause ((c, d), i) = fmap (\pair -> Clause pair (NormalB (litI i)) []) casePairs
  where 
    litC = LitP . CharL
    litI = LitE . IntegerL . fromIntegral

    casingFunctions :: [Char -> Pat]
    casingFunctions = [litC . toUpper, litC . toLower]

    casePairs :: [[Pat]]
    casePairs = [ [f c, g d] | f <- casingFunctions, g <- casingFunctions ]

