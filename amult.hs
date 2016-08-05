-- multinomial analysis that considers Auger
-- for H2O with states 2a1, 1b2, 3a1, 1b1
import Data.List
import Text.Printf

-------- ======== global variable ======== --------
bA   = False :: Bool  -- toggle Auger analysis

-- example input lines for testing
--ds = [20.0,0.4,0.1269480396,5.143592415e-2,0.2023602794,9.899421803e-2,5.018767105e-3,4.689024479e-2,0.3519103128,6.847626998e-2,2.598374102e-2,0.7829925189]
--ps =  (sumOrbOccp . drop 2) ds


-------- ======== input related ======== --------

-- check if p is in valid range 0<p<=1
pRngChk :: (Num a, Ord a) => a -> a
pRngChk p
    | p < 0    =error "negative probability"
    | p > 1    =error "probability exceeds unity"
    | p == 0   =error "zero probability not considered"
    | otherwise = p


-- locate and sum occupation of each initial condition
sumOrbOccp :: (Num b) => [b] -> [b]
sumOrbOccp os  = map ($os) [o2a1,o1b2,o3a1,o1b1]
    where o2a1  = (f 3  1)  -- hard coded column positions
          o1b2  = (f 3  4)
          o3a1  = (f 3  7)
          o1b1  = (f 1 10)
          -- get  n elements starting at the m-th position of a list, sum
          f n m = (sum . take n .snd . splitAt (m-1))


-- prepare list of target occupation ps
-- check if probabilities are valid
prpPs :: (Num a, Ord a) => [a] -> [a]
prpPs =  map pRngChk . sumOrbOccp . map pRngChk 


-------- ======== Auger probabilities ======== --------

-- mathcal P (a,n,m) [for bA = True]
pA :: (Integral a, Fractional b) => a -> a -> b -> b
pA 0 0 _ =1
pA 0 1 m =m
pA 1 1 m =1-m/6
pA 0 2 m =1/36*m^2
pA 1 2 m =1/18*(6*m-m^2)
pA 2 2 m =1/36*(6-m)^2
pA _ _ _ =0


-- no Auger probs, bA = False
--pNoA :: (Integral a, Fractional b) => a -> a -> b -> b
pNoA 0 _ _ = 1
pNoA _ _ _ = 0


-------- ======== Multinomial Analysis (Eq. 3) ======== --------
-- permutate indices, step one,
-- make lists with correct number of zeros, ones, and twos
prmSeeds :: (Eq a, Num t, Num a) => a -> [[t]]
prmSeeds q  | q==0  = f [[0,0,0]]
            | q==1  = f [[1,1,0]]
            | q==2  = f [[2,2,0],[1,0,1]]
            | q==3  = f [[3,3,0],[2,1,1]]
            | q==4  = f [[4,4,0],[3,2,1],[2,0,2]]
            | q==5  = f [[5,5,0],[4,3,1],[3,1,2]]
            | q==6  = f [        [5,4,1],[4,2,2],[3,0,3]]
            | q==7  = f [                [5,3,2],[4,1,3]]
            | q==8  = f [                        [5,2,3],[4,0,4]]
            | otherwise = error "q must not exceed 8"
            where zs = repeat 0 -- empty list
                  os = repeat 1
                  ts = repeat 2
                  f = map g
                  g [k,l,m] = (take (5-k)) zs ++ (take l) os ++ (take m) ts


-- permutate indices step two
-- permutate, filter impossible configurations, remove duplicates
prmInds :: (Integral a) => a -> [[Int]]
prmInds  = concat . map (nub . filter f . permutations) . prmSeeds
    where f (a:n:_) = n >= a  -- filter out permutations where a exceeds n


multProd q (a:n:ms)
    | q == sum (a:n:ms) = (*) (fA a n m) . prod
    | otherwise = \xs -> 0
    where   nms = n:ms
            m = (fromIntegral . sum) ms
            fA | bA        = pA    -- Auger prob. function
               | otherwise = pNoA
            prod = product . (zipWith f nms)
                where f 0 = \x -> x^^2   -- where f is π in Eq. (3)
                      f 1 = \x -> 2*(1-x)*x
                      f 2 = \x -> (1-x)^^2
                      f _ = \x -> 0


-- for every q use a list containing permutations of all
-- electron configurations to get the correct products
-- and sum them up
multSum :: (Eq c, Fractional c) => Int -> [c] -> c
multSum q ps = (sum . map (`f` ps)) is
    where f  = multProd q
          is = prmInds q


-------- ======== Formating Output ======== --------
--nceOt :: (Num a) => (a,a) -> [a] -> String
niceOut [keV,b] = (unwords .  (se:) . (sb:) . map (printf "%14.7e") )
    where se = printf  "%6.1f" keV; sb = printf  "%6.2f" b



-------- ========        main      ======== --------
-- pipe input data to this script
-- use `grep -Ev "^$|^#"` to filter comments from input
main =  do
    line <- getLine
    if null line
        then return ()
        else do
            let rawData = (map read .words) line :: [Double]
            let eb =take 2 rawData
            let ps = (prpPs . tail . tail) rawData
            (putStrLn . niceOut eb . map (`multSum`  ps)) [1..8]
            main


{- Input files ought to be organised like this example:
# 1     2       3                   4                   5                   6                   7                   8                   9                   10                  11                  12
#>> param. <<   >> transition probabilities                                                                                                                                                                       <<
# E     b       >> initial condition 2a1 orbital occupied             <<    >> 1b2 orbital                                        <<    >> 3a1 orbital                                        <<    >> 1b1 orbtl. <<
#[keV]  [au]    2a1->2a1            2a1->1b2            2a1->3a1            1b2->2a1            1b2->1b2            1b2->3a1            3a1->2a1            3a1->1b2            3a1->3a1            1b1->1b1
 0020    0.40   0.1269480396e+00    0.5143592415e-01    0.2023602794e+00    0.9899421803e-01    0.5018767105e-02    0.4689024479e-01    0.3519103128e+00    0.6847626998e-01    0.2598374102e-01    0.7829925189e+00
 0020    0.60   0.1648643287e+00    0.5618833391e-01    0.1706703322e+00    0.1480957416e+00    0.2789072993e-01    0.6727443586e-01    0.3493231849e+00    0.6161736900e-01    0.6928626147e-01    0.8186759488e+00

where the first two rows contain energy and collision
parameter, the consecutive columns probabilities 
to occupy a specific state (0<p<=1), sorted in set of
3, 3, 3 and 1.  Comments ought to be stripped off,
eg with grep, before piping to stdin. 
-}

-- vim: set ts=4 sw=4 sts=4 nowrap et ai:
