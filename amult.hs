-- multinomial analysis that considers Auger
-- version 2, trinomial analysis, adds capture probs to input
-- for H2O with states 2a1, 1b2, 3a1, 1b1
import Data.List
import Text.Printf

-------- ======== global variable ======== --------
bA    = True :: Bool  -- toggle Auger analysis
excdUt = False  :: Bool -- error when probability exceeds unity, or set to 1

-- example input lines for testing
--ds = [20.0,0.4,0.1269480396,5.143592415e-2,0.2023602794,9.899421803e-2,5.018767105e-3,4.689024479e-2,0.3519103128,6.847626998e-2,2.598374102e-2,0.7829925189]
--ps =  (sumOrbOccp . drop 2) ds


-------- ======== input related ======== --------

-- check if p is in valid range 0<p<=1
pRngChk :: (Num a, Ord a) => a -> a
pRngChk p
    | p < 0               = error "negative probability"
    | excdUt && p > 1     = error "probability exceeds unity"
    | not excdUt && p > 1 = 1
    | otherwise           = p


-- locate and sum occupation of each initial condition
-- followed by capture, same order
sumOrbOccp :: (Num b) => [b] -> [b]
sumOrbOccp os  = map ($os) [o2a1,o1b2,o3a1,o1b1,c2a1,c1b2,c3a1,c1b1]
    where o2a1  = (f 3  1)  -- hard coded column positions
          o1b2  = (f 3  4)
          o3a1  = (f 3  7)
          o1b1  = (f 1 10)
          c2a1  = (f 1 11)
          c1b2  = (f 1 12)
          c3a1  = (f 1 13)
          c1b1  = (f 1 14)
          -- get  n elements starting at the m-th position of a list, sum
          f n m = (sum . take n .snd . splitAt (m-1))


-- prepare list of target occupation ps
-- and capture probabilities
-- check if probabilities are valid
prpPs :: (Num a, Ord a) => [a] -> [a]
prpPs =  map pRngChk . sumOrbOccp . map pRngChk


-------- ======== Auger probabilities ======== --------

-- mathcal P (a,n,m) [for bA = True]
pA :: (Integral a, Fractional b) => a -> a -> b -> b
pA 0 0 _ =1
pA 0 1 m =m/6
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


-- seeds to build list of numbers to permutate,
-- of length l, distributing q electrons
-- where the 3-tuple shows the numbers of
-- recursive
-- (zeros,ones,twos) in the list
prmsds :: (Integral b) => b ->b -> [(b,b,b)]
prmsds 0 _ = []
prmsds 1 0 = [(1,0,0)]
prmsds 1 1 = [(0,1,0)]
prmsds l 0 = [(l,0,0)]
prmsds l 1 = [(l-1,1,0)]
prmsds l q
    | 2*l > q = (:) (l-q,q,0) . map f .  prmsds l $ (q-2)
    | otherwise  = error "q must not exceed 2*l in prmsds"
    where f (a,b,c) = (a-1,b,c+1)


--  builds seeds for permutation
--  builds impossible values too, have to be filtered out
prmBL (k,l,m) = (take k) zs ++ (take l) os ++ (take m) ts
    where zs = repeat 0 -- empty list
          os = repeat 1
          ts = repeat 2


-- permutate indices step two
--  for Binomial analysis only!
-- permutate, filter impossible configurations, remove duplicates
--prmBiInds :: (Integral a) => a -> [[Int]]
prmBiInds q
    | q <= 8 = (concat . map (nub . filter f . permutations) . filter g . map prmBL . prmsds l) q
    | otherwise = error "prmBiInds"
    where l = 5 -- length of index list
          f (a:n:_) = n >= a  -- filter out permutations where a exceeds n
          g as = length as == l



-- lists of the form [[c,i,o],...] are needed for trinomial
-- this forms it form intput lists of form [o...,c...]
trinProbList ocs
    |length ocs /= 8 = error "trinProdList: list lenght not 8"
    |otherwise = zipWith f os cs
    where os = take 4 ocs
          cs = drop 4 ocs
          f o c
            | o+c > 1.0 = error "trinProbList: negative ionization probability"
            |otherwise    = [c, 1.0-o-c, o]


trinIndList kls
    |length kls /= 8 = error "trinIndList: list lenght not 8"
    |otherwise  = zipWith f ks ls
    where ks = take 4 kls
          ls = drop 4 kls
          f k l
            | 2-k-l < 0 = error "trinIndList: negative index m"
            | otherwise = [k, l, 2-k-l]

-- compared to binomial, it requires twice the number of orbital
-- structure index list (a:k2a1:k1b2:...:k1b1:l2a1:l1a1:...)
-- indices
trinProd k l (a:kls)
    | k == sum ks && l == sum (a:ls)  = (*) (fA a nkl m')  . trinOrbProd iss . trinProbList
    | otherwise = \xs -> 0
    where ks = take 4 kls  -- indices corresponding to captured electrons
          ls = drop 4 kls
          nkl = (+) (head ks) (head ls)  -- number of 2a1 electron that are removed
          m'   = fromIntegral $ sum ((tail ks) ++ (tail ls)) -- number of electrons lost from outer orbitals
          iss = trinIndList kls
          fA | bA        = pA    -- Auger prob. function
             | otherwise = pNoA


--trinOrbProd :: (Fractional c, Integral b) => [[b]] -> [[c]] -> c
trinOrbProd klms cios   = (product . (zipWith f klms)) cios
    where f ps is
            | length ps /= 3  = error "trinOrbProd: wrong list lenght"
            | length is /= 3  = error "trinOrbProd: wrong list lenght"
            | otherwise       = (product . zipWith (^^) is) ps
            -- lists of three elements, each required
            -- calculates c^k * i^l * o^m


binProd q (a:n:ms)
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
binSum :: (Eq c, Fractional c) => Int -> [c] -> c
binSum q ps = (sum . map (`f` ps)) is
    where f  = binProd q
          is = prmBiInds q

-- net removal probabilty calculated from sum q*p_q
netRecPQ  = \x -> x++f x
    where f    =  (:[]) . sum . zipWith (*) [1..]

-- compare both results for pnet
compPnet pn pqs | (f pn pqs) < (2*eps) = pqs
                | otherwise            = error "failed"
                where f y = abs . (-) y . last
                      eps = 1.11e-16

-------- ======== Formating Output ======== --------
--nceOt :: (Num a) => (a,a) -> [a] -> String
niceOut [keV,b] pn = (unwords .  (se:) . (sb:) . (++ [spn]) . map (printf "%14.7e") )
    where se  = printf  "%6.1f" keV; sb = printf  "%6.2f" b
          spn = printf "%14.7e" pn



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
            let os = (take 4 . prpPs . tail . tail) rawData
            let pnet = ((*2) . (4-) . sum ) os
            (putStrLn . niceOut eb pnet . netRecPQ . map (`binSum`  os)) [1..8]
            main


{- Input files ought to be organised like this example:
#1      2       3                   4                   5                   6                   7                   8                   9                   10                  11                  12                  13                 14                 15                 16
#               >> inito 1 (2a1)                                      <<    >> inito 2 (1b2)                                      <<    >> inito 3 (3a1)                                      <<    >> i. 4 (1b1) <<    >> capture into projectile, originating from orbital:                  <<
#ELAB   B       2a1->2a1            2a1->1b2            2a1->3a1            1b2->2a1            1b2->1b2            1b2->3a1            3a1->2a1            3a1->1b2            3a1->3a1            1b1->1b1            2a1                1b2                3b1                1b1
 0020    0.40   0.1269480396e+00    0.5143592415e-01    0.2023602794e+00    0.9899421803e-01    0.5018767105e-02    0.4689024479e-01    0.3519103128e+00    0.6847626998e-01    0.2598374102e-01    0.7829925189e+00    0.3450854978D+00   0.5630833588D+00   0.3127868638D+00   0.1439064381D+00
 0020    0.60   0.1648643287e+00    0.5618833391e-01    0.1706703322e+00    0.1480957416e+00    0.2789072993e-01    0.6727443586e-01    0.3493231849e+00    0.6161736900e-01    0.6928626147e-01    0.8186759488e+00    0.3306837587D+00   0.5359625416D+00   0.2898823084D+00   0.1354759864D+00

where the first two rows contain energy and collision
parameter, the consecutive columns probabilities
to occupy a specific state (0<p<=1), sorted in set of
3, 3, 3 and 1.  Followed by capture probabilities, for
4 orbitals.

Comments ought to be stripped off,
eg with grep, before piping to stdin.
-}

-- vim: set ts=4 sw=4 sts=4 nowrap et ai:
