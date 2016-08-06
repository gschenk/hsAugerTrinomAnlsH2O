-- multinomial analysis that considers Auger
-- version 2, trinomial analysis, adds capture probs to input
-- for H2O with states 2a1, 1b2, 3a1, 1b1
import Data.List
import Text.Printf
import Data.Time
--import Data.Array.Parallel  --doesn't work with base 4.8 yet

-------- ======== data type definitions ======== --------
data TriInds = TrI Int [Int] [Int] deriving (Show, Eq)

-------- ======== global variable ======== --------
bA    = True :: Bool  -- toggle Auger analysis
excdUt = False  :: Bool -- error when probability exceeds unity, or set to 1, negative probs set to zero when small
corrTld = 1e-7  -- depending on excdUt, correct when error smaller than this, else throw exception

-- example input lines for testing
--ds = [20.0,0.4,0.1269480396,5.143592415e-2,0.2023602794,9.899421803e-2,5.018767105e-3,4.689024479e-2,0.3519103128,6.847626998e-2,2.598374102e-2,0.7829925189]
--ps =  (sumOrbOccp . drop 2) ds


-------- ======== input related ======== --------

-- check if p is in valid range 0<p<=1
pRngChk :: Double -> Double
pRngChk p
    | p < 0               = error "negative probability"
    | excdUt && p > 1     = error "probability exceeds unity"
    | not excdUt && p > 1  = 1
    | p-1 >= corrTld =  error "probability exceeds unity, considerably"
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
prpPs :: [Double] -> [Double]
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
pNoA :: (Integral a, Fractional b) => a -> a -> b -> b
pNoA 0 _ _ = 1
pNoA _ _ _ = 0


-------- ======== Multinomial Analysis (Eq. 3) ======== --------
-- permutate indices, step one,

-- creates a list of `n` element lists that contain indices i=0,1,2
--   such that all combinations that sum up to `q` are included
--   prmr function recursivly constructs rules to build these lists
--   bldr function builds the actual lists
--   when q>n false lists are built, that have to be filtered out again
prmSeeds :: Integral a => Int -> Int -> [[a]]
prmSeeds n = filter (\xs -> n == length xs) . map bldr .  prmr n
    where   prmr 0 _ = []
            prmr 1 0 = [(1,0,0)]
            prmr 1 1 = [(0,1,0)]
            prmr l 0 = [(l,0,0)]
            prmr l 1 = [(l-1,1,0)]
            prmr l q
                | 2*l > q = (:) (l-q,q,0) . map f .  prmr l $ (q-2)
                | otherwise  = error "q must not exceed 2*l in prmsds"
                where f (a,b,c) = (a-1,b,c+1)
            bldr (k,l,m) = (take k) zs ++ (take l) os ++ (take m) ts
                where zs = repeat 0 -- zeros
                      os = repeat 1 -- ones
                      ts = repeat 2 -- twos


-- permutate indices step two for Binomial analysis (only!)
-- permutate, filter impossible configurations, remove duplicates
prmBiInds :: Integral a => Int -> [[a]]
prmBiInds q
    | q <= 8 = (concat . map (nub . filter f . permutations) . prmSeeds l) q
    | otherwise = error "prmBiInds"
    where l = 5 -- length of index list
          f (a:n:_) = n >= a  -- filter out permutations where a exceeds n
          g as = length as == l



-- creates and filters a lit of all possible permutations of indices that lead to
-- k electrons at the projectile and l electrons in the continuum
-- the elements of that list have a structure of (a,ks,ls)
-- where 'a' is the number of Auger electron, ks the list of indices for capture
-- and 'l' the list for e in the continuum; these lists have n=9 elements
triInds :: Int -> Int -> [TriInds]
triInds k l  = (map toTrI . nub . filters . map strtr . concat . map permutations . prmSeeds n) (k+l)
    where   n = 9
            filters = filter f3 . filter f2 . filter f1
            strtr (a:kls) = (a,take 4 kls,drop 4 kls)
            f1 (_,ks,_) = sum ks == k  -- fits count of captured electrons
            f2 (_,ks,ls) = foldl1 (&&) . map (<=2) $ zipWith (+) ks ls  -- remove only 2 or less electrons per orbital
            f3 (a,ks,ls) = (head ks) + (head ls) >= a  -- Auger electrons must not exceed number of holes
            toTrI (a,ks,ls) = TrI a ks ls



-------- ======== prepare associated list of indices ======== --------


-- returns associated list with octal key
klTriInds :: [(Int, [TriInds])]
klTriInds = map (\x -> (f x,  triInds (fst x) (snd x) ) )   klTuples
    where f (a,b) = 10*a+b
          g (a,ks,ls) = TrI a ks ls

-- makes a list of (k,l) tuples
klTuples :: [(Int, Int)]
klTuples =  concat $ map (\x -> (map (\y -> (x,y)) [0..8-x])) [0..8]

-- lookup, error when nothing is returned
myLookup :: Eq a => a -> [(a, b)] -> b
myLookup  x xs = case lookup x xs of
    Just y -> y
    Nothing -> error "Lookup failed"


-------- ======== Trinomial ======== --------

-- Input in the form
--   probabilities, a list of 8 single particle probabilities
--     [os]++[cs]
--     os, four, occupation of target orbitals
--     cs, four, capture probabilities
--   nine integers, each in 0,1,2
--     organised in a 3-tuple,
--     first, Auger index a
--     second, ks, four, number of electrons at the projectile
--     third, ls, four, number of electrons in the continuum
--trinProd :: (Ord a, Fractional a, Integral b) => [b] -> [b] -> [a] -> [a] -> a
trinProd ks ls os cs = tc * (f cs ks) * (f is ls) * (f os ns)
    where is = map (1-) $ zipWith (+) os cs  -- ionisation probability
          ns = map (2-) $ zipWith (+) ks ls
          tc = product $ zipWith g ks ls
          g 0 1 = 2 --trinomial coefficients
          g 1 0 = 2  -- only a few trivial cases needed
          g 1 1 = 2
          g _ _ = 1
          f rs is  -- powers (0,1,2) of probs, checks input
            | minimum is < 0   = error "trinProd, negative index"
            | maximum is > 2   = error "trinProd, forbidden index"
            | excdUt && minimum rs < 0   = error "trinProd, negative prob."
            | (minimum rs)*(-1) >= corrTld = error "trinProd, exceedingly negative prob."
            | maximum rs > 1.0 = error "trinProd, probability > 1"
            |otherwise = product $ zipWith (^^) (map (\x -> if x < 0 then 0 else x) rs) is
            -- negative probabilities are corrected to zero
            -- probs smaller than -corrTld cause an exception still


-- sum of trinomial products
--   that lead to configuration k, l
--   takes as argument one set of probabilities
--   in the form (os++cs), where os ar four single-particle
--   occupation probabilities and cs capture probabilities
--   the sum is realised by calling a function that
--   generates a list with all possible permutations of
--   9 indices, each in {0,1,2} where the first is
--   the number of Auger electrons a, the rest
--   are lists of four, ks and ls, capture and ionisation
--   probabilities respectively.
trinSum :: [Double] -> Int -> Int -> Double
trinSum ocs k l  =  sum $ map f is
    where is =  myLookup (10*k+l) klTriInds
          os = take 4 ocs
          cs = drop 4 ocs
          f (TrI a ks ls) =(*) (gA a n m) $ trinProd ks ls os cs
            where n = head ks + head ls
                  m = fromIntegral . sum $ (tail ks)++(tail ls)
          gA a n m -- Auger correction factor
             | not bA && a == 0 = 1
             | not bA && a /= 0 = 0
             | otherwise        = pA a n m


-- gets trinomial probabilities for all k, l combinations
klTrinProbs :: [Double] -> [Double]
klTrinProbs ps = concat $ map (\x-> (map (trinSum ps x) [0..8-x])) [0..8]


-------- ======== Binomial ======== --------
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
--niceOut :: (Num a, Floating b) => [a] -> b -> [b] -> String
niceOut :: (PrintfArg a, PrintfArg b) => [a] -> b -> [b] -> String
niceOut [keV,b] pn = (unwords .  (se:) . (sb:) . (++ [spn]) . map (printf "%14.7e") )
    where se  = printf  "%6.1f" keV; sb = printf  "%6.2f" b
          spn = printf "%14.7e" pn

-- function usefull for data-file comments/headers,
-- only for interactive use
klStr =  concat $ map (\x-> (map (f x) [0..8-x])) [0..8]
    where f a b = show (a,b)


-------- ========        main      ======== --------
-- pipe input data to this script
-- use `grep -Ev "^$|^#"` to filter comments from input
main =  do
    line <- getLine
    if null line
        then return ()
        else do
            let rawData = (map read .words) line :: [Double]
            let eb = take 2 rawData
            let ps = prpPs . drop 2 $ rawData
            let os = take 4 ps
            let pnet = ((*2) . (4-) . sum ) os
            let pkl = klTrinProbs ps
            putStrLn $ niceOut eb pnet pkl
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
