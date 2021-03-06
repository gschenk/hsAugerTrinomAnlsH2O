-- this programme performs a multinomial analysis
-- for H_2O that considers Auger processes
--
-- version 2, trinomial analysis, adds capture probs to input
-- for H2O with states 2a1, 1b2, 3a1, 1b1
--
-- This is a Haskell source code that works correctly
-- with The Glorious Glasgow Haskell Compilation System,
-- version 8.10
--
-- Use cabal or stack to build
--
-- This code is released to the public domain under a
-- Creative Commons CC0 1.0 dedication
-- gerald schenk 2016
-- gschenk@yorku.ca (perm email: gschenk@gmx.net)
--
import           Data.List
import           Data.Time
import           Text.Printf
--import Data.Array.Parallel  --doesn't work with base 4.8 yet


-------- ======== data type definitions ======== --------
data TriInds = TrI Int [Int] [Int] deriving (Show, Eq)

-------- ======== global variable ======== --------
bA    = True :: Bool  -- toggle Auger analysis: True on, False off
excdUt = False  :: Bool -- error when probability exceeds unity, or set to 1, negative probs set to zero when small
corrTld = 1e-7  -- depending on excdUt, correct when error smaller than this, else throw exception
epsilon = 1e-12


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
inptMask :: (Num a, Eq a) => [a] -> [a]
inptMask os  = (`map` fs) ($os)
  where fs    = [o2a1,o1b2,o3a1,o1b1,c2a1,c1b2,c3a1,c1b1]
        f 1 m = (!!m)
        f n m = sum . take n . drop m -- sum n elements, start at m-th pos
        o2a1  = f 3  0  -- hard coded column positions m and lengths n
        o1b2  = f 3  3
        o3a1  = f 3  6
        o1b1  = f 1  9
        c2a1  = f 1 10
        c1b2  = f 1 11
        c3a1  = f 1 12
        c1b1  = f 1 13


-- prepare list of target occupation ps
-- and capture probabilities
-- check if probabilities are valid
prpPs :: [Double] -> [Double]
prpPs =  map pRngChk . inptMask . map pRngChk


-------- ======== Auger probabilities ======== --------

-- mathcal P (a,n,m) [for bA = True]
-- probabilities for a Auger electrons due to n 2a_1 holes,
-- when already m 1b_2, 3a_1, and 1b_1 electrons have been removed
--- T. Spranger and T. Kirchner, J. Phys B 37, 4159 (2004).
fAugerProbTobi :: (Integral b, Fractional a) => b -> b -> b -> a
fAugerProbTobi a n m = pA a n m'
    where m'        = fromIntegral m
          pA 0 0 _ =1
          pA 0 1 m =m/6
          pA 1 1 m =1-m/6
          pA 0 2 m =1/36*m^2
          pA 1 2 m =1/18*(6*m-m^2)
          pA 2 2 m =1/36*(6-m)^2
          pA _ _ _ =0


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
            bldr (k,l,m) = take k zs ++ take l os ++ take m ts
                where zs = repeat 0 -- zeros
                      os = repeat 1 -- ones
                      ts = repeat 2 -- twos


-- creates and filters a lit of all possible permutations of indices that lead to
-- k electrons at the projectile and l electrons in the continuum
-- the elements of that list have a structure of (a,ks,ls)
-- where 'a' is the number of Auger electron, ks the list of indices for capture
-- and 'l' the list for e in the continuum; these lists have n=9 elements
triInds :: Int -> Int -> [TriInds]
triInds k l  = (map toTrI . nub . filters . map strtr . concatMap permutations . prmSeeds n) (k+l)
    where   n = 9
            filters = filter f3 . filter f2 . filter f1
            strtr (a:kls) = (a,take 4 kls,drop 4 kls)
            f1 (_,ks,_) = sum ks == k  -- fits count of captured electrons
            f2 (_,ks,ls) = all (<=2) $ zipWith (+) ks ls  -- remove only 2 or less electrons per orbital
            f3 (a,ks,ls) = head ks + head ls >= a  -- Auger electrons must not exceed number of holes
            toTrI (a,ks,ls) = TrI a ks ls



-------- ======== prepare associated list of indices ======== --------
-- this ought to speed up computation, as the index list has to
-- be built only once, consecutive lines get it from memory

-- returns associated list with octal key
klTriInds :: [(Int, [TriInds])]
klTriInds = map (\x -> (f x, uncurry triInds x ) )   klTuples
    where f (a,b) = 10*a+b
          g (a,ks,ls) = TrI a ks ls

-- makes a list of all possible (k,l) tuples
klTuples :: [(Int, Int)]
klTuples =  concatMap (\x -> map (\y -> (x,y)) [0..8-x]) [0..8]

-- lookup, error when nothing is returned
myLookup :: Eq a => a -> [(a, b)] -> b
myLookup  x xs = case lookup x xs of
    Just y  -> y
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
trinProd ks ls os cs = tc * f cs ks * f is ls * f os ms
    where is = map (1-) $ zipWith (+) os cs  -- ionisation probability
          ms = map (2-) $ zipWith (+) ks ls
          tc = product $ zipWith g ks ls
          g _ 1 = 2 --trinomial coefficients
          g 1 _ = 2  -- when neither k,l,m equal 2
          g _ _ = 1  -- otherwise
          f rs is  -- powers (0,1,2) of probs, checks input
            | minimum is < 0   = error "trinProd, negative index"
            | maximum is > 2   = error "trinProd, forbidden index"
            | excdUt && minimum rs < 0   = error "trinProd, negative prob."
            | minimum rs *(-1) >= corrTld = error "trinProd, exceedingly negative prob."
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
trinSum :: Int -> Int -> [Double] -> Double
--trinSum :: [Double] -> Int -> Int -> Double
trinSum k l ocs =  sum $ map f is
    where is =  myLookup (10*k+l) klTriInds
          os = take 4 ocs
          cs = drop 4 ocs
          f (TrI a ks ls) =(*) (gA a n m) $ trinProd ks ls os cs
            where n = head ks + head ls
                  m = sum $ tail ks ++ tail ls
          gA a n m -- Auger correction factor
             | not bA && a == 0 = 1
             | not bA && a /= 0 = 0
             | otherwise        = fAugerProbTobi a n m


-- gets trinomial probabilities for all k, l combinations
klTrinProbs :: [Double] -> [Double]
klTrinProbs ps = map (\x -> uncurry trinSum x ps) klTuples

-- net removal probabilty calculated from sum q*p_q
netRecPQ  = \x -> x++f x
    where f    =  (:[]) . sum . zipWith (*) [1..]

-- calculate net ionization from k,l probabilities
-- and append it to the list with k,l probs
apdNetRecKL :: [Double] -> [Double]
apdNetRecKL pkls = (pkls++) . (:[]) . sum $ zipWith (*) (map q klTuples) pkls
    where q(k,l) = fromIntegral $ k+l


-------- ======== Formating Output ======== --------
--niceOut :: (Num a, Floating b) => [a] -> b -> [b] -> String
--niceOut :: (PrintfArg a, PrintfArg b) => [a] -> [b] -> String
niceOut [keV,b] = unwords .  (se:) . (sb:) .  map (printf "%14.7e") 
    where se  = printf  "%6.1f" keV; sb = printf  "%6.2f" b


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
            let pkl = klTrinProbs ps
            putStrLn $ niceOut eb $ apdNetRecKL pkl
--            putStrLn $ niceOut eb pnet pkl
            main


{- Input files ought to be organised like this example (strip the comments though):
--1      2       3                   4                   5                   6                   7                   8                   9                   10                  11                  12                  13                 14                 15                 16
--ELAB   B        | inito 1 (2a1)                                        |    |  inito 2 (1b2)                                       |    |  inito 3 (3a1)                                       |    |  i. 4 (1b1)  |    |  capture into projectile, originating from orbital:                   |
--[keV]  [au]       2a1->2a1            2a1->1b2            2a1->3a1            1b2->2a1            1b2->1b2            1b2->3a1            3a1->2a1            3a1->1b2            3a1->3a1            1b1->1b1            2a1                1b2                3b1                1b1
0020    0.40   0.1269480396e+00    0.5143592415e-01    0.2023602794e+00    0.9899421803e-01    0.5018767105e-02    0.4689024479e-01    0.3519103128e+00    0.6847626998e-01    0.2598374102e-01    0.7829925189e+00    0.3450854978e+00   0.5630833588e+00   0.3127868638e+00   0.1439064381e+00
0020    0.60   0.1648643287e+00    0.5618833391e-01    0.1706703322e+00    0.1480957416e+00    0.2789072993e-01    0.6727443586e-01    0.3493231849e+00    0.6161736900e-01    0.6928626147e-01    0.8186759488e+00    0.3306837587e+00   0.5359625416e+00   0.2898823084e+00   0.1354759864e+00

where the first two rows contain energy and collision
parameter, the consecutive columns probabilities
to occupy a specific state (0<p<=1), sorted in set of
3, 3, 3 and 1.  Followed by capture probabilities, for
4 orbitals.

Comments ought to be stripped off,
eg with grep, before piping to stdin.
-}

-- vim: set ts=4 sw=4 sts=4 nowrap et ai:
