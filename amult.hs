-- multinomial analysis with Auger
-- for H2O with states 2a1, 1b2, 3a1, 1b1
import Data.List

-------- ======== input related ======== --------

-- get the n elements starting at the m-th position of a list 
takeAt :: Int -> Int -> [a] -> [a]
takeAt n m = (take n .snd . splitAt (m-1))


-- check if p is in valid range 0<p<=1
pRngChk :: (Num a, Ord a) => a -> a
pRngChk p
    | p < 0    =error "negative probability"
    | p > 1    =error "probability exceeds unity"
    | p == 0   =error "zero probability not considered"
    | otherwise = p


-- locate and sum occupation of each initial condition
sumOrbOccupation :: (Num b) => [b] -> [b]
sumOrbOccupation os  = map ($os) [o2a1,o1b2,o3a1,o1b1]
    where   o2a1 = (sum . takeAt 3  1)
            o1b2 = (sum . takeAt 3  4)
            o3a1 = (sum . takeAt 3  7)
            o1b1 = (sum . takeAt 1 10)


-- organize data of one line into (e,b) and occupations
-- orgData :: [a] -> ((a, a), [a])
orgData :: (Num b, Ord b) => [b] -> ((b, b), [b])
orgData ds = ((e,b),os)
    where   os  = (map pRngChk . sumOrbOccupation . map pRngChk . snd . splitAt 2) ds
            e   = head ds
            b   = (head . tail) ds



-------- ======== Auger probabilities ======== --------
-- mathcal P (k,a,m)
pA :: (Integral a, Fractional b) => a -> a -> b -> b
pA 0 0 _ =1
pA 0 _ _ =0
pA 1 0 m =m
pA 1 1 m =1-m/6
pA 1 _ _ =0
pA 2 0 m =1/36*m^2
pA 2 1 m =1/18*(6*m-m^2)
pA 2 2 m =1/36*(6-m)^2
pA _ _ _ =0


-- no Auger probs
pNoA :: (Integral a, Fractional b) => a -> a -> b -> b
pNoA _ 0 _ = 1
pNoA _ _ _ = 0


-- ratio of removal prob and occupation
po :: (Fractional c) => c -> c
po = (+) (-1) . ((/) 1)

-- binomial coefficients
binom :: Int -> Int -> Int
binom n 0 = 1
binom 0 k = 0
binom n k = binom (n-1) (k-1) * n `div` k

-- probability not to remove any electrons
noRem :: Num c => [c] -> c
noRem = sum . map (^2) 

-- repetitive elements of calculation
facPi :: (Fractional c) => Int -> c -> c
facPi 0 = ((+1) . (*0)) -- neutral element
facPi 1 = ((*) 2 . po )
facPi 2 = ((^^2) . po )
facPi k= ((*) b . (^^k) . po )
    where b = (fromIntegral . (binom n) ) k
          n = 2

-------- ======== Multinomial Products ======== --------
-- is indices [n,m1b2,m3a1,m1b1,a]
-- ps probability array, corresponding to indices
-- q  final target charge state
-- n  electrons lost from 2a1
-- m1b2,m3a1,m1b1
--    electrons lost from respective orbitals
-- a  number of Auger electrons
-- bA boolean turning Auger analysis on and off
--expression bA q n ms a p
multProd :: (Eq a, Fractional a) => Bool -> Int -> [Int] -> [a] -> a
multProd bA q is ps
    | q == (sum is) && isA /= 0 = ((*isA) . prod) ps
    | not bA && a /= 0          =0
    | length is /= 1+length ps    = error "Unfortunately, array lenghts of indices and probs in function 'sumTerm' do not match."
    | otherwise                 =0
    where   a = last is
            n = head is
            m = (sum . tail) nms
            nms = init is  -- removes a from list
            isA = f  n a (fromIntegral m)  -- caclulate Auger probs
            f | bA          = pA
              | otherwise   = pNoA
            prod = product . (zipWith facPi is) 

-- permutate indices
prmSeeds :: (Eq a, Num t, Num a) => a -> [[t]]
prmSeeds q  | q==0  = [zs]
            | q==1  = map cstr [[1,1,0]]
            | q==2  = map cstr [[2,2,0],[1,0,1]]
            | q==3  = map cstr [[3,3,0],[2,1,1]]
            | q==4  = map cstr [[4,4,0],[3,2,1],[2,0,2]]
            | q==5  = map cstr [[5,5,0],[4,3,1],[3,1,2]]
            | q==6  = map cstr [        [5,4,1],[4,2,2],[3,0,3]]
            | q==7  = map cstr [                [5,5,1],[4,1,3]]
            | q==8  = map cstr [                        [5,2,3],[4,0,4]]
            | otherwise = error "q must not exceed 8"
            where zs = [0,0,0,0,0] -- empty list
                  os = [1,1,1,1,1]
                  ts = [2,2,2,2,2]
                  cstr [k,l,m] = (drop k) zs ++ (take l) os ++ (take m) ts

-- permutate indices step two
prmInds :: (Integral a) => a -> [[Int]]
prmInds  = concat . map (nub . filter f . permutations) . prmSeeds
    where f xs = head xs >= last xs  -- filter out permutations where a exceeds n


-- use all possible permutations of the indices [n,m...,a]
multSum q ps = ((*pnot) . sum . map (`f` ps) . prmInds) q
    where   f    = multProd bA q
            pnot = ((^^2) . product)   ps
            bA   = False


-- example input line
--ds = [20.0,0.4,0.1269480396,5.143592415e-2,0.2023602794,9.899421803e-2,5.018767105e-3,4.689024479e-2,0.3519103128,6.847626998e-2,2.598374102e-2,0.7829925189]
--ps =  (sumOrbOccupation . drop 2) ds

-- use `grep -Ev "^$|^#" < $1` to filter comments from input
main =  do
    line <- getLine
    if null line
        then return ()
        else do
            let rawData = (map read .words) line :: [Double]
            let netOcpData = orgData rawData
            let ps =snd netOcpData
            (print . fst ) netOcpData
            (print . map (`multSum`  ps)) [1..8]
            main


{- Input files ought to be organised like this example:
#1      2       3                   4                   5                   6                   7                   8                   9                   10                  11                  12
#               >> inito 1 (2a1)                                      <<    >> inito 2 (1b2)                                      <<    >> inito 3 (3a1)                                      <<    >> i. 4 (1b1) <<
#ELAB   B       2a1->2a1            2a1->1b2            2a1->3a1            1b2->2a1            1b2->1b2            1b2->3a1            3a1->2a1            3a1->1b2            3a1->3a1            1b1->1b1
 0020    0.40   0.1269480396e+00    0.5143592415e-01    0.2023602794e+00    0.9899421803e-01    0.5018767105e-02    0.4689024479e-01    0.3519103128e+00    0.6847626998e-01    0.2598374102e-01    0.7829925189e+00
 0020    0.60   0.1648643287e+00    0.5618833391e-01    0.1706703322e+00    0.1480957416e+00    0.2789072993e-01    0.6727443586e-01    0.3493231849e+00    0.6161736900e-01    0.6928626147e-01    0.8186759488e+00

where the first two rows contain energy and collision
parameter, the consecutive columns probabilities 
to occupy a specific state (0<p<=1), sorted in set of
3, 3, 3 and 1.  Comments ought to be stripped off,
eg with grep, before piping to stdin. 
-}

-- vim: set ts=4 sw=4 sts=4 nowrap et ai:
