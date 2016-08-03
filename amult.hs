-- multinomial analysis with Auger
-- for H2O with states 2a1, 1b2, 3a1, 1b1


-- mathcal P (k,m,a)
pA :: (Integral a, Fractional b) => a -> a -> a -> b
pA k a n 
    | k == 0 && a == 0  =1
    | k == 0 && a == 1  =0
    | k == 1 && a == 0  =m
    | k == 1 && a == 1  =1-m/6
    | k == 2 && a == 0  =1/36*m^2
    | k == 2 && a == 1  =1/18*(6*m-m^2)
    | k == 2 && a == 2  =1/36*(6-m)^2
    where m = fromIntegral n

-- no Auger probs
pNoA :: (Integral a, Fractional b) => a -> a -> a -> b
pNoA k 0 m = 1
pNoA k 1 m = 0

-- ratio of removal prob and occupation
po :: (Fractional c) => c -> c
po = (+) (-1) . ((/) 1)

-- get the n elements starting at the m-th position of a list 
takeAt :: Int -> Int -> [a] -> [a]
takeAt n m = (take n .snd . splitAt (m-1))

-- organize data of one line into (e,b) and occupations
orgData :: [a] -> ((a, a), [a])
orgData ds = ((e,b),os)
    where   os  = (snd . splitAt 2) ds
            e   = head ds
            b   = (head . tail) ds

-- locate and sum occupation of each initial condition
sumOrbOccupation :: (Num b) => [b] -> [b]
sumOrbOccupation os  = map ($os) [o2a1,o1b2,o3a1,o1b1]
    where   o2a1 = (sum . takeAt 3  1)
            o1b2 = (sum . takeAt 3  4)
            o3a1 = (sum . takeAt 3  7)
            o1b1 = (sum . takeAt 1 10)

--main: (sumOrbOccupation . snd . orgData ) line

-- binomial coefficients
binom :: Int -> Int -> Int
binom n 0 = 1
binom 0 k = 0
binom n k = binom (n-1) (k-1) * n `div` k

-- probability not to remove any electrons
noRem :: Num c => [c] -> c
noRem = sum . map (^2) 

-- ns ::(Integral a) => [a,a,a,a]
--probPi ns 


-- repetitive elements of calculation
facPi :: (Fractional c) => Int -> c -> c
facPi 0 p = 1 
facPi 1 p = ((*) 2 . (^^n) . po ) p
facPi 2 p = ((*) 1 . (^^n) . po ) p
facPi n p = ((*) b . (^^n) . po ) p
    where b = (fromIntegral . (binom 2) ) n


-- example input line
ds = [20.0,0.4,0.1269480396,5.143592415e-2,0.2023602794,9.899421803e-2,5.018767105e-3,4.689024479e-2,0.3519103128,6.847626998e-2,2.598374102e-2,0.7829925189]

-- vim: set ts=4 sw=4 sts=4 et ai:
