Compiling
---------
Make sure you have g++ and libboost-thread-dev, then:
cd ann_1.1.2
make linux-g++
cd ../
make

Running
-------
  avd [-d dimensions] [-n points] [-r representatives] [-t threads] [-s side multiplier]      

  -d: number of dimensions, at least 3
  Specify the number of dimensions for the input points 

  -n: number of points, at least 2
  Specify the number of input points

  -r: number of representatives per leaf, at least 1
  Specify the number of representatives per leaf. A bigger number produces a smaller tree, but the queries are a bit slower.
  Don't use 1, as the tree will be quite large.

  -t: number of threads, at least 1
  Specify the number of threads to use for constructing the AVD.

  -s: side_multiplier, 1 for vanilla squares
  Specify 1 to run Method B. Specify a bigger number to run Method B'

  -f: fill leaves with reprs completely
  Specify this flag to fill leaves with representatives. Needed with Method B'.

  -q: just construct, no queries
  Just test the construction, don't run any queries

Examples
--------
Method B : ./avd -d 3 -n 128 -r 4 -s 1 -t 4 
Method B': ./avd -d 3 -n 128 -r 8 -s 4 -t 4 -f

Input sets
----------
Two input sets are included:
"random_numbers" contains 40.000 randomly generated numbers
"real_data" contains 10.000, 7-dimensional points representing coordinates and radiation emission of galaxy objects

Default choices and how to change them
--------------------------------------
Default input set: "real_data"
To change to "random_numbers", open scale.cpp, uncomment line 38, comment lines 39 and 48 to 50

Default number of query points: 100000
To change the number of query points, open query.cpp and line 135 and change the exponent of pow(10,4)

Default method of calculating the query points: random in each leaf (fixed number of queries for each leaf, random inside the area of the leaf)
To change the method to random in the hypercube, open query.cpp, comment line 138 and uncomment line 140

In all cases, don't forget to recompile!

Interpreting results
--------------------
Unfortunately the current format is not very user friendly. Here is a typical output and the breakdown of each column.

3 128 4 1 1626 1332 9318 54 53 14 
100000 50000 10000 65513 31387 4542 841 179 36 0 102194 304 0 0 0 0 0 

1st line:
Dimensions of input points     -> 3
Number of input points         -> 128
Max Representatives per leaf   -> 4
Side multiplier                -> 1
Number of well separated pairs -> 1626
Quadtree nodes                 -> 1332
Quadtree leaves                -> 9318
Total construction time        -> 54
Merging time                   -> 53
Max depth of quadtree          -> 14

2nd line
Total query time for KD exact        -> 100000
Total query time for KD approx       -> 50000
Total query time for AVD             -> 10000
Exact results of KD approx           -> 65513
Results KD approx with 0%<e<10%      -> 31387
Results KD approx with 10%<e<20%     -> 4542
Results KD approx with 20%<e<30%     -> 841
Results KD approx with 30%<e<40%     -> 179
Results KD approx with 40%<e<50%     -> 36
Results KD approx with e>60% (wrong) -> 0
Results AVD with 0%<e<10%            -> 102194
Results AVD with 10%<e<20%           -> 304
Results AVD with 20%<e<30%           -> 0
Results AVD with 30%<e<40%           -> 0
Results AVD with 40%<e<50%           -> 0
Results AVD with e>60% (wrong)       -> 0

Future Work
-----------
Provide the option to save the tree
Provide more runtime flags
