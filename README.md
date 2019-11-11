# STCov

STCov is a program for calculating statistics (mean, variance, covariance, shared branch length) of pairwise coalescence times given a rooted, bifurcating species tree using the multi-species coalescent.

## Installation Instructions

STCov is implemented in C++ (a precompiled binary is not yet added to this github, but will be added soon). 

If you choose to build the software from scratch, follow the build instructions below

## Build Instructions

STCov requires the following libraries and executables in order to compile and run: 

* A C++-11 compiler (gcc 4.8 or later, for example)
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) - Used for linear algebra and matrix calculations

On OSX, from inside the STCov folder, build STCov by running: 

```
clang++ -w -std=c++11 -stdlib=libc++ -I /usr/local/include/ main.cpp -o STCov
```
Failure to compile usually has to do with the location of the Eigen program, ensure that it is able to be located in the pathway ```/usr/local/include/```.

## Usage

STCov takes only a configuration file "ExampleConfig.txt" and can be executed using: 
```
./STCov ExampleConfig.txt
```
where ExampleConfig.txt is in the same directory as the executable. 

### File format and required arguments

The configuration file describes the species tree to be read in to STCov, there are three types of arguments that must appear, in order. 

* The size of the tree, indicated by the integer number of present day species, denote this by ```N```.
```
nspecies = N 
```
* The population sizes of all present day species (leaf node population sizes), in ```N``` consecutive lines. (Note that present day species must labeled by numbers ```1```...```N```)
```
-n 1 popsize_1
-n 2 popsize_2
...
-n N popsize_N
```
* In order of occurence (most recent to most ancient), the split events on the tree parameterized by a split time, species involved (total of 2) and the population size for the new merged branch. There must be a total of ```N-1``` of these events.
```
-t split_time i j popsize_i,j
```
Where ```i``` and ```j``` are two of the modern day species, where all lineages of ```i``` are merged into the lineages of ```j``` at time ```split_time``` and the population size of this new ancient branch becomes ```popsize_i,j```. In later split events, population ```j``` represents the subtree of all branches that have merged into ```j``` (meaning we should NOT use the symbol ```i``` in any later merge events). 

### A note on units and notation

We use tree parameterization in a fashion similar to (and inspired by) the simulation program ```ms``` (Hudson 2012). However our units on time and population sizes are not identical to this program. STCov measures time in units of 2N generations, whereas ```ms``` measures time in 4N generations. As well, STCov uses haploid population sizes, while ```ms``` uses diploid population sizes. So to convert between the two, divide times in STCov by 2 and mulitply population sizes by 2 when  porting them to ```ms```.

### Example configuration file

Here we demonstrate an example configuration file for a tree of 4 species. 

We input the topology (4,(3,(2,1)) with the following parameters: 
* Leaf population sizes (4 in total): N_1 = 1.0, N_2 = 0.75, N_3 = 2.4, N_4 = 0.25
* Split times (4-1 = 3 in total): T_1,2 = 0.4, T_(1,2),3 = 0.8, T_(1,2,3),4 = 1.2
* Internal population sizes (4-1 =3 in total): N_1,2 = 2.1, N_(1,2),3 = 0.15, N_(1,2,3),4 = 1.8

The configuration file then needs to look as follows: 
```
nspecies = 4
-n 1 1.0
-n 2 0.75
-n 3 2.4
-n 4 0.25
-t 0.4 1 2 2.1
-t 0.8 2 3 0.15
-t 1.2 3 4 1.8
```
(Note that any line starting with a ```#``` will be ignored and can be used for commenting)

## Output 

STCov prints on the command line the tree read in from the configuration file. Verify that this tree accurately represents the one provided in the config file. 

For a provided number of species, ```N```, STCov assumes a total of 4 individuals from each species have been sampled. These individuals are labeled using ```a```,```b```,```c```,```d``` as subscripts of their species of origin. For example, the first individual from species 2 has the label ```2_a```. We assume all individuals are interchangeable. 

As STCov calculates statistics of pairwise coalescence times, each pair is indexed using a comma separated list of their names, e.g. ```1_b, 5_c```. 

There are 4 files output by STCov

* ```Labels.txt``` used to read the next 3 files, is a list of all pairs of individuals and their row/column index in the data files. There are (```N choose 2```) possible pairings of indiviudals, which is the number of rows in this file.  
* ```Means.txt``` outputs the expected time to coalescence for a each pair of individuals conditional on the species tree from the configuration file. 
* ```Covariance.txt``` outputs the variance/covariance matrix between time to coalesce for two pairs of individuals. The diagonal corresponds to the variance in coalescence time for a given pair, and the off diagonal corresponds to covariances. Use the ```Labels.txt``` file as the index for both the rows and columns of this file. 
* ```SharedBranchLength.txt``` outputs the expected amount of shared branch length between two pairs of indiviudals. The diagonal corresponds to 2 times the expected time to coalescence for an individual pair. Use the ```Labels.txt``` file as the index for both the rows and columns of this file. 

### Interpreting the output file, and example. 

Suppose we are interested in the covariance between the coalescence time T_```1_a, 5_c``` and T_```1_b,4_d```. To find this value in ```Covariance.txt```, we look up the index values for ```1_a, 5_c``` and for ```1_b, 4_d``` from ```Labels.txt``` and then find the corresponding entry in ```Covariance.txt```. Note the the covariance and shared branch length matrices are symmetric, so the order of (row,column) does not matter. 





