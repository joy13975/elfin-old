# Elfin
A computational protein design tool based on repeat unit construction (RUC). This is my MEng final year project and is still being updated.

The idea of Elfin is to use repeatable protein parts as if they are rigid construction blocks (think Lego) and build a shape towards the user's desire. A poster I presented for my final year project can be found [here](res/pdf/poster.pdf).

![alt tag](res/png/ProteinBristol.png)
Figure: actual output of the Elfin GA with "Bristol" written by hand as input. The visuals were created using PyMol.

Jump to: 
[Python Setup](#python-setup)

[Compile GA](#compile-ga)

[Usage](#usage)

### Project Status

31 Mar 2017: Stage 4/4:
   * **GA complete and shown to work on custom inputs**
      * Passes all benchmarks
      * Two grid search passes done - find results in GSOutV1 and GSOutV2
      * Implemented diversity enforcenement using CRC
      * PyMol script ./src/Python/Compare.py available for showing user input shape versus Elfin output
      * Found OMP_SCHEDULE=dynamic,16 to be best for most machines
      * Update Makefile to account for clang++
      * Made Setup, Compile and Usage sections more complete in README.md
   * **TODO: OpenMP 4.5 - might need to turn classes into structs to properly map them...**
   * **Thesis and journal paper (Structural Biology) are currently being prepared**

16 Mar 2017: Stage 3/4:
   * **GA mostly implemented (~90% as far as this project is concerned)**
      * Includes unit tests of critical code (ideally should be every bit of the code, but I have time constraints)
         * Geomtry & Maths
         * Randomiser & Parallel randomiser
         * Evolution functions
      * Features multiple evolution operators, namely:
         * Crossing (combine 2 highly ranked parents at a random legal point)
         * Point mutation (firstly inherit from one parent, then apply either of {insert, swap, delete})
         * Limb mutation (sever from a random point in sequence, then grow a new "limb")
         * Fully randomise (self explanatory)
      * Implemented straightforward parallelism
         * Simple #pragma omp parallel for
         * Using parallel version of stdlib functions (-D_GLIBCXX_PARALLEL) so I don't need to re-implement stuff like sorting.
         * TODO: further investigate memory access patterns and optimisation opportunities
         * TODO: test on GPUs using OpenMP 4.5
         * TODO: benchmark on a couple of different platforms.
      * TODO: grid search to find good GA parameters - DONE
      * TODO: a more complete documentation - N/A
      * TODO: write a more visually helpful version of Synth.py to use Kabsch and align two sequences so they can appear to be more aligned in PyMol (instead of fixing one end to origin) - DONE

10 Feb 2017: Stage 2/4:
   * **Main Algorithm Formulation**
      * In progress, thinking of a Genetic Algorithm that operates on domain constrained genes.
         * Genes will be sequences of repeat modules restricted by interaction relationship.
         * Main kernel looking to include Kabsch (shape sampling still not solved...) and nearest-point poly-line scoring.
      * Fabio proposed Dead-end Elimination, still need to look into the literature.
         * Consider as alternative but no time to implement every possible algorithm 

8 Feb 2017: Stage 1/4 complete:
   * **Pair and single module library**
      * Cleaned and generated, available as ./src/Python/GenXDB.py
      * Produced a file called xDB.json in ./res/ containing the necessary matrix transformation data.
   * **Benchmark protein (PDB) generator**
      * Done, available as ./src/Python/GenBench.py
      * Generated 10 fully random, length-20 benchmarks. They are found in in ./bm/
   * **Benchmark protein validation**
      * Done using Rosetta minimisation, results in ./bm/l20/scores/
      * Simple output extraction script available as ./src/Python/RMSDStat.py
      * Be ware of lever-effect - in our case shapes are mostly elongated rather than globular, hence the RMSD of around 5 should be acceptable.
      * Out of 10 random benchmark shapes only 1 slightly exceeded RMSD of 5, which is supposedly the international competition standard.
   * **Shape specification method**
      * Done, using a Matlab script available as ./lib/matlab/drawPoints.m
      * It spawns a figure that lets you draw a shape, and auto-scales to the appropriate size (protein-like scale using average module center distances).
      * Copy the point data into a .csv file without the square brackets.
         * TODO: auto export... - SCRAPED
      * Some "fun" examples e.g. alphabets and "simple" shapes are found in ./bm/fun/*.csv
   * **Naive greedy algorithm**
      * Done, available as ./src/Python/Greedy.py
      * Test scripts available as ./tests/Test*.py
      * Passed positive control test
      * Fails badly in Lv2 (oversampled) and Lv3 (under-sampled) inputs due to spiral dead-end modules, 
         * This is expected because once the naïve algorithm makes a mistake due to its "myopia", it can run into a "spiral" type module that causes collision and be forced to terminate prematurely before meeting length criteria.
      * Can take either JSON (e.g. those in ./bm/l20/) or CSV (e.g. those in ./bm/fun/) as input target shape
   * **PyMol integration***
      * Benchmark generation supports drag-and-drop, found in ./src/Python/GenBenchmark.py
      * Axes and CSV drawing functions available in ./src/PyMolUtils/LineUtils.py
         * TODO: replace line with cylinders to make width actually work - DONE


### [Python Setup](#python-setup) 
This is for setting up the Python virtual environment to use scripts in ./src/Python.

1. Install all of [must-have tools](#must-have-tools) and optionally, [complementaries](#optional-tools)
2. To set up the Python environment:
```
cd elfin                         #if not already in the directory
virtualenv venv                  #the name 'venv' is required
. ./activate                     #get into the virtual environment
pip install -r requirements.txt  #automatically fetch and install the libraries needed (locally)
```

### [Compile GA](#compile-ga)
This is for compiling the C++ source code of the Elfin GA.

1. Install submodule jutil - a generic utility library I wrote in another repository.
```
git submodule init 
git submodule update
```

2. Compile the GA
```
cd src/GA
make
```

Notes:
-You can specify your compiler of choice by doing e.g. for clang++: ```make CXX=clang++```.
-For clang, you will need to specify the C++ standard library include path, which depends on where you installed GCC. See .src/GA/Makefile for details (the INCS variable).
-For clang, you also need libiomp include and library paths specified (again see Makefile).
-To speed up the compilation, use the -j flag with the number of cores on your system.
-To build without OpenMP, you can do ```make OMP=no```

### [Usage](#usage)
Once you have compiled the GA successfully, you can test run it with:
```
./bin/elfin
```

To get help:
```
./bin/elfin -h
```

Normally you would do something like this:
```
./bin/elfin -gps <POPULATION_SIZE> -i <INPUT_FILE>
```

Notes:
-Default configuration is in src/GA/config.json. 
-Command-line arguments override whichever configuration file is used.

### [Must-have Tools](#must-have-tools): 
1. Python 2.7.9+
2. [VirtualEnv](https://virtualenv.pypa.io/en/stable/) for separation of python environment

### [Optional Tools](#optional-tools):
1. [PyMol]() for 3D visualisation of PDB and CSV
2. [Rosetta](https://www.rosettacommons.org/software/license-and-download) for minimisation and confirmation of results
