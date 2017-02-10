# elfin
A computational protein design tool based on repeat module construction. This is my MEng final year project and is under progress.

The main idea is to use repeatable protein parts like rigid construction blocks and build a shape towards the user's defined shape, very much like Lego. However, there is a non-trivial interaction relationship between the repeat modules: not only are some not compatible with others, they also combine in different but rigid angles.

### Status

10 Feb 2017: Stage 2/4 ongoing:
   * **Main Algorithm Formulation**
      * In progress, thinking of a Genetic Algorithm that operates on domain constrained genes.
         * Genes will be sequences of repeat modules restricted by interaction relationship.
         * Main kernel looking to include Kabsch (shape sampling still not solved...) and nearest-point poly-line scoring.
      * Fabio proposed Dead-end Elimination, still need to look into the literature.
         * TODO: Decide whether to include/replace with this method 

8 Feb 2017: Stage 1/4 complete:
   * **Pair and single module library**
      * Cleaned and generated, available as ./lib/python/GenXDB.py
      * Produced a file called xDB.json in ./res/ containing the necessary matrix transformation data.
   * **Benchmark protein (PDB) generator**
      * Done, available as ./lib/python/GenBench.py
      * Generated 10 fully random, length-10 benchmarks. They are found in in ./bm/
   * **Benchmark protein validation**
      * Done using Rosetta minimisation, results in ./bm/l10/scores/
      * Simple output extraction script available as ./lib/python/RMSDStat.py
      * Be ware of lever-effect - in our case shapes are mostly elongated rather than globular, hence the RMSD of around 5 should be acceptable.
      * Out of 10 random benchmark shapes only 1 slightly exceeded RMSD of 5, which is supposedly the international competition standard.
   * **Shape specification method**
      * Done, using a Matlab script available as ./lib/matlab/drawPoints.m
      * It spawns a figure that lets you draw a shape, and auto-scales to the appropriate size (protein-like scale using average module center distances).
      * Copy the point data into a .csv file without the square brackets.
         * TODO: auto export... maybe?
      * Some "fun" examples e.g. alphabets and "simple" shapes are found in ./bm/fun/*.csv
   * **Naïve greedy algorithm**
      * Done, available as ./lib/python/Greedy.py
      * Test scripts available as ./tests/Test*.py
      * Passed positive control test
      * Fails badly in Lv2 (oversampled) and Lv3 (under-sampled) inputs due to spiral dead-end modules, 
         * This is expected because once the naïve algorithm makes a mistake due to its "myopia", it can run into a "spiral" type module that causes collision and be forced to terminate prematurely before meeting length criteria.
      * Can take either JSON (e.g. those in ./bm/l10/) or CSV (e.g. those in ./bm/fun/) as input target shape
   * **PyMol integration***
      * Benchmark generation supports drag-and-drop, found in ./lib/python/GenBenchmark.py
      * Axes and CSV drawing functions available in ./pymol/elfin.py
         * TODO: replace line with cylinders to make width actually work


### Setup 
1. Install all of [must-have tools]() and optionally [complementaries]()

### Must-have Tools: 
1. [BioPython](http://biopython.org/wiki/Download) library for PDB IO
2. [VirtualEnv](https://virtualenv.pypa.io/en/stable/) for separation of python environment

### Optional Complementaries:
1. [PyMol]() for 3D visualisation of PDB and CSV
2. [Rosetta](https://www.rosettacommons.org/software/license-and-download) for minimisation and confirmation of results
