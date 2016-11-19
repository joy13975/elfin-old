pdbs contains two chains named A and B that are two modules from which extract the relative orientation, by calculating teh center of mass using the CA atoms. It is not necessary to use all atoms.

each file is named module1(chain A)-module2(chainB).pdb The dash separates the name of two modules, the _ is part of the name of a single module.
some files contain the interaction between a module and another copy of itself, e.g. D14-D14.pdb

when building a whole protein chain the names carry the information of what modules are compatible.
For example, to build a three module pieceformed by D79, D79_j1_D54 and D54 we will use the two pdbs D79-D79_j1_D54.pdb and D79_j1_D54-D54.pdb
