# SWORD

SWORD (Smith Waterman On Reduced Database) is a fast and sensitive software for protein sequence alignment. It is implemented as a module for SW# library which is based on CUDA enabled GPUs. SWORD consists of two steps, first being a heuristic and the second being optimal alignment phase. In the first step, for each query, it reduces the database to A sequences that score best with the query, which are in the second step aligned to the query using SW# library. SWORD utilizes multithreading.

## DEPENDENCIES

### LINUX and MAC OS

Application uses following software:

1. SW# library - freely available from http://sourceforge.net/projects/swsharp/
2. gcc 4.*+
3. nvcc 2.*+

## INSTALLATION

### LINUX

Makefile is provided in the project root folder. SW# must be downloaded and compiled before compiling SWORD. To generate an executable,
first move the SWORD folder to SW# root directory. Inside SWORD root, run:
    
    make

After running make, an executable named sword will appear in the current directory.

The software should also run on MAC OS, but it was not tested in that environment.

## EXAMPLES

All examples assume that make has been run and that SWORD was successfully compiled.
Simplest protein fasta database search can be executed using the following command:

    ./sword -i <query> -j <database>

This will run a search using the default, sensitive mode.

\*note: By default, SWORD will try to cache the database used if no cache file exists to speed up future runs on the same database. This behavior is not mandatory and can be disabled with the --nocache flag.


For the complete list of parameters and their descriptions run the following command:

    ./sword -h (or ./sword --help)

To remove SWORD executable, run:

    make clean
    
