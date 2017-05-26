# SWORD

SWORD (Smith Waterman On Reduced Database) is a fast and sensitive software for protein sequence alignment. SWORD consists of two steps, first being a heuristic and the second being optimal alignment phase. In the first step, for each query, it reduces the database to A sequences that score best with the query, which are in the second step aligned to the query using OPAL library (https://github.com/Martinsos/opal). SWORD utilizes multithreading.

## DEPENDENCIES

### LINUX and MAC OS

Application uses following software:

1. SSE4.1 or higher
2. gcc 4.8+

## INSTALLATION

### LINUX

To build SWORD run the following commands from your terminal:

    git clone --recursive https://github.com/rvaser/sword.git sword
    cd sword/
    make

After running make, an executable named sword will appear in the current directory.

#### Troubleshooting

If you have cloned the repository without '--recursive', run the following commands before running 'make':

    git submodule init
    git submodule update

The software should also run on MAC OS, but it was not tested in that environment.

## EXAMPLES

All examples assume that make has been run and that SWORD was successfully compiled.
Simplest protein FASTA database search can be executed using the following command:

    ./sword -i <query> -j <database>

This will run a search using the default, sensitive mode.

For the complete list of parameters and their descriptions run the following command:

    ./sword -h (or ./sword --help)

To remove SWORD executable, run:

    make clean

## Contact information

For additional information, help and bug reports please send an email to: robert.vaser@fer.hr.

## Acknowledgement

This work has been supported in part by Croatian Science Foundation under the project UIP-11-2013-7353 and in part
by the Foundation of the Croatian Academy of Sciences and Arts.
