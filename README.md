## About Haploflow
Haploflow is a strain-aware viral genome assembler from short read sequence data. 
It uses a flow algorithm on a deBruijn graph data structure to resolve viral strains. Haploflow is still actively under development and release 0.1 is freely available under the 
_license_. It was written entirely in C++ and currently only works on UNIX systems.\
This README lists the requirements, installation information and a short tutorial on how to use Haploflow and its parameters.

## Requirements
 - CMake >= 2.8
 - Boost >= 1.55
 
## Installation
First, clone this repository using `git clone` _address_, then enter the directory which you cloned Haploflow to and create a build folder,
e.g. `mkdir build`. Enter this new directory and run cmake with `cd build; cmake ..`. This will create a Makefile which you can then run
to create the Haploflow executable: `make`
This should create a `haploflow` executable in your build directory.

## Short tutorial
Using the executable you can show the help and parameters using `./haploflow --help`. This lists the parameters as follows:
~~~~
HaploFlow parameters:
  --help                                Produce this help message
  -  [ --read-file ] arg                read file (fastq)
  -  [ --dump-file ] arg                deBruijn graph dump file produced by 
                                        HaploFlow
  --log arg                             log file (default: standard out)
  -  [ --k ] arg (=41)                  k-mer size, default 41, please use an 
                                        odd number
  -  [ --out ] arg                      folder for output, will be created if 
                                        not present. WARNING: Old results will 
                                        get overwritten
  -  [ --error-rate ] arg (=0.0199999996)
                                        percentage filter for erroneous kmers -
                                        kmers appearing less than relatively e%
                                        will be ignored
  --create-dump arg                     create dump of the deBruijn graph. 
                                        WARNING: This file may be huge
  --from-dump arg                       run from a Haploflow dump of the 
                                        deBruijn graph.
  -  [ --two-strain ] arg (=0)          mode for known two-strain mixtures
  -  [ --strict ] arg (=1)              more strict error correction, should be
                                        set to 5 in first run on new data set 
                                        to reduce run time. Set to 0 if low 
                                        abundant strains are expected to be 
                                        present
  -  [ --filter ] arg (=500)            filter contigs shorter than value
  -  [ --thresh ] arg (=-1)             Provide a custom threshold for 
                                        complex/bad data 
~~~~

The input reads are given with the `--read-file` option and the output directory with `--out`, which are the only required options. 
Haploflow will then run with default parameters.\
The most important other parameters are `k`, the *k*-mer size of the deBruijn graph. This is 41 by default, increasing this value might
improve assembly for large read lengths or very deep sequencing runs.
`error-rate` is the next parameter, which determines a lower bound of coverage or detection limit of different strains and 
is a percentage value. By default this value is set to `0.02`, because Illumina data is expected to have less than 2% errors. 
Setting this value too low can cause Haploflow to run far slower; setting it too high will prevent Haploflow from finding lower abundant
strains.\
The `strict` parameter is complementary in the sense that it determines an overall lower bound for read coverage. Setting it to `-1` 
imposes no constraints, `0` will use the inflection point of the coverage histogram and every value `≥1` will result in use of a sliding window over the coverage histogram to determine the lower bound.\
Finally, the last error correction parameter is `thresh`: it is mutually exclusive with the `strict` parameter and will overwrite its
value if set. This parameter sets a fixed threshold below which *k*-mers are ignored.
Finally, Haploflow by default filters contigs shorter than 500 bp. This value can be changed using the `filter` option. 
The parameters `create-dump`, `from-dump` and `dump-file` are just needed if the deBruijn graph is supposed to be written to a file to be
reused in another run. This file is possibly huge (because uncompressed), so use with caution.
