Steered-MD reweight
-------------------

smd-reweight is a tool to perform reweighting of steered MD simulations
using the formalism introduced in Colizzi and Bussi, J. Am. Chem. Soc., 134, 5173 (2012).

First compile the source smd-reweight.cpp with a suitable compiler, e.g.

    ./configure CXX=g++
    make

or

    ./configure CXX=icpc
    make

Then ask for help

    ./smd-reweight --help

The program just reads a text file with columns and a simple format.
Using straightforward file manipulation it can be made compatible with
any molecular dynamics engine.


