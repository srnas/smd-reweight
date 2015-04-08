Steered-MD reweight
-------------------

smd-reweight is a tool to perform reweighting of steered MD simulations
using the formalism introduced in
[Colizzi and Bussi, JACS (2012)](http://dx.doi.org/10.1021/ja210531q).
Its usage is described in detail in
[Di Palma, Colizzi, and Bussi, Methods Enzymol. (2015)](http://dx.doi.org/10.1016/bs.mie.2014.10.055).
An earlier version of this software was also used to analyze simulations in
[Di Palma, Colizzi, and Bussi, RNA (2013)](http://dx.doi.org/10.1261/rna.040493.113).

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


