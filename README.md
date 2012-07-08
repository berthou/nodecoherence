Simulating Grover's algorithm with BECs (Bose-Einstein condensate) with no decoherence.
=============

Download
--------

Run `git clone https://github.com/berthou/nodecoherence.git` or download & extract `https://github.com/berthou/nodecoherence/zipball/master`

Go to the main folder :

`cd nodecoherence/`


Before compiling and running the programm you need to install dependencies.

Dependencies
------------

The following dependencies are required :

* [GNU Scientific Library](http://www.gnu.org/software/gsl/) -- run `make dependencies` to install the library (linux only) or install it by yourself
* [gnuplot](www.gnuplot.info/) -- (Optional)


Compilation & execution
-----------------------

To make sure that everything went fine run `make -s testing`. The output should be :

`Np = 1
Done.
`
If you have gnuplot installed you can view the output graph by runnning `./run-for-each-value.sh 1 1 1`.

You can change the arguments but the 2nd must be greater than the 1st one and the last one should'nt be greater than 4 (too much computation needed).
