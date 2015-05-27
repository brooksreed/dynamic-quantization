# dynamic-quantization

Implements and compares a few dynamic quantization algorithms from the literature.  

The 08 TAC LP-based method requires an installation of AMPL/CPLEX.  

'azuma08fwd';           % Azuma 08 Automatica

'azuma08feedback';      % uses Azuma 08 Automatica in feedback

'azumaLPfwd';           % uses Azuma 08 TAC finite-horizon LP

'azumaLPfeedback';      % uses Azuma 08 TAC in feedback

'minami07';             % uses Minami 07 CDC (specific to feedback)

The quantizerCompare script runs a test using one of the above methods and one of a couple simple test systems.  

The azumaTest function implements the quantizer and runs the simulation.  

The designAzuma08TAC function designs the LP-based quantizer descriped in the 08 TAC paper.  

The staticNearestNeighbor function is a simple nearest-neighbor quantizer.  

The writeAMPLRunOptions function is used to setup the AMPL/CPLEX solver call.  



