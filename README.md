# mironovs_automodel
This is a 1d simulation source to model Cauchi problem for reduced non-linear light propagation equation in automodel case (thus 1d). The equation is modified KdV with variable koefficients.
Repository includes article with plot sources and scripts to postprocess results.
eq/build - builds test.cpp
eq/run - runs simulation (examples are provided)
plot/plot2d.asy - postprocesses dat files written bu simulation binary. Postprocessing includes center of mass, packet width and soliton extraction analysis. Semi-analytic data for solitions is plotted as well. Ot also dumps signals to singnal.dat file
plot/breeser.asy - processes signal.dat file in order to detect breezer parameters and dumps them to another breezer_params.dat file
plot/plot_breezer_params.asy - plots detected breeser parameters and semianalytinc data as well.
