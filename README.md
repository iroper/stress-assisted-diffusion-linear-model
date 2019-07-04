# stress-assisted-diffusion-linear-model
MATLAB code to solve the model and produce the figures in the paper 'The effect of mechanical stress on lithium distribution and geometry optimisation for multi-material lithium-ion anodes'.

To produce the plots as seen in the paper, add '%' to the start of the line under each heading (lines 31, 45, 131, 195, 302, 416) in 'Produce_plots.m' to uncomment different sections of the code (so not everything is calculated each time). Then run this code.

Parameters are set in Produce_plots.m but can easily be changed by the user to further explore the model presented in this paper. If other materials are to be used, OCV data must be produced in a .mat file and the mechanical and chemical constants must be updated in constants.m.

The extra .m files are functions used by Produce_plots.m. The .mat files are OCV data points for silicon and graphite used by load_OCV_curves.m.

This code was written in MATLAB R2018a.

Any comments or questions please contact i.roper1993@gmail.com.
