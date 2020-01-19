# stress-assisted-diffusion-linear-model
This programme solves the model and produce the figures in the research paper 'The effect of mechanical stress on lithium distribution and geometry optimisation for multi-material lithium-ion anodes'. It was written in 2019 using MATLAB R2018a.

To produce the plots as seen in the paper, add '%' to the start of the line under each heading in the file Produce_plots.m to uncomment different sections of the code. This is so that not everything is calculated each time the program is run. Executing the Produce_plots.m script runs the program.

Constants for silicon and graphite (example materials used in the paper) are set in constants.m and the OCV data for each material is suppled in a .mat files (Silicon_points_interp.mat and Graphite_points_interp.mat). The materials used in the model can be changed by changing the mechanical and chemical constants in constants.m and by supplying the OCV data in a new .mat file and change the core_OCV_data and shell_OCV_data parameters in Produce_plots.m.

The parameters used in each plot (e.g. the core volume when calculating the effective OCV) can be altered at the top of each section in Produce_plots.m. The functions used to calculate the results to be plotted are then called, and any equations used are labelled in reference to the equation number in the paper.

The other .m files in this repository are functions used by Produce_plots.m and all include a description of the function, the inputs and the outputs at the beginning of the .m file.

Any comments or questions please contact i.roper1993@gmail.com.
