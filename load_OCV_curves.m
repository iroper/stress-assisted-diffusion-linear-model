function [OCV_1, OCV_2, mu_SF_1, mu_SF_2] = load_OCV_curves(OCV_data_1, OCV_data_2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Imports the OCV data measured as volts against some measure of SOC 
% (either c/c_max or mAh/g) and rescales the x-axis to be between zero and 
% one. Also creates nondimensional interpolation functions of the
% stress-free chemical potential of both materials.
%
% Inputs
% ------
% OCV_data_1 : string of the name of the .mat file containing the x and y
%              values of the OCV data (voltage on the y axis, x axis
%              irrelevant) for the core material. Should only have one 
%              field.
% OCV_data_2 : string of the name of the .mat file containing the x and y
%              values of the OCV data (voltage on the y axis, x axis
%              irrelevant) for the core material. Should only have one 
%              field.
%
% Outputs
% -------
% OCV_1   : OCV data points for the core material, with a normalised x-axis
%           such that it is between zero and one.
% OCV_2   : OCV data points for the shell material, with a normalised 
%           x-axis such that it is between zero and one.
% mu_SF_1 : Interpolation function of the stress-free chemical potential of
%           the core material against state of charge.
% mu_SF_2 : Interpolation function of the stress-free chemical potential of
%           the shell material against state of charge.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global fund_charge R_g T

% Load the data
load_OCV_data_1 = load(OCV_data_1);
load_OCV_data_2 = load(OCV_data_2);
% Extract field names
StructNames1 = fieldnames(load_OCV_data_1);
StructNames2 = fieldnames(load_OCV_data_2);
% Convert fields into matrices
eval(['OCV_1 = load_OCV_data_1.' StructNames1{1} ';']);
eval(['OCV_2 = load_OCV_data_2.' StructNames2{1} ';']);
% Scale the data so both concentrations are in [0,1]
OCV_1(:,1) = (OCV_1(:,1)-min(OCV_1(:,1)))/(max(OCV_1(:,1)) - min(OCV_1(:,1)));
OCV_2(:,1) = (OCV_2(:,1)-min(OCV_2(:,1)))/(max(OCV_2(:,1)) - min(OCV_2(:,1)));

% Define non-dimensional stress-free chemical potential interpolation
% functions
mu_SF_1 = @(x) interp1(OCV_1(:,1), OCV_1(:,2)*((-fund_charge)/(R_g*T)), x, 'linear'); % Eq. 20
mu_SF_2 = @(x) interp1(OCV_2(:,1), OCV_2(:,2)*((-fund_charge)/(R_g*T)), x, 'linear'); % Eq. 20