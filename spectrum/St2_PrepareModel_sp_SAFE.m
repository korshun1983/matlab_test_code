function [CompStruct] = St2_PrepareModel_sp_SAFE(InputParam)
      
%   Part of the toolbox for solving problems of wave propagation
% in arbitrary anisotropic inhomogeneous waveguides.
% For details see User manual 
% and comments to the main script gen_aniso.m
%
%   St2_PrepareModel_sp_SM.m M-file      
%      St2_PrepareModel_sp_SM.m prepares the specific set of the parameters and 
% the methods, which will be used for the computations.
% This preparation and the selection of the parameters and methods
% is based on the model and the problem, selected as the input.
% For clarity and convenience, the definition is broken into parts:
%
%   - data preparation -  organizes input data in the form,
% which is convenient for particular problems
%
%   - methods selection - defines handles to the function,
% which will be used to prepare the data, solve the problem, etc.
%
% NB!!! Specific implementation of this step is chosen at the
% step of methods selection.
%
%   [T.Zharnikov, D.Syresin, SMR v0.12_12.2012]
%
% function [CompStruct] = St2_PrepareModel_sp_SM(InputParam)
%
%  Inputs - 
%
%       InputParam - structure containing initial parameters of the model.
%
%  Outputs -
%
%       CompStuct - structure containing the parameters of the model 
%                   in the form, which is convenient and ready for
%                   computations
%
%  M-files required-
%
% Last Modified by Timur Zharnikov SMR v0.1_01.2012

%###############################################################################
%
%   Code for St2_PrepareModel_sp_SM
%
%###############################################################################
%===============================================================================
% Initialization
%===============================================================================

% Retain all necessary information from Input structure
CompStruct = InputParam;

%===============================================================================
% Preparing structure for the computations
%===============================================================================

% Prepare input data for the computations. 
% Construct structure with parameters.
[CompStruct] = CompStruct.Methods.St2_1_PrepareModelParams(CompStruct);

% Assign the methods, which will be used to solve the problem
[CompStruct] = CompStruct.Methods.St2_2_PrepareModelMethods(CompStruct);

end