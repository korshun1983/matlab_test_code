function [InputParam] = St1_SetModel()

%   Part of the toolbox for solving problems of wave propagation
% in arbitrary anisotropic inhomogeneous waveguides.
% For details see User manual 
% and comments to the main script gen_aniso.m
%
%   St1_SetModel.m M-file      
%      St1_SetModel.m sets default values for the parameters, 
% which describe the model. The parameters can be modified 
% to define the model of the user's choice.
% For clarity and convenience, the definition is broken into 
% four parts:
%
%   - configuration - contains general information, 
% such as type of the problem (spectrum, source, etc.),
% method used (spectral method, SAFEM, matrix Riccati, etc.), etc.
%
%   - methods - contains the handles of the methods to solve the problem
%
%   - user - parameters, which are supposed to be available to the user,
% such as model geometry, types of layers, interface conditions, 
% frequency range, etc.
%
%   - advanced - parameters, which are not obvious, less understood,
% and generally which require more care in handling, such as
% the number of approximation points, outer radius of the model, etc.
%
%   [T.Zharnikov, D.Syresin, SMR v0.12_12.2012]
%
% function [InputParam] = St1_SetModel()
%
%  Inputs - 
%   no inputs for this script
%
%  Outputs -
%       InputParam - structure containing input parameters 
%                    for the computation. The details are not listed
%                    because it can vary depending on the tool
%                    configuration
%
%  M-files required-
%
% Last Modified by Timur Zharnikov SMR v0.1_01.2012

%###############################################################################
%
%   Code for St1_SetModel
%
%###############################################################################
%===============================================================================
% Initialization
%===============================================================================

%===============================================================================
% Setting the configuration of the problem and method 
% Setting parameters of the model 
%===============================================================================

%   Set the configuration parameters
% such as type of the problem (spectrum, source, etc.),
% method used (spectral method, SAFEM, matrix Riccati, etc.), etc.
% NB! Typically not for the user modification
[InputParam] = St1_1_SetModelConfig(); 
%InputParam = St1_1_SetModelConfig; 

%   Set the methods, which will be used to define and solve
% the model, which is specified by the configuration.
% NB! Typically not for the user modification
[InputParam] = St1_2_PrepareModelMethods(InputParam);

%   Set the user model parameters 
% such as model geometry, types of layers, interface conditions, 
% frequency range, etc.
% NB! Typically the parameters, which are defined by the user.
[InputParam] = InputParam.Methods.St1_3_SetModelUser(InputParam); 
InputParam.Model.N_disp=length(InputParam.Model.f_array); % number of frequencies

%   Set the advanced parameters 
% such as the number of approximation points, 
% outer radius of the model, etc.
% NB! Typically not for the user modification
[InputParam] = InputParam.Methods.St1_4_SetModelAdvanced(InputParam); 
%InputParam = St1_3_SetModelAdvanced; 

end
