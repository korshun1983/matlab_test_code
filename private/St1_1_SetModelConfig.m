function [InputParam] = St1_1_SetModelConfig()

%   Part of the toolbox for solving problems of wave propagation
% in arbitrary anisotropic inhomogeneous waveguides.
% For details see User manual 
% and comments to the main script gen_aniso.m
%
%   St1_1_SetModelConfig.m M-file      
%      St1_1_SetModelConfig.m sets the configuration of toolbox,
% e.g. default values for the parameters, 
% which contain general information, such as type of the problem 
% (spectrum, source, etc.), method used (spectral method, SAFEM, 
% matrix Riccati, etc.), etc.
% The parameters can be modified to define the model 
% of the user's choice.
%
% NB! Typically this script will not be modified by the user.
%
%   [T.Zharnikov, D.Syresin, SMR v0.12_12.2012]
%
% function [] = St1_1_SetModelConfig(InputParam)
%
%  Inputs - 
%   no inputs for this script
%
%  Outputs -
%       InputParam - structure containing Config structure 
%                    with the configuration parameters of the model.
%
%  M-files required-
%
% Last Modified by Timur Zharnikov SMR v0.1_01.2012

%###############################################################################
%
%   Code for St1_1_SetModelConfig
%
%###############################################################################
%===============================================================================
% Initialization
%===============================================================================

%===============================================================================
% Setting the configuration parameters
%===============================================================================

% Choice of the problem to be solved 

InputParam.Config.ProblemType = 'spectrum'; 
%   types of problems: 
% - 'spectrum' - find the spectrum of the waveguide;
% - 'source' - compute the problem with the source;
% - 'ExFun' - compute the excitation function;
% - etc.

% Choice of the method 

InputParam.Config.NumMethod = 'SAFE'; 
%   types of methods: 
% - 'SM' - spectral method;
% - 'SAFE' - semianalytical FEM;
% - 'Riccati' - matrix Riccati approach;
% - etc. (e.g., route search)

% Choice of the method 

InputParam.Config.SpeedUp = 'no'; 
%   types of speed up: 
% - 'no' - standard method - no sparsity, no eigs, no parallelization;
% - 'yes' - try to speed up the computation - sparsity, 'eigs', parallel;
% - etc. (e.g., route search)
%
% N.B.!!! In the case of SAFE the computations are dones automatically on
% sparse matrices, so the SpeedUp option should be put to 'no'.

% Choice of whether to save intermediate data to disk 

InputParam.Config.SaveData = 'yes'; 
% 'yes' - save data, requires more time, necessary for the development;
% 'no' - no intermediate data, faster but no intermediate control

% Choice of the BC at outer boundary 

InputParam.Config.OuterBC = 'fixed'; 
%   types of conditions: 
% - 'fixed' - fixed outer radius for all range of parameters;
% - 'adjust' - adjust outer radius according to the frequency 
%       and wavelength;
% - 'PML' - PML boundary conditions (perfectly matched layer);
% - 'ABC' - ABC boundary conditions (absorbing boundary conditions);
% - 'TTBC' - truncated transparent boundary conditions;
% - etc.

% Choice of the PML at outer boundary 

InputParam.Config.PML = 'r2'; 
%   types of conditions: 
% - 'none' - no PML;
% - 'r2' - gamma = 1 + 1i*a*r^2/omega;
% - etc.

% Choice of eccentricity option 

InputParam.Config.Eccentricity = 'no'; 
% 'yes' - present, 'no' - absent

% Choice of the symmetry

InputParam.Config.Symmetry = 'none'; 
%   types of variables: 
% - 'none' - no symmetry, compute full series of basis functions;
% - 'plane0' - plane of mirror symmetry at theta=0;
% - 'plane_pi2' - plane of mirror symmetry at theta=pi/2;
% - etc. 

% Choice of the eigenvalue variable 

InputParam.Config.EigenVar = 'k'; 
%   types of variables: 
% - 'omega' - frequency;
% - 'k' - wavevector (e.g., necessary for viscoelasticity);
% - etc. 

InputParam.Config.SloUnits = 'us/ft'; 
%   types of units to use for slowness input and output: 
% - 'us/ft' - mus/ft; - use this one for now
% - 'us/m' - mus/m;

InputParam.Config.FreqUnits = 'kHz'; 
%   types of units to use for frequency input and output: 
% - 'kHz' - kHz; - use this one for now
% - 'Hz' - Hz;

InputParam.Config.PressureUnits = 'GPa'; 
%   types of units to use for pressure and elastic moduli input and output: 
% - 'GPa' - GPa; - use this one for now
% - 'Pa' - Pa;

InputParam.Config.CheckAsymptote = 'yes'; 
%   'yes' - check and compare against know asymptotes, 
% like mode search, anisotropic Stoneley, etc.
%   'no' - do not check against asymptotes

InputParam.Config.DisplayAttenuation = 'yes'; 
%   'yes' - display attenuation plot
%   'no' - do not display attenuation plot

end