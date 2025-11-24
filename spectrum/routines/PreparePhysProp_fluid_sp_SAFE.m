function [PhysProp] = PreparePhysProp_fluid_sp_SAFE(CompStruct,ii_l)

%   Part of the toolbox for solving problems of wave propagation
% in arbitrary anisotropic inhomogeneous waveguides.
% For details see User manual 
% and comments to the main script gen_aniso.m
%
% PreparePhysProp_fluid_sp_SAFE M-file      
%      PreparePhysProp_fluid_sp_SAFE, by itself, prepares structure 
%      with physical properties of the homogeneous ideal fluid layer.
%      It computes Fourier components of the elastic moduli tensor
%      in the borehole coordinate system.
%      Since the fluid is ideal, elastic moduli are isotropic and the same 
%      in any coordinate system and there is no need to rotate it.
%
%   [T.Zharnikov SMR v0.3_08.2014]
%
% function [PhysProp] = PreparePhysProp_fluid_sp_SAFE(CompStruct,ii_l)
%
%  Inputs - 
%
%       CompStuct - structure containing the parameters of the model;
%
%       ii_l - number of the layer, for which the physical properties are
%              introduced
%
%  Outputs -
%
%       PhysProp - structure, containing physical properties of the ii_l-th
%       layer
%
%  M-files required-
%
% Last Modified by Timur Zharnikov SMR v0.3_08.2014

%###############################################################################
%
%   Code for PreparePhysProp_fluid_sp_SAFE
%
%###############################################################################

%===============================================================================
% Extracting ideal fluid layer parameters and inclination
%===============================================================================

PhysProp.rho = CompStruct.Model.DomainParam{ii_l}(1);
PhysProp.lambda = CompStruct.Model.DomainParam{ii_l}(2);

% %===============================================================================
% % Assembling radial density distibution.
% %   In case of homogeneity radial dependence is absent.
% %===============================================================================
% PhysProp.rho_dist = transpose(linspace(PhysProp.rho,PhysProp.rho,CompStruct.Model.LayerNr(ii_l)));

end