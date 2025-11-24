function [PhysProp] = PreparePhysProp_HTTI_sp_SAFE(CompStruct,ii_l)

%   Part of the toolbox for solving problems of wave propagation
% in arbitrary anisotropic inhomogeneous waveguides.
% For details see User manual 
% and comments to the main script gen_aniso.m
%
% PreparePhysProp_HTTI_sp_SAFE M-file      
%      PreparePhysProp_HTTI_sp_SAFE, by itself, prepares structure 
%      with physical properties of the homogeneous TTI (HTTI) layer.
%      It computes Fourier components of the elastic moduli tensor
%      in the borehole coordinate system.
%      The computation assumes that borehole is vertical and VTI 
%      axis is rotated to angle phi in horizontal plane 
%      and then inclined to angle theta. 
%      For this reason back rotation is employed here.
%
%   [T.Zharnikov, D.Syresin, SMR v0.3_08.2014]
%
% function [PhysProp] = PreparePhysProp_HTTI_sp_SAFE(CompStruct,ii_l)
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
%   Code for PreparePhysProp_HTTI_sp_SAFE
%
%###############################################################################

% %==============================================================================
% % DS!!!
% %===============================================================================
% for i=1:size(Input.SolidSubdomains,2)
%     layer_num=Input.SolidSubdomains(i);
%     Inc_angle=Input.Elastic_prop(layer_num,6);
%     
%     C11=Input.Elastic_prop(layer_num,1);
%     C33=Input.Elastic_prop(layer_num,2);
%     C44=Input.Elastic_prop(layer_num,3);
%     C66=Input.Elastic_prop(layer_num,4);
%     C13=Input.Elastic_prop(layer_num,5);
%     C12=C11-2*C66;    
%     
%     C_VTI=[C11 C12 C13 0 0 0; C12 C11 C13 0 0 0; C13 C13 C33 0 0 0; 0 0 0 C44 0 0; 0 0 0 0 C44 0; 0 0 0 0 0 C66];
%     C_TTI=C_tensor_rotation(C_VTI,Inc_angle );
%      
%          
%     % Use only 5 meanining values (increse speed of eigs);
%     power=fix(log10(max(max(abs(C_TTI)))));
%     C_TTI=round(C_TTI/10^(power-6))*10^(power-6);
%     Input.C(Input.SolidSubdomains(i),:,:) = C_TTI;
% end

%===============================================================================
% Extracting VTI layer parameters and inclination
%===============================================================================

PhysProp.rho = CompStruct.Model.DomainParam{ii_l}(1);
PhysProp.c_VTI = CompStruct.Model.DomainParam{ii_l}(2:6);
PhysProp.theta = CompStruct.Model.DomainParam{ii_l}(7);
if max(size(CompStruct.Model.DomainParam{ii_l})) == 8
    PhysProp.phi = CompStruct.Model.DomainParam{ii_l}(8);
else
    PhysProp.phi = 0;
end;

%===============================================================================
% Computing rotated elastic moduli tensor
%===============================================================================

    [c_ij c_ijkl] = CompStruct.Methods.em_tensor_VTI(PhysProp.c_VTI); 
    [rot_m] = CompStruct.Methods.rot_matrix(PhysProp.theta,PhysProp.phi);
    % rotating back to get elastic moduli tensor in the borehole system
    [c_ij_rot c_ij_rot_back] = CompStruct.Methods.rot_c_ij(c_ij, rot_m); 
    
%     % Use only 5 leading digits (increse speed of eigs) - suggested by DS;
%     power = fix(log10(max(max(abs(c_ij_rot_back)))));
%     c_ij_rot_back = round(c_ij_rot_back/10^(power-6))*10^(power-6);

    PhysProp.c_ij = c_ij_rot_back;
    
% %===============================================================================
% % Assembling radial density distibution.
% %   In case of homogeneity radial dependence is absent.
% %===============================================================================
% PhysProp.rho_dist = transpose(linspace(PhysProp.rho,PhysProp.rho,CompStruct.Model.LayerNr(ii_l)));

end