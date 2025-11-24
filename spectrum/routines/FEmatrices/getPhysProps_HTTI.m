function ElPhysProps = getPhysProps_HTTI(CompStruct,BasicMatrices,...
                FEMatrices,ii_d,ii_el)

%   Part of the toolbox for solving problems of wave propagation
% in arbitrary anisotropic inhomogeneous waveguides.
% For details see User manual 
% and comments to the main script gen_aniso.m
%
%   getPhysProps_HTTI function constructs the interpolation for the elastic moduli
%  and the density inside the element. It computes the matrices, which indicates the expansion
%  and the coefficients of the expansion of C_ij matrices' elements and Rho
%  into the intepolating functions L_j (L1, L2, L3).
%  For HTTI solid.
% The expansion is of the form C_ij (:::) = sum_n C_ij(n)*N_n(:::)
% For each node of the element the value is the same - c_ij{ii_l}
% The expansion is of the form Rho (:::) = sum_n Rho(n)*N_n(:::)
% For each node of the element the value is the same - rho{ii_l}
%
%  Current implementation is for cubic elements.
% NB!!! This implementation is specific for FE computations
%
%   [T.Zharnikov, SMR v0.3_08.2014]
%
% function ElPhysProps = getPhysProps_HTTI(ii_l,BasicMatrices,...
%                DEMeshProps,ii_el)
%
%  Inputs - 
%
%       CompStruct - structure containing parameters of the model
%
%  Outputs -
%
%       Pos - structure array, indicating positions of the blocks corresponding
%               to the various hierarchy levels (layers, harmonics, variables) 
%               inside the full matrix representation matrix.
%
%  M-files required-
%
% Last Modified by Timur Zharnikov SMR v0.3_08.2014

%###############################################################################
%
%   Code for getPhysProps_HTTI
%
%###############################################################################
%===============================================================================
% Initialization
%===============================================================================

% Identify number of nodes 
N_nodes = CompStruct.Advanced.N_nodes;

% Allocate memory for the matrices.
% The implementation is for cubic elements, hence expansion is up to the
% 3rd power in L_i

ElPhysProps.CijMatrix = zeros(N_nodes,6,6); 
ElPhysProps.RhoVec = zeros(N_nodes,1); 

% retrieve elastic moduli tensor and the density for the domain
Cij_6x6 = FEMatrices.PhysProp{ii_d}.c_ij.*1.e9;
RhoValue = FEMatrices.PhysProp{ii_d}.rho*1.e3;

%===============================================================================
% Compute the elements of CijMatrix expansion into intepolating functions
% The expansion is of the form C_ij (:::) = sum_n C_ij(n)*N_n(:::)
% For each node of the element the value is the same - c_ij{ii_l}
% The expansion is of the form Rho (:::) = sum_n Rho(n)*N_n(:::)
% For each node of the element the value is the same - rho{ii_l}
%===============================================================================

% process all of the nodes of the element
for ii = 1:N_nodes
    ElPhysProps.RhoVec(ii) = RhoValue;
    
    % process all of the elastic moduli tensor components
    ElPhysProps.CijMatrix(ii,:,:) = Cij_6x6;
%     for aa = 1:6
%         for bb = 1:6
%             ElPhysProps.CijMatrix(ii,aa,bb) = Cij_6x6(aa,bb);
%         end;
%     end;
end;
    
end
