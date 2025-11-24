function [FEMatrices,FullMatrices] = ...
        RemoveRedundantVariables_SAFE(CompStruct,BasicMatrices,FEMatrices,FullMatrices)
% function [Fluid_bound_matrix,Solid_bound_matrix]=Fluid_solid_boundary_cubic( p,e,s,bnd_arr,Mud_subdomain)
%         %p nodes e edges s CompStruct, bnd_arr Mud_subdomain
        
%   Part of the toolbox for solving problems of wave propagation
% in arbitrary anisotropic inhomogeneous waveguides.
% For details see User manual 
% and comments to the main script gen_aniso.m
%
% MatricesParts_HTTI_sp_SAFE_cubic M-file      
%      MatricesParts_HTTI_sp_SAFE_cubic, by itself, 
%      prepares the blocks of the full stiffness and mass matrices 
%      for the particular domain.
%      Sparse matrix representation is used for speed and storage (the suggestion of Denis Syresin).
%      Ready and can handle inhomogeneous anisotropic solid case, but currently is used for the homogeneous case 
%       (more of an issue with the parameter input)
%
%   [T.Zharnikov SMR v0.3_09.2014]
%
% function [MMatrix_d,K1Matrix_d,K2Matrix_d,K3Matrix_d] = ...
%        MatricesParts_HTTI_sp_SAFE_cubic(DElements, DEMeshProps, DNodes,CompStruct)
%
%  Inputs - 
%
%       BasicMatrices - structure containing the basic matrices (differentiation, etc.);
%
%       CompStuct - structure containing the parameters of the model;
%
%       ii_l - number of the layer, for which the basic matrices are
%              assembled
%
%  Outputs -
%
%       L_parts - matrix representation for the governing operator
%
%       Srr_parts - matrix representation for the stress (s_rr) operator
%
%       D_parts - matrix representation for the displacement operator
%
%       M_parts - matrix representation for the left-hand side of the governing
%                   equations' operator
%
%  M-files required-
%
% Last Modified by Timur Zharnikov SMR v0.3_09.2014

%###############################################################################
%
%   Code for MatricesParts_HTTI_sp_SAFE_cubic
%
%###############################################################################

%===============================================================================
% Initialization
%===============================================================================

FullRemoveNodesArr = [];
FullRemoveVarArr = [];
FullVarArr = [1:BasicMatrices.Pos{CompStruct.Data.N_domain}(2)];

%===============================================================================
% Assemble array of the variables, which should be removed from 
% the problem 
%===============================================================================

% finding positions of the rows, which correspond to the nodes, 
% which will be removed
% NB! Nodes in DNodes and BNodesFull are ordered! 
% It is important condition for this script to work.

for ii_d = 1:CompStruct.Data.N_domain

    if ~isempty(FEMatrices.DNodesRem{ii_d})
%        FullRemoveNodesArr = [ FullRemoveNodesArr FEMatrices.DNodesRem{ii_d} ];
        
%         DNodesRemSize = size(FEMatrices.DNodesRem{ii_d});
%         DNodesRemSize = DNodesRemSize(2);
%         DRemoveVarPos = [];
%         for ii_n = 1:DNodesRemSize
%             DNodeRemPos = find(ismember(FEMatrices.DNodes{ii_d},FEMatrices.DNodesRem{ii_d}(ii_n)));
%             DRemoveVarPos(( CompStruct.Data.DVarNum(ii_d)*( ii_n - 1 ) + 1 ):...
%                 (CompStruct.Data.DVarNum(ii_d)*ii_n)) = ...
%                 BasicMatrices.Pos{ii_d}(1) - 1 + ...
%                 [ (CompStruct.Data.DVarNum(ii_d)*( DNodeRemPos - 1 ) + 1 ):...
%                 (CompStruct.Data.DVarNum(ii_d)*DNodeRemPos) ];
%         end;

        VarArr = [1:CompStruct.Data.DVarNum(ii_d)]';
        sizeVarArr = size(FEMatrices.DNodesRem{ii_d});
        LargeVarArr = repmat(VarArr,1,sizeVarArr(2));
        RemoveNodesPos = find(ismember(FEMatrices.DNodes{ii_d},FEMatrices.DNodesRem{ii_d}));        
        DRemoveVarPos = CompStruct.Data.DVarNum(ii_d)*( RemoveNodesPos - 1 );
        DRemoveVarPos = repmat(DRemoveVarPos,CompStruct.Data.DVarNum(ii_d),1);
        DRemoveVarPos = DRemoveVarPos + LargeVarArr;
        DRemoveVarPos = BasicMatrices.Pos{ii_d}(1) - 1 + reshape(DRemoveVarPos,1,[]);
        
        %     DVarPos = CompStruct.Data.DVarNum(ii_d)*( FEMatrices.DNodes{ii_d} - 1 );
        %     DVarPos = repmat(DVarPos,CompStruct.Data.DVarNum(ii_d),1);
        %     DVarPos = DVarPos + LargeVarArr;
        %     DVarPos = BasicMatrices.Pos{ii_d}(1) - 1 + reshape(DVarPos,1,[]);
        
        FullRemoveVarArr = [ FullRemoveVarArr DRemoveVarPos ];
        %     FullVarArr = [ FullVarArr DVarPos ];
    end;
    
end;

FullKeepVarArr = find(~ismember(FullVarArr,FullRemoveVarArr));

%===============================================================================
% Sum the rows and columns, which correspond to the coincident nodes
% in both domains
%===============================================================================

FullMatrices.K1Matrix = FullMatrices.K1Matrix(FullKeepVarArr,FullKeepVarArr);
FullMatrices.K2Matrix = FullMatrices.K2Matrix(FullKeepVarArr,FullKeepVarArr);
FullMatrices.K3Matrix = FullMatrices.K3Matrix(FullKeepVarArr,FullKeepVarArr);
FullMatrices.MMatrix = FullMatrices.MMatrix(FullKeepVarArr,FullKeepVarArr);
FullMatrices.PMatrix = FullMatrices.PMatrix(FullKeepVarArr,FullKeepVarArr);

end