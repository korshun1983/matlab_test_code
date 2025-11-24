function [FEMatrices] = ...
        ICMatrices_ff_ss_SAFE_cubic(CompStruct,BasicMatrices,FEMatrices,ii_int,ii_d1,ii_d2)
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

% identify, what is the type of the first domain - fluid or solid
% put fluid domain first for the computation of the IC matrices and later
% adjust accordingly

D1NodesNum = size(FEMatrices.DNodes{ii_d1},2);
D1FEVarNum = CompStruct.Data.DVarNum(ii_d1)*D1NodesNum;
D2NodesNum = size(FEMatrices.DNodes{ii_d2},2);
D2FEVarNum = CompStruct.Data.DVarNum(ii_d2)*D2NodesNum;

% create zero sparse matrix for the 
% fluid-solid boundary condition matrix
FEMatrices.ZeroD12 = sparse(D1FEVarNum,D2FEVarNum);
FEMatrices.ZeroD21 = sparse(D1FEVarNum,D2FEVarNum);

% prepare interface elements
BElements = FEMatrices.BoundaryEdges(1:2,ismember(FEMatrices.BoundaryEdges(3,:),ii_int));
NBElements = size(BElements,2);

BNodesFull = [];

%===============================================================================
% Process all of the edges, which belong to the boundary
%===============================================================================

for ii_ed = 1:NBElements
    
    % find the start and the end nodes of the edge
    % depending on whether the fluid domain is the first or the last,
    % the orientation of the edge changes.
    % In order to identify the proper mid-side nodes, the order of the edge
    % nodes should be adjusted
    EdgeNodes = [];
    EdgeNodes(1) = BElements(1,ii_ed); % the first node of the boundary element (edge)
    EdgeNodes(2) = BElements(2,ii_ed); % the second node of the boundary element (edge)
    
    % find the elements, which contain these nodes
    FNodes = ismember(FEMatrices.DElements{ii_d1}(1:10,:), EdgeNodes(1:2));
    NNodes = sum(FNodes);
    [~,D1EdgeEl] = max(NNodes);
    
    % read the nodes, which belong to the adjoint (to the edge) elements
    TriNodes = FEMatrices.DElements{ii_d1}(1:10,D1EdgeEl);
    
    % find the nodes in between (mid-side nodes for cubic interpolation)
    EdgeNodesPos.D1(1) = find(ismember(TriNodes,EdgeNodes(1)));
    EdgeNodesPos.D1(2) = find(ismember(TriNodes,EdgeNodes(2)));
    % find the interior nodes of the edge element
    if ( EdgeNodesPos.D1(2) > mod(EdgeNodesPos.D1(1),3) )
        % find the third and the fourth nodes
        EdgeNodesPos.D1(3) = ( EdgeNodesPos.D1(1) + 1 )*2;
        EdgeNodesPos.D1(4) = ( EdgeNodesPos.D1(1) + 1 )*2 + 1;
    else
        % find the third and the fourth nodes
        EdgeNodesPos.D1(3) = ( EdgeNodesPos.D1(2) + 1 )*2 + 1;
        EdgeNodesPos.D1(4) = ( EdgeNodesPos.D1(2) + 1 )*2;
    end;
    EdgeNodes(3) = TriNodes(EdgeNodesPos.D1(3));
    EdgeNodes(4) = TriNodes(EdgeNodesPos.D1(4));
    
    BNodesFull = [BNodesFull EdgeNodes];
        
end;

%===============================================================================
% Assemble the full list of the nodes, which belong to the interface
%===============================================================================

[FEMatrices.BNodesFull{ii_int}] = sort(unique(BNodesFull));

end