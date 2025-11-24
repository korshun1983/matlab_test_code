function [FEMatrices,FullMatrices] = ...
        AssembleFullMatrices_fs_SAFE_cubic(CompStruct,BasicMatrices,FEMatrices,FullMatrices,ii_int,ii_d1,ii_d2)
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

% D1NodesNum = size(FEMatrices.DNodes{ii_d1},2);
% D1FEVarNum = CompStruct.Data.DVarNum(ii_d1)*D1NodesNum;
% D2NodesNum = size(FEMatrices.DNodes{ii_d2},2);
% D2FEVarNum = CompStruct.Data.DVarNum(ii_d2)*D2NodesNum;
% % create zero sparse matrix for the 
% % fluid-solid boundary condition matrix
% FEMatrices.ZeroD12 = sparse(BFEVarNum1,BFEVarNum2);
% FEMatrices.ZeroD21 = sparse(BFEVarNum1,BFEVarNum2);

%===============================================================================
% Assembling full matrices by adding the mass and stiffness matrices of the
% next block
% and by inserting the matrices, which are responsible for the boundary
% conditions
%===============================================================================

% % example of script for inserting the new matrix blocks using sparse
% % matrices procedures (potentially can increase the speed)
% FullMatrices.K1Matrix = sparse(CompStruct.Data.DVarNum{N_domain}(2),CompStruct.Data.DVarNum{N_domain}(2));
% [ DVecRepRow, DVecRepCol, DVecRepV ] = find(K1Matrix_d{ii_d});
% DVecRepRow = DVecRepRow + CompStruct.Data.DVarNum{ii_d}(1) - 1;
% DVecRepCol = DVecRepCol + CompStruct.Data.DVarNum{ii_d}(1) - 1;
% DSize = size(K1Matrix_d{ii_d});
% K1Matrix = K1Matrix + sparse(DVecRepRow, DVecRepCol, DVecRepV, DSize(1), DSize(2));

% adding the matrices for the next block
FullMatrices.K1Matrix = blkdiag(FullMatrices.K1Matrix,FEMatrices.K1Matrix_d{ii_d2});
FullMatrices.K2Matrix = blkdiag(FullMatrices.K2Matrix,FEMatrices.K2Matrix_d{ii_d2});
FullMatrices.K3Matrix = blkdiag(FullMatrices.K3Matrix,FEMatrices.K3Matrix_d{ii_d2});
FullMatrices.MMatrix = blkdiag(FullMatrices.MMatrix,FEMatrices.MMatrix_d{ii_d2});
curFullMatSize = size(FullMatrices.PMatrix);
curMatD12Size = size(FEMatrices.PMatrixD12{ii_int});
curZeroMat12Size = [ ( curFullMatSize(1) - curMatD12Size(1) ), curMatD12Size(2) ];
ZeroMatrix12 = sparse(curZeroMat12Size(1),curZeroMat12Size(2));
InsertPMatrixD12 = [ ZeroMatrix12 ; FEMatrices.PMatrixD12{ii_int} ];
curMatD21Size = size(FEMatrices.PMatrixD21{ii_int});
curZeroMat21Size = [ curMatD21Size(1), ( curFullMatSize(2) - curMatD21Size(2) ) ];
ZeroMatrix21 = sparse(curZeroMat21Size(1),curZeroMat21Size(2));
InsertPMatrixD21 = [ ZeroMatrix21 FEMatrices.PMatrixD21{ii_int} ];
FullMatrices.PMatrix = [ FullMatrices.PMatrix	InsertPMatrixD12; ...
                         InsertPMatrixD21       FEMatrices.PMatrix_d{ii_d2}];
%FullMatrices.PMatrix = FullMatrices.PMatrix.*0.; % wrong results !!!

% identify nodes, which should be removed from the computation

FEMatrices.DNodesRem{ii_d2} = [ FEMatrices.DNodesRem{ii_d2} ];
FEMatrices.DNodesRem{ii_d2} = sort(unique(FEMatrices.DNodesRem{ii_d2}));
RemoveNodesPos = ismember(FEMatrices.DNodes{ii_d2},FEMatrices.DNodesRem{ii_d2});
FEMatrices.DNodesComp{ii_d2} = FEMatrices.DNodes{ii_d2}(~RemoveNodesPos);
                     
%===============================================================================
% Adjust the list of variables (remove duplicates, if any, etc.)
% and the positions of the variables in the full vector of variables
%===============================================================================

% DfNodesNum = size(FEMatrices.DNodes{ii_df},2);
% DfFEVarNum = DfNodesNum;
% DsNodesNum = size(FEMatrices.DNodes{ii_ds},2);
% DsFEVarNum = CompStruct.Data.DVarNum(ii_ds)*DsNodesNum;
% BFEVarNumFluid = DfFEVarNum;
% BFEVarNumSolid = DsFEVarNum;
% 
% % identify positions of the element nodes and the variables in the
% % respective matrices
% for ii_n = 1:4
%     DfNodePos = find(ismember(FEMatrices.DNodes{ii_df},EdgeNodes(ii_n)));
%     DsNodePos = find(ismember(FEMatrices.DNodes{ii_ds},EdgeNodes(ii_n)));
%     VarRowArr( ( CompStruct.Data.DVarNum(ii_df)*( ii_n - 1 ) + 1 ):( CompStruct.Data.DVarNum(ii_df)*ii_n ) ) = ...
%         [( CompStruct.Data.DVarNum(ii_df)*( DfNodePos - 1 ) + 1 ):( CompStruct.Data.DVarNum(ii_df)*DfNodePos )];
%     VarColArr( ( CompStruct.Data.DVarNum(ii_ds)*( ii_n - 1 ) + 1 ):( CompStruct.Data.DVarNum(ii_ds)*ii_n ) ) = ...
%         [( CompStruct.Data.DVarNum(ii_ds)*( DsNodePos - 1 ) + 1 ):( CompStruct.Data.DVarNum(ii_ds)*DsNodePos )];
% end;

end

