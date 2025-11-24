function [FEMatrices] = MatricesParts_HTTI_sp_SAFE_cubic(CompStruct,...
            BasicMatrices,FEMatrices,ii_d)

%         CompStruct.Methods.MatricesParts_sp_SAFE{ii_d}(FEMatrices.DElements{ii_d},...
%             FEMatrices.DEMeshProps{ii_d},FEMatrices.DNodes{ii_d},...
%             CompStruct,BasicMatrices,ii_d); %,integration_acc);
        
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

NDomainElements = size(FEMatrices.DElements{ii_d},2);
DNodesNum = size(FEMatrices.DNodes{ii_d},2);
DFEVarNum = CompStruct.Data.DVarNum(ii_d)*DNodesNum;
VarVecArr = zeros(1,CompStruct.Data.DVarNum(ii_d)*CompStruct.Advanced.N_nodes);

% create zero sparse matrices
FEMatrices.K1Matrix_d{ii_d} = sparse(DFEVarNum,DFEVarNum);
FEMatrices.K2Matrix_d{ii_d} = sparse(DFEVarNum,DFEVarNum);
FEMatrices.B1tCB2Matrix_d{ii_d} = sparse(DFEVarNum,DFEVarNum);
FEMatrices.B2tCB1Matrix_d{ii_d} = sparse(DFEVarNum,DFEVarNum);
FEMatrices.K3Matrix_d{ii_d} = sparse(DFEVarNum,DFEVarNum);
FEMatrices.MMatrix_d{ii_d} = sparse(DFEVarNum,DFEVarNum);
FEMatrices.PMatrix_d{ii_d} = sparse(DFEVarNum,DFEVarNum);

% el_x_nodes = zeros(1,3);
% el_y_nodes = zeros(1,3);
% el_x_nodes_full = zeros(1,10);
% el_y_nodes_full = zeros(1,10);

% % define the array with the numbers of nodes, belonging to the domain
% DNodesEl = DElements(1:10,:);
% % transform it into the vector, sort, and find the nodes, belonging to the domain 
% DNodesVec = reshape(DNodesEl,1,[]);
% % find the indices (ordered) of the nodes, which belong to the domain
% DNodes = sort(unique(DNodesVec));
% % OrderDNodesVec = sort(DNodesVec);
% % DiffDNodesVec = diff(OrderDNodesVec);
% % DNodesIndex = [1, ( 1 + find(DiffDNodesVec) )];
% % % find the indices (ordered) of the nodes, which belong to the domain
% % DNodes = OrderDNodesVec(NodesIndex);

%===============================================================================
% Process all of the elements, which belong to the domain
% Compute the K_i and M matrices for the particular element
%===============================================================================

for ii_el = 1:NDomainElements
    
    % read the numbers of nodes, which belong to the element
    TriNodes = FEMatrices.DElements{ii_d}(1:10,ii_el);
    
    % identify positions of the element nodes and the variables in the
    % respective matrices
    for ii_n = 1:10
        DNodePos = find(ismember(FEMatrices.DNodes{ii_d},TriNodes(ii_n)));
        VarVecArr(( CompStruct.Data.DVarNum(ii_d)*( ii_n - 1 ) + 1 ):...
            (CompStruct.Data.DVarNum(ii_d)*ii_n)) = ...
                [ (CompStruct.Data.DVarNum(ii_d)*( DNodePos - 1 ) + 1 ):...
                  (CompStruct.Data.DVarNum(ii_d)*DNodePos) ];
    end;

%     % read x and coordinates of the corner nodes in the element
%     el_x_nodes = MeshNodes(1,MeshTri(1:3,ii_el));
%     el_y_nodes = MeshNodes(2,MeshTri(1:3,ii_el));  
%     el_x_nodes_full = MeshNodes(1,MeshTri(1:10,ii_el));
%     el_y_nodes_full = MeshNodes(2,MeshTri(1:10,ii_el));  

    % read physical properties of the element
    % and prepare the matrices for the expansion of these elements
    % into interpolating functions L_i
    ElPhysProps = CompStruct.Methods.getPhysProps{ii_d}...
        (CompStruct,BasicMatrices,FEMatrices,ii_d,ii_el);
    
    % read the properties of the element (a, b, c, delta as defined in 
    % the book of Zienkiewicz (v.1 p.181))
    TriProps.a = FEMatrices.DEMeshProps{ii_d}.a(:,ii_el);
    TriProps.b = FEMatrices.DEMeshProps{ii_d}.b(:,ii_el);
    TriProps.c = FEMatrices.DEMeshProps{ii_d}.c(:,ii_el);
    TriProps.delta = FEMatrices.DEMeshProps{ii_d}.delta(:,ii_el);

    % Compute the vectors, which enter into the expressions for the 
    % matrices of derivatives of the shape functions
    % dxNLMatrix(nN,ii,jj,kk) = d/dxL(1)*dNLMatrix(1,nN,ii,jj,kk) +
    % d/dxL(2)*dNLMatrix(2,nN,ii,jj,kk) + d/dxL(3)*dNLMatrix(3,nN,ii,jj,kk)
    % d/dx L_i = 1/(2*delta) * b_i
    % the factor 1/delta is omitted for the simplification of computations
    TriProps.dxL = 1/2*TriProps.b;
    TriProps.dyL = 1/2*TriProps.c;
    
    % compute the contributions of the element to the stiffness and mass
    % matrices
    [ElMatrices] = CompStruct.Methods.KM_matrix{ii_d}(BasicMatrices, CompStruct, FEMatrices, ...
                       ii_d, ii_el, ElPhysProps, TriProps);
                         
    % prepare sparse matrices to insert into the matrices for the domain
    % compute rows and columns, to which to insert
    VarRows = repmat(VarVecArr,1,30);
    VarCols = reshape(repmat(VarVecArr,30,1),1,[]);
    % create vectors of elements, which to insert
    strK1Matrix_el = reshape(ElMatrices.K1Matrix,1,[]);
    strK2Matrix_el = reshape(ElMatrices.K2Matrix,1,[]);
    strB1tCB2Matrix_el = reshape(ElMatrices.B1tCB2Matrix,1,[]);
    strB2tCB1Matrix_el = reshape(ElMatrices.B2tCB1Matrix,1,[]);
    strK3Matrix_el = reshape(ElMatrices.K3Matrix,1,[]);
    strMMatrix_el = reshape(ElMatrices.MMatrix,1,[]);
    % create sparse matrices
    spK1Matrix_el = sparse(VarRows,VarCols,strK1Matrix_el,DFEVarNum,DFEVarNum);
    spK2Matrix_el = sparse(VarRows,VarCols,strK2Matrix_el,DFEVarNum,DFEVarNum);
    spB1tCB2Matrix_el = sparse(VarRows,VarCols,strB1tCB2Matrix_el,DFEVarNum,DFEVarNum);
    spB2tCB1Matrix_el = sparse(VarRows,VarCols,strB2tCB1Matrix_el,DFEVarNum,DFEVarNum);
    spK3Matrix_el = sparse(VarRows,VarCols,strK3Matrix_el,DFEVarNum,DFEVarNum);
    spMMatrix_el = sparse(VarRows,VarCols,strMMatrix_el,DFEVarNum,DFEVarNum);
                         
    % insert the matrices for the element into the matrices for the domain
    FEMatrices.K1Matrix_d{ii_d} = FEMatrices.K1Matrix_d{ii_d} + spK1Matrix_el;
    FEMatrices.K2Matrix_d{ii_d} = FEMatrices.K2Matrix_d{ii_d} + spK2Matrix_el;
    FEMatrices.B1tCB2Matrix_d{ii_d} = FEMatrices.B1tCB2Matrix_d{ii_d} + spB1tCB2Matrix_el;
    FEMatrices.B2tCB1Matrix_d{ii_d} = FEMatrices.B2tCB1Matrix_d{ii_d} + spB2tCB1Matrix_el;
    FEMatrices.K3Matrix_d{ii_d} = FEMatrices.K3Matrix_d{ii_d} + spK3Matrix_el;
    FEMatrices.MMatrix_d{ii_d} = FEMatrices.MMatrix_d{ii_d} + spMMatrix_el;
                         
end;

%===============================================================================
% Make the K_i and M matrices sparse - for the speed and the storage space
%===============================================================================

% K1Matrix_d = sparse(K1Matrix_d);
% K2Matrix_d = sparse(K2Matrix_d);
% K3Matrix_d = sparse(K3Matrix_d);
% MMatrix_d = sparse(MMatrix_d);

end

