function [FEMatrices] = ...
        ICMatrices_fluid_HTTI_SAFE_cubic(CompStruct,BasicMatrices,FEMatrices,ii_int,ii_d1,ii_d2)
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

D1NodesNum = size(FEMatrices.DNodes{ii_d1},2);
D1FEVarNum = CompStruct.Data.DVarNum(ii_d1)*D1NodesNum;
D2NodesNum = size(FEMatrices.DNodes{ii_d2},2);
D2FEVarNum = CompStruct.Data.DVarNum(ii_d2)*D2NodesNum;

% create zero sparse matrix for the 
% fluid-solid boundary condition matrix
FEMatrices.ZeroD12 = sparse(D1FEVarNum,D2FEVarNum);
FEMatrices.ZeroD21 = sparse(D2FEVarNum,D1FEVarNum);

% identify, what is the type of the first domain - fluid or solid
% put fluid domain first for the computation of the IC matrices and later
% adjust accordingly
switch CompStruct.Model.DomainType{ii_d1} 
    case 'fluid'
        ii_df = ii_d1;
        ii_ds = ii_d2;
    case 'HTTI'
        ii_df = ii_d2;
        ii_ds = ii_d1;
    otherwise
end;

% boundary elements in 2d formulation are edges 
BElements = FEMatrices.BoundaryEdges(1:2,ismember(FEMatrices.BoundaryEdges(3,:),ii_int));
NBElements = size(BElements,2);

BNodesFull = [];

% % the number of nodes should be multiplied by 3, because each edge
% % contains two more internal nodes, which are added during the
% % computation of the IC matrices (cubic elements)
% % for the closed boundary it results in the multiplication factor of 3
% % for the total number of boundary nodes
% BNodesNum = 3*size(FEMatrices.BNodes{ii_int},2);
% BFEVarNumFluid = BNodesNum;
% BFEVarNumSolid = 3*BNodesNum;
% VarVecArr = zeros(1,10);

DfNodesNum = size(FEMatrices.DNodes{ii_df},2);
DfFEVarNum = DfNodesNum;
DsNodesNum = size(FEMatrices.DNodes{ii_ds},2);
DsFEVarNum = CompStruct.Data.DVarNum(ii_ds)*DsNodesNum;
BFEVarNumFluid = DfFEVarNum;
BFEVarNumSolid = DsFEVarNum;

% create zero sparse matrix for the 
% fluid-solid boundary condition matrix
 BMatrixDfs = sparse(BFEVarNumFluid,BFEVarNumSolid);

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
% Process all of the edges, which belong to the boundary
% Compute the interface conditions' matrices for the particular element
%===============================================================================

for ii_ed = 1:NBElements
    
    % find the start and the end nodes of the edge
    % depending on whether the fluid domain is the first or the last,
    % the orientation of the edge changes.
    % In order to identify the proper mid-side nodes, the order of the edge
    % nodes should be adjusted
    EdgeNodes = [];
    switch CompStruct.Model.DomainType{ii_d1}
        case 'fluid'
            EdgeNodes(1) = BElements(1,ii_ed); % the first node of the boundary element (edge)
            EdgeNodes(2) = BElements(2,ii_ed); % the second node of the boundary element (edge)
        case 'HTTI'
            % the orientation of the edge is reversed, together with the
            % direction of the normal
            EdgeNodes(1) = BElements(2,ii_ed); % the first node of the boundary element (edge)
            EdgeNodes(2) = BElements(1,ii_ed); % the second node of the boundary element (edge)
        otherwise
    end;
    
    % find the elements, which contain these nodes
    FNodes = ismember(FEMatrices.DElements{ii_df}(1:10,:), EdgeNodes(1:2));
    NNodes = sum(FNodes);
    [~,DfEdgeEl] = max(NNodes);
    FNodes = ismember(FEMatrices.DElements{ii_ds}(1:10,:), EdgeNodes(1:2));
    NNodes = sum(FNodes);
    [~,DsEdgeEl] = max(NNodes);
    
    % read the nodes, which belong to the adjoint (to the edge) elements
    TriNodes_f = FEMatrices.DElements{ii_df}(1:10,DfEdgeEl);
    TriNodes_s = FEMatrices.DElements{ii_ds}(1:10,DsEdgeEl);
    
%     % read the properties of the element (a, b, c, delta as defined in 
%     % the book of Zienkiewicz (v.1 p.181))
%     TriProps.a = DEMeshProps.a(:,D1EdgeEl);
%     TriProps.b = DEMeshProps.b(:,D1EdgeEl);
%     TriProps.c = DEMeshProps.c(:,D1EdgeEl);
%     TriProps.delta = DEMeshProps.delta(:,D1EdgeEl);    
    
    % find the nodes in between (mid-side nodes for cubic interpolation)
    EdgeNodesPos.Df(1) = find(ismember(TriNodes_f,EdgeNodes(1)));
    EdgeNodesPos.Df(2) = find(ismember(TriNodes_f,EdgeNodes(2)));
    % find the interior nodes of the edge element
    if ( EdgeNodesPos.Df(2) > mod(EdgeNodesPos.Df(1),3) )
        % find the third and the fourth nodes
        EdgeNodesPos.Df(3) = ( EdgeNodesPos.Df(1) + 1 )*2;
        EdgeNodesPos.Df(4) = ( EdgeNodesPos.Df(1) + 1 )*2 + 1;
    else
        % find the third and the fourth nodes
        EdgeNodesPos.Df(3) = ( EdgeNodesPos.Df(2) + 1 )*2 + 1;
        EdgeNodesPos.Df(4) = ( EdgeNodesPos.Df(2) + 1 )*2;
    end;
    EdgeNodes(3) = TriNodes_f(EdgeNodesPos.Df(3));
    EdgeNodes(4) = TriNodes_f(EdgeNodesPos.Df(4));
    
    EdgeNodesPos.Ds(1) = find(ismember(TriNodes_s,EdgeNodes(1)));
    EdgeNodesPos.Ds(2) = find(ismember(TriNodes_s,EdgeNodes(2)));
    EdgeNodesPos.Ds(3) = find(ismember(TriNodes_s,EdgeNodes(3)));
    EdgeNodesPos.Ds(4) = find(ismember(TriNodes_s,EdgeNodes(4)));

    BNodesFull = [BNodesFull EdgeNodes];
        
    % find the coordinates of the boundary element (edge)
    %  x coord of the first and the second node of the boundary element
    %  (edge)
    EdgeNodesX = squeeze(FEMatrices.MeshNodes(1,EdgeNodes(1:2)));  
    %  y coord of the first and the second node of the boundary element
    %  (edge)
    EdgeNodesY = squeeze(FEMatrices.MeshNodes(2,EdgeNodes(1:2)));  
    
%     % read x and coordinates of the corner nodes in the element
%     el_x_nodes = MeshNodes(1,MeshTri(1:3,ii_el));
%     el_y_nodes = MeshNodes(2,MeshTri(1:3,ii_el));  
%     el_x_nodes_full = MeshNodes(1,MeshTri(1:10,ii_el));
%     el_y_nodes_full = MeshNodes(2,MeshTri(1:10,ii_el));  

    % read physical properties of the element and prepare the matrices 
    % for the expansion of these elements into interpolating functions L_i
    % Df - fluid
    ElPhysProps.DfEl = CompStruct.Methods.getPhysProps{ii_df}...
        (CompStruct,BasicMatrices,FEMatrices,ii_df,DfEdgeEl);
    % Ds - HTTI
    ElPhysProps.DsEl = CompStruct.Methods.getPhysProps{ii_ds}...
        (CompStruct,BasicMatrices,FEMatrices,ii_ds,DsEdgeEl);

    % find the components of the edge vector
    EdgeProps.DxEdge = EdgeNodesX(2) - EdgeNodesX(1); %dx=(x(2)-x(1));
    EdgeProps.DyEdge = EdgeNodesY(2) - EdgeNodesY(1); %dy=(y(2)-y(1));
    EdgeProps.Dl = sqrt( EdgeProps.DxEdge^2 + EdgeProps.DyEdge^2 ); %length of the boundary element 
    
    % define the oriented normal to the boundary edge by taking the cross product
    % of the edge vector with the (-e_z) basis vector
    EdgeProps.Normal = transpose(cross([EdgeProps.DxEdge/EdgeProps.Dl,...
        EdgeProps.DyEdge/EdgeProps.Dl,0],[0,0,-1]));   
    % cross product of the boundary edge vector and z axis vector ( required to define the normal to the edge)
    
    % compute the contributions of the element to the stiffness and mass
    % matrices
    [BMatrixDfs_edge, BMatrixDsf_edge]...
            = CompStruct.Methods.IC_matrix{ii_int}(BasicMatrices, CompStruct, ii_int, ii_df, ii_ds, ...
                             ElPhysProps, EdgeProps, EdgeNodesPos);
                         

    % identify positions of the element nodes and the variables in the
    % respective matrices
    for ii_n = 1:4
        DfNodePos = find(ismember(FEMatrices.DNodes{ii_df},EdgeNodes(ii_n)));
        DsNodePos = find(ismember(FEMatrices.DNodes{ii_ds},EdgeNodes(ii_n)));
        VarRowArr( ( CompStruct.Data.DVarNum(ii_df)*( ii_n - 1 ) + 1 ):( CompStruct.Data.DVarNum(ii_df)*ii_n ) ) = ...
            [( CompStruct.Data.DVarNum(ii_df)*( DfNodePos - 1 ) + 1 ):( CompStruct.Data.DVarNum(ii_df)*DfNodePos )];
        VarColArr( ( CompStruct.Data.DVarNum(ii_ds)*( ii_n - 1 ) + 1 ):( CompStruct.Data.DVarNum(ii_ds)*ii_n ) ) = ...
            [( CompStruct.Data.DVarNum(ii_ds)*( DsNodePos - 1 ) + 1 ):( CompStruct.Data.DVarNum(ii_ds)*DsNodePos )];
    end;

    % prepare sparse matrices to insert into the matrices for the domain
    % compute rows and columns, to which to insert
    VarRows = repmat(VarRowArr,1,CompStruct.Data.DVarNum(ii_ds)*4);
    VarCols = reshape(repmat(VarColArr,CompStruct.Data.DVarNum(ii_df)*4,1),1,[]);
    % create vectors of elements, which to insert
    strBMatrixDfs_el = reshape(BMatrixDfs_edge,1,[]);
    % create sparse matrices
    spBMatrixDfs_el = sparse(VarRows,VarCols,strBMatrixDfs_el,BFEVarNumFluid,BFEVarNumSolid);
                         
    % insert the matrices for the element into the matrices for the domain
    BMatrixDfs = BMatrixDfs + spBMatrixDfs_el;

%     % calculation of boundary full matrix in fluid
%     Fluid_bound_matrix(rows,colons)=Fluid_bound_matrix(rows,colons)+Int_fluid;

end;

%===============================================================================
% Assemble the full list of the nodes, which belong to the interface
%===============================================================================

[FEMatrices.BNodesFull{ii_int}] = sort(unique(BNodesFull));

%===============================================================================
% Compute the solid-fluid boundary condition matrix - it is transpose of
% the fluid-solid one according to the formulation (see the papers and the reports)
%===============================================================================

BMatrixDsf = transpose(BMatrixDfs);
% K1Matrix_d = sparse(K1Matrix_d);

% assign the rectangular blocks of the matrices responsible for the
% boundary conditions according to which domain is the first - solid or
% fluid 
% this is necessary to account for the proper size of the matrix blocks and
% for the ease of their proper placement into the full matrix
switch CompStruct.Model.DomainType{ii_d1} 
    case 'fluid'
        FEMatrices.PMatrixD12{ii_int} = BMatrixDfs; %++ or -- is Ok
        FEMatrices.PMatrixD21{ii_int} = BMatrixDsf;
    case 'HTTI'
        % change the sign because the normal vector should be directed from
        % the fluid into the solid
        FEMatrices.PMatrixD12{ii_int} = - BMatrixDsf;
        FEMatrices.PMatrixD21{ii_int} = - BMatrixDfs;
    otherwise
end;

end

