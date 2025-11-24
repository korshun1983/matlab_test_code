function [CompStruct,BasicMatrices,FEMatrices,FullMatrices] = St3_1_PrepareBasicMatrices_sp_SAFE(CompStruct,InputParam)
      
%   Part of the toolbox for solving problems of wave propagation
% in arbitrary anisotropic inhomogeneous waveguides.
% For details see User manual 
% and comments to the main script gen_aniso.m
%
%   St3_1_PrepareBasicMatrices_sp_SAFE.m M-file      
%      St3_1_PrepareBasicMatrices_sp_SAFE.m 
% prepares the basic matrices (blocks), which will be used 
% to construct spectral method operators during computations.
% Such matrices are radial, frequency grids, 
% positions of blocks corresponding to various layers and harmonics
% in the big spectral method matrix,
% C(n), A(n)
%
% NB!!! This implementation is specific for spectrum computation 
% with the spectral method.
%
%   [T.Zharnikov, D.Syresin, SMR v0.3_08.2014]
%
% function [CompStruct,BasicMatrices] =
% St3_1_PrepareBasicMatrices_sp_SAFE(CompStruct)
%
%  Inputs - 
%
%       CompStuct - structure containing the parameters of the model 
%                   in the form, which is convenient and ready for
%                   computations
%
%Outputs -
%
%       CompStuct - structure containing the parameters of the model 
%                   in the form, which is convenient and ready for
%                   computations. It is updated with the additional
%                   infromation.
%
%       BasicMatrices - structure containing the basic matrices (differentiation, 
%                   matrix representation blocks, etc.);
%
%  M-files required-
%
% Last Modified by Timur Zharnikov SMR v0.3_08.2014

%###############################################################################
%
%   Code for St3_1_PrepareBasicMatrices_sp_SAFE
%
%###############################################################################
%===============================================================================
% Initialization
%===============================================================================

% Retaining all necessary information 
% CompStruct = InputParam;

%===============================================================================
% Construction of various grids, matrices and 
% position arrays
%===============================================================================

%===============================================================================
% Introduce basic matrices, like L_x, L_y, L_z
%===============================================================================

[BasicMatrices] = CompStruct.Methods.AssembleBasicMatrices(CompStruct);
    
%==============================================================================
% f - grid, Rx, Ry
%==============================================================================
BasicMatrices.f_grid=CompStruct.f_grid;
num_layers=length(InputParam.Model.DomainType); nl=num_layers;
var_vel=CompStruct.Asymp.V_SH*1e3/(CompStruct.f_grid(CompStruct.if_grid)*CompStruct.Misc.F_conv);
if strcmpi(CompStruct.Model.LDomain_in_LSH,'yes') && strcmpi(CompStruct.Model.AddDomain_Exist,'yes')
   var_wl=var_vel;
   CompStruct.Model.AddDomainL_m=CompStruct.Model.AddDomainL*var_wl; % in meters
   if strcmpi(CompStruct.Model.AddDomainLoc,'ext') 
       CompStruct.Model.DomainRx(nl-1)=InputParam.Model.DomainRx(nl-2)+InputParam.Model.DomainRx(nl-1)*var_wl;
       CompStruct.Model.DomainRy(nl-1)=InputParam.Model.DomainRy(nl-2)+InputParam.Model.DomainRy(nl-1)*var_wl;
       CompStruct.Model.DomainRx(nl)=CompStruct.Model.DomainRx(nl-1)+CompStruct.Model.AddDomainL_m;
       CompStruct.Model.DomainRy(nl)=CompStruct.Model.DomainRy(nl-1)+CompStruct.Model.AddDomainL_m;

   end
   if strcmpi(CompStruct.Model.AddDomainLoc,'int') 
       CompStruct.Model.DomainRx(nl)=InputParam.Model.DomainRx(nl-1)+InputParam.Model.DomainRx(nl)*var_wl;
       CompStruct.Model.DomainRy(nl)=InputParam.Model.DomainRy(nl-1)+InputParam.Model.DomainRy(nl)*var_wl;
   end
end

if strcmpi(CompStruct.Model.LDomain_in_LSH,'none') && strcmpi(CompStruct.Model.AddDomain_Exist,'yes')
    CompStruct.Model.AddDomainL_m=CompStruct.Model.AddDomainL; % in meters
    if strcmpi(CompStruct.Model.AddDomainLoc,'ext') 
       CompStruct.Model.DomainRx(nl)=CompStruct.Model.DomainRx(nl-1)+CompStruct.Model.AddDomainL_m;
       CompStruct.Model.DomainRy(nl)=CompStruct.Model.DomainRy(nl-1)+CompStruct.Model.AddDomainL_m;
   end
   if strcmpi(CompStruct.Model.AddDomainLoc,'int')
      varx=CompStruct.Model.DomainRx(nl)-CompStruct.Model.AddDomainL_m;
      vary=CompStruct.Model.DomainRy(nl)-CompStruct.Model.AddDomainL_m;
      if varx<=CompStruct.Model.DomainRx(nl-1) || vary<=CompStruct.Model.DomainRy(nl-1)
         fprintf(1,'Error in specifying model geometry!!!\n\n'); return;
      end
   end
end

if strcmpi(CompStruct.Model.LDomain_in_LSH,'yes') && strcmpi(CompStruct.Model.AddDomain_Exist,'none')
   var_wl=var_vel;
   CompStruct.Model.DomainRx(nl)=InputParam.Model.DomainRx(nl-1)+InputParam.Model.DomainRx(nl)*var_wl;
   CompStruct.Model.DomainRy(nl)=InputParam.Model.DomainRy(nl-1)+InputParam.Model.DomainRy(nl)*var_wl;
end

if strcmpi(CompStruct.Model.LDomain_in_LSH,'none') && strcmpi(CompStruct.Model.AddDomain_Exist,'none')
   ;
end

%==============================================================================
% Assembling the mesh
%===============================================================================
[FEMatrices.MeshNodes,FEMatrices.BoundaryEdges, FEMatrices.MeshTri,...
    FEMatrices.MeshProps,CompStruct] = CompStruct.Methods.PrepareMesh(CompStruct);
    
%===============================================================================
% Preparing the matrices, which are necessary for the description of each
% domain
%===============================================================================
    
for ii_d = 1:CompStruct.Data.N_domain
    
%===============================================================================
% Preparing matrix, which characterizes variation of physical properties of the layers.
%   E.g. construction of C(n)(r) for TTI layers
%===============================================================================

    [FEMatrices.PhysProp{ii_d}] = CompStruct.Methods.PreparePhysProp{ii_d}(CompStruct,ii_d);

%===============================================================================
% Finding the elements and the nodes, which belong to the particular domains
%===============================================================================

    % finding elements
    [FEMatrices.DElements{ii_d}] = ...
        FEMatrices.MeshTri(:,ismember(FEMatrices.MeshTri(11,:),ii_d)); 
    % finding their properties
    [FEMatrices.DEMeshProps{ii_d}.DS] = ...
        FEMatrices.MeshProps.DS(:,ismember(FEMatrices.MeshTri(11,:),ii_d));
    [FEMatrices.DEMeshProps{ii_d}.delta] = ...
        FEMatrices.MeshProps.delta(:,ismember(FEMatrices.MeshTri(11,:),ii_d));
    [FEMatrices.DEMeshProps{ii_d}.a] = ...
        FEMatrices.MeshProps.a(:,ismember(FEMatrices.MeshTri(11,:),ii_d));
    [FEMatrices.DEMeshProps{ii_d}.b] = ...
        FEMatrices.MeshProps.b(:,ismember(FEMatrices.MeshTri(11,:),ii_d));
    [FEMatrices.DEMeshProps{ii_d}.c] = ...
        FEMatrices.MeshProps.c(:,ismember(FEMatrices.MeshTri(11,:),ii_d));
    % finding the nodes
    % define the array with the numbers of nodes, belonging to the domain
    DNodesEl = FEMatrices.DElements{ii_d}(1:10,:);
    % transform it into the vector, sort, and find the nodes, belonging to
    % the domain
    DNodesVec = reshape(DNodesEl,1,[]);
    % find the indices (ordered) of the nodes, which belong to the domain
    [FEMatrices.DNodes{ii_d}] = sort(unique(DNodesVec));
    % OrderDNodesVec = sort(DNodesVec);
    % DiffDNodesVec = diff(OrderDNodesVec);
    % DNodesIndex = [1, ( 1 + find(DiffDNodesVec) )];
    % % find the indices (ordered) of the nodes, which belong to the domain
    % DNodes = OrderDNodesVec(NodesIndex);

    [FEMatrices.DNodesRem{ii_d}] = [];
    [FEMatrices.DNodesComp{ii_d}] = [];
    [FEMatrices.DTakeFromVarPos{ii_d}] = [];
    [FEMatrices.DPutToVarPos{ii_d}] = [];
    [FEMatrices.DZeroVarPos{ii_d}] = [];

    
%===============================================================================
% Computing the blocks of mass and stiffness matrices for FE computations 
%===============================================================================
%ii_d;
%tic;
%input('Do you want more 1? Y/N [Y]: ', 's');
    % computing the blocks of the mass and stiffness matrices for each domain
    [FEMatrices] = CompStruct.Methods.MatricesParts_sp_SAFE{ii_d} ...
                    (CompStruct,BasicMatrices,FEMatrices,ii_d); 
%toc;

%cutting digits to speed up the computations - suggested by Denis Syresin
%
% power=fix(log10(max(max(abs(M)))));
% M=round(M/10^(power-15))*10^(power-15);

end;

%===============================================================================
% Preparing the matrices, which are necessary for the description of each boundary
% and inserting them into the boundary conditions matrices
%===============================================================================

for ii_int = 1:( CompStruct.Data.N_domain - 1 )

%===============================================================================
% Finding the elements and the nodes, which belong to the particular domains
%===============================================================================

    % finding the nodes
    % define the array with the numbers of nodes, belonging to the boundary
    BNodesEl = FEMatrices.BoundaryEdges(1:2,ismember(FEMatrices.BoundaryEdges(3,:),ii_int));
    % transform it into the vector, sort, and find the nodes, belonging to
    % the boundary
    BNodesVec = reshape(BNodesEl,1,[]);
    % find the indices (ordered) of the nodes, which belong to the domain
    [FEMatrices.BNodes{ii_int}] = sort(unique(BNodesVec)); 
    ii_d1 = ii_int;
    ii_d2 = ii_int + 1;

%===============================================================================
% Computing the blocks of the matrices, which describe the interface and
% bounary conditions ( for FE computations )
%===============================================================================

    % computing the blocks of the mass and stiffness matrices for each
    % domain (D1 - ii_l, D2 - ii_l)
    [FEMatrices] = ...
        CompStruct.Methods.IC_Matrices_sp_SAFE{ii_int}(CompStruct,BasicMatrices,FEMatrices,...
                ii_int,ii_d1,ii_d2); %,integration_acc);

end;

%===============================================================================
% Calculating indices, which describe start positions of various blocks of
% matrices
%===============================================================================

BasicMatrices = CompStruct.Methods.FindPos(CompStruct,FEMatrices,BasicMatrices);

%===============================================================================
% Construct the full matrices
%===============================================================================

% start assembling full matrices - insert the blocks for the first domain
FullMatrices.K1Matrix = FEMatrices.K1Matrix_d{1};
FullMatrices.K2Matrix = FEMatrices.K2Matrix_d{1};
FullMatrices.K3Matrix = FEMatrices.K3Matrix_d{1};
FullMatrices.MMatrix = FEMatrices.MMatrix_d{1};
FullMatrices.PMatrix = FEMatrices.PMatrix_d{1};
FEMatrices.DNodesRem{1} = [];
FEMatrices.DNodesComp{1} = FEMatrices.DNodes{1};

for ii_int = 1:( CompStruct.Data.N_domain - 1 )

%===============================================================================
% Assembling the blocks of the matrices together and inserting the boundary
% conditions at the domain interfaces
%===============================================================================

    ii_int;
    ii_d1 = ii_int;
    ii_d2 = ii_int + 1;
    % assembling the mass and stiffness matrices for domain by domain (D1 - ii_l, D2 - ii_l) from the precomputed blocks
    % besides, the structure of variables vectors and the positions of variables are adjusted as well
    [FEMatrices,FullMatrices] = ...
        CompStruct.Methods.AssembleFullMatrices{ii_int}(CompStruct,BasicMatrices,FEMatrices,FullMatrices,...
                ii_int,ii_d1,ii_d2); %,integration_acc);

end;

%===============================================================================
% Assembling the blocks of the matrices together and inserting the boundary
% conditions at the outer interface
%===============================================================================

ii_int = CompStruct.Data.N_domain;
ii_d1 = ii_int;
ii_d2 = ii_int + 1;

% define the array with the numbers of nodes, belonging to the outer boundary
BNodesEl = FEMatrices.BoundaryEdges(1:2,ismember(FEMatrices.BoundaryEdges(3,:),ii_int));
% transform it into the vector, sort, and find the nodes, belonging to
% the boundary
BNodesVec = reshape(BNodesEl,1,[]);
% find the indices (ordered) of the nodes, which belong to the domain
[FEMatrices.BNodes{ii_int}] = sort(unique(BNodesVec));

% assembling the mass and stiffness matrices for the outer boundary (D1 - ii_l, D2 - ii_l) from the precomputed blocks
% besides, the structure of variables vectors and the positions of variables are adjusted as well
[FEMatrices,FullMatrices] = ...
    CompStruct.Methods.AssembleFullMatrices{ii_int}(CompStruct,BasicMatrices,FEMatrices,FullMatrices,...
    ii_int,ii_d1,ii_d2); %,integration_acc);

%===============================================================================
% Remove the rows and columns of the full matrices corresponding to the nodes, 
% which should be removed and do not participate in the computation 
% (e.g. rigid boundary condition, FF or SS interface, etc.)
%===============================================================================

[FEMatrices,FullMatrices] = ...
    CompStruct.Methods.RemoveRedundantVariables(CompStruct,BasicMatrices,FEMatrices,FullMatrices);

end
