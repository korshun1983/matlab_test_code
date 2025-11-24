function [FullMatrices] = ...
                    AssembleFullFEMatrices_sp_SAFE()

% Currently this procedure is not used
% It contains several possible leads to modify other procedures
                
%   Part of the toolbox for solving problems of wave propagation
% in arbitrary anisotropic inhomogeneous waveguides.
% For details see User manual 
% and comments to the main script gen_aniso.m
%
% AssembleFullFEMatrices_sp_SAFE M-file      
%      AssembleFullFEMatrices_sp_SAFE, by itself, prepares 
%      full matrices (mass, stiffness, etc.) for SAFE computations
%
%   [T.Zharnikov, SMR v0.3_08.2014]
%
% function [FullMatrices] = ...
%                    AssembleFullFEMatrices_sp_SAFE()
%
%  Inputs - 
%
%       CompStuct - structure containing the parameters of the model;
%
%       ii_l - number of the layer, for which the basic matrices are
%              assembled
%
%  Outputs -
%
%       DDR1 - first derivative differentiation matrix 
%
%       DDR2 - second derivative differentiation matrix 
%
%       ZM - zero matrix
%
%       IdM - identity matrix
%
%       DiagR1 - diagonal matrix, with the elements on the diagonal,
%                  representing r^(-1) vector on the radial grid
%
%       DiagR2 - diagonal matrix, with the elements on the diagonal,
%                  representing r^(-2) vector on the radial grid
%
%       RMatrix - radial grid vector
%
%  M-files required-
%
% Last Modified by Timur Zharnikov SMR v0.3_08.2014

%###############################################################################
%
%   Code for AssembleFullFEMatrices_sp_SAFE
%
%###############################################################################

FullMatrices.K1Matrix = sparse(CompStruct.Data.DVarNum{N_domain}(2),CompStruct.Data.DVarNum{N_domain}(2));
FullMatrices.K2Matrix = sparse(CompStruct.Data.DVarNum{N_domain}(2),CompStruct.Data.DVarNum{N_domain}(2));
FullMatrices.K3Matrix = sparse(CompStruct.Data.DVarNum{N_domain}(2),CompStruct.Data.DVarNum{N_domain}(2));
FullMatrices.MMatrix = sparse(CompStruct.Data.DVarNum{N_domain}(2),CompStruct.Data.DVarNum{N_domain}(2));
BCMatrices.ICBCMatrix = sparse(CompStruct.Data.DVarNum{N_domain}(2),CompStruct.Data.DVarNum{N_domain}(2));

for ii_d = 1:CompStruct.Data.N_domains
    
%===============================================================================
% Inserting matrix blocks into the full matrices
%===============================================================================

    [ DVecRepRow, DVecRepCol, DVecRepV ] = find(K1Matrix_d{ii_d});
    DVecRepRow = DVecRepRow + CompStruct.Data.DVarNum{ii_d}(1) - 1;
    DVecRepCol = DVecRepCol + CompStruct.Data.DVarNum{ii_d}(1) - 1;
    DSize = size(K1Matrix_d{ii_d});
    K1Matrix = K1Matrix + sparse(DVecRepRow, DVecRepCol, DVecRepV, DSize(1), DSize(2));
    [ DVecRepRow, DVecRepCol, DVecRepV ] = find(K2Matrix_d{ii_d});
    DVecRepRow = DVecRepRow + CompStruct.Data.DVarNum{ii_d}(1) - 1;
    DVecRepCol = DVecRepCol + CompStruct.Data.DVarNum{ii_d}(1) - 1;
    DSize = size(K2Matrix_d{ii_d});
    K2Matrix = K2Matrix + sparse(DVecRepRow, DVecRepCol, DVecRepV, DSize(1), DSize(2));
    [ DVecRepRow, DVecRepCol, DVecRepV ] = find(K3Matrix_d{ii_d});
    DVecRepRow = DVecRepRow + CompStruct.Data.DVarNum{ii_d}(1) - 1;
    DVecRepCol = DVecRepCol + CompStruct.Data.DVarNum{ii_d}(1) - 1;
    DSize = size(K3Matrix_d{ii_d});
    K3Matrix = K3Matrix + sparse(DVecRepRow, DVecRepCol, DVecRepV, DSize(1), DSize(2));
    [ DVecRepRow, DVecRepCol, DVecRepV ] = find(MMatrix_d{ii_d});
    DVecRepRow = DVecRepRow + CompStruct.Data.DVarNum{ii_d}(1) - 1;
    DVecRepCol = DVecRepCol + CompStruct.Data.DVarNum{ii_d}(1) - 1;
    DSize = size(MMatrix_d{ii_d});
    MMatrix = MMatrix + sparse(DVecRepRow, DVecRepCol, DVecRepV, DSize(1), DSize(2));

%cutting digits to speed up the computations - suggested by Denis Syresin
%
% power=fix(log10(max(max(abs(M)))));
% M=round(M/10^(power-15))*10^(power-15);

end;

% symmetrization of the matrices by going to u_z -> i*u_z as described in
% e.g. Bartoli_Marzani_DiScalea_Viola_JSoundVibration_v295_p685_2006
% Makes sense for the lossless waveguides, as the matrices become symmetric.
% Does not work for lossy waveguides.
% Suggested by Denis Syresin, 2013. Not used now
% [M_full,M2_full,K1_full,K2_full,K3_full,M_Source] = symmetrization_new( M_full,M2_full,K1_full,K2_full,K3_full,M_Source,nodes);


% solid external boundary condition and symmetrization procedure
nodes=[1:lsm sort(repmat(1:lsm,1,3))];
[ M_full ] = Rigid_boundary( p,e,M_full,rigid_solid_boundaries);
% reduce not intersecting nodes (all zero rows in volume matrice)
[M_full,M2_full,K1_full,K2_full,K3_full,M_Source,nodes] = reduce_matrix(M_full,M2_full,K1_full,K2_full,K3_full,M_Source,lsm);

%===============================================================================
% Inserting interface and boundary conditions matrices blocks into the full matrices
%===============================================================================

[FullMatrices] = CompStruct.Methods.InsertBC(FullMatrices,BasicMatrices,CompStruct);

for ii_int = 1:(CompStruct.Data.N_domains - 1)
    ii_D1 = ii_int;
    ii_D2 = ii_D1 + 1;
    [FullMatrices] = CompStruct.Methods.InsertIC(...
        FullMatrices,BasicMatrices,CompStruct,ii_D1,ii_D2,ii_int);
    
%cutting digits to speed up the computations - suggested by Denis Syresin
%
% power=fix(log10(max(max(abs(M)))));
% M=round(M/10^(power-15))*10^(power-15);

end;

end