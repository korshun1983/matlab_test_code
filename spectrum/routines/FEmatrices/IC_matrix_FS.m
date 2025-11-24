function [BMatrixDfs_edge, BMatrixDsf_edge]...
            = IC_matrix_FS(BasicMatrices, CompStruct, ii_int, ii_df, ii_ds, ...
                             ElPhysProps, EdgeProps, EdgeNodesPos)
     
%   Part of the toolbox for solving problems of wave propagation
% in arbitrary anisotropic inhomogeneous waveguides.
% For details see User manual 
% and comments to the main script gen_aniso.m
%
%   IC_matrix_FS function computes the interface conditions' matrices, 
%  which enter into the FE formulation.
%  They are obtained by multiplying the matrix of the the coefficients of the expansion 
%  of the integrands for Nfl_t_rho_(nN) matrix elements into the intepolating functions 
%  L_j (L1, L2) defined on the edge onto the values of the integrals 
%  L*\int_{0}^{1} (1 - xi)^a xi^b.
%  For HTTI solid - fluid contact.
%  The K_i and M matrices are defined according to Bartoli_Marzani_DiScalea_Viola_JSoundVibration_v295_p685_2006.
%  N_m = (N1*E_N, ..., Nn*E_N) = kron(N_vector,E_N) (3 x 3Nn matrix)
%  (nN) = (N1*nx, N1*ny, 0, N2*nx, ..., Nn*nx, Nn*ny, 0) (1 x 3Nn matrix)
%  N_fl = (N1, ..., Nn) =  (Nn matrix)
%
%  Current implementation is for cubic elements.
% NB!!! This implementation is specific for FE computations
%
%   [T.Zharnikov, SMR v0.3_08.2014]
%
% function [BMatrixD1_edge]...
%            = IC_matrix_FS(BasicMatrices, ElPhysProps, EdgeProps)
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
%   Code for IC_matrix_FS
%
%###############################################################################
%===============================================================================
% Initialization
%===============================================================================

% Identify number of nodes 

% E_matrix = eye(3);
NEdge_nodes = CompStruct.Advanced.NEdge_nodes;

% get the length of the edge (the element of the boundary) 
Dl = EdgeProps.Dl;

% Allocate memory for the matrices.
% The implementation is for cubic elements, hence expansion is up to the
% 3rd power in L_i

BMatrixDfs_edge = zeros(CompStruct.Data.DVarNum(ii_df)*NEdge_nodes,CompStruct.Data.DVarNum(ii_ds)*NEdge_nodes); %multiply by Dl afterwards
BMatrixDsf_edge = zeros(CompStruct.Data.DVarNum(ii_ds)*NEdge_nodes,CompStruct.Data.DVarNum(ii_df)*NEdge_nodes); %multiply by Dl afterwards

%===============================================================================
% Compute the elements of K_i and M matrices according to the above
% presented formulas
%===============================================================================

[ElMatrices] = CompStruct.Methods.IC_el_matrix{ii_int}...
    (BasicMatrices, CompStruct, ElPhysProps, EdgeProps, EdgeNodesPos);

% process all of the elements of the 2d 
% NtRhoNMatrix, etc. matrices (3*N_nodes x 3*N_nodes)
BMatrixDfs_edge = Dl*ElMatrices.NfltRhonNMatrix;
BMatrixDsf_edge = Dl*ElMatrices.NfltRhonNMatrix';
        
end