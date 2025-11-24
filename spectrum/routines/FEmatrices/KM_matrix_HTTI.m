function [ElMatrices] = KM_matrix_HTTI(BasicMatrices, CompStruct,  FEMatrices, ...
                                       ii_d, ii_el, ElPhysProps, TriProps)
     
%   Part of the toolbox for solving problems of wave propagation
% in arbitrary anisotropic inhomogeneous waveguides.
% For details see User manual 
% and comments to the main script gen_aniso.m
%
%   KM_el_matrix_HTTI function computes the matrices, which indicates the expansion
%  and the coefficients of the expansion of the integrands for K_i and M matrices' elements
%  into the intepolating functions L_j (L1, L2, L3).
%  For HTTI solid.
%  The K_i and M matrices are defined according to Bartoli_Marzani_DiScalea_Viola_JSoundVibration_v295_p685_2006.
%  N_m = (N1*E_N, ..., Nn*E_N) = kron(N_vector,E_N) (3 x 3Nn matrix)
%  B2_m = Lz_m * N_m = Lz_m * kron(N_vector,E_N) (6 x 3Nn matrix)
%  B1_m = Lx_m*d/dx N_m + Ly_m*d/dy N_m = ((Lx_m*d/dx + Ly_m*d/dy)*N1, ..., (Lx_m*d/dx + Ly_m*d/dy)*Nn)
%
%  K1_m = B1_adjoint * C_m * B1 
%  K1_m ab :::  = conv(conj(B1_m ia :::), conv(C_ij :::, B1_m jb :::) ) (3Nn x 3Nn matrix)
%  K2_m = B1_adjoint * C_m * B2 - B2_adjoint * C_m * B1
%  K2_m ab :::  = conv(conj(B1_m ia :::), conv(C_ij :::, B2_m jb :::) ) 
%                   - conv(conj(B2_m ia :::), conv(C_ij :::, B1_m jb :::) ) (3Nn x 3Nn matrix)
%  K3_m = B2_adjoint * C_m * B2 
%  K3_m ab :::  = conv(conj(B2_m ia :::), conv(C_ij :::, B2_m jb :::) ) (3Nn x 3Nn matrix)
%  M_m = N_adjoint * Rho_m * N 
%  M_m ab :::  = conv(conj(N_m ia :::), conv(C_ij :::, N_m jb :::) ) (3Nn x 3Nn matrix)
%
%  Current implementation is for cubic elements.
% NB!!! This implementation is specific for FE computations
%
%   [T.Zharnikov, SMR v0.3_08.2014]
%
% function [B1tCB1Matrix, B1tCB2Matrix, B2tCB1Matrix, B2tCB2Matrix, NtRhoNMatrix]...
%            = KM_el_matrix_HTTI(CMatrix, RhoMatrix, BasicMatrices, ...
%                            NLMatrix, dxNLMatrix, dyNLMatrix)
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
%   Code for KM_el_matrix_HTTI
%
%###############################################################################
%===============================================================================
% Initialization
%===============================================================================

% Identify number of nodes 

% E_matrix = eye(3);
N_nodes = CompStruct.Advanced.N_nodes;

% get the triangle area 
delta = TriProps.delta(1);

% Allocate memory for the matrices.
% The implementation is for cubic elements, hence expansion is up to the
% 3rd power in L_i

ElMatrices.MMatrix = zeros(CompStruct.Data.DVarNum(ii_d)*N_nodes,CompStruct.Data.DVarNum(ii_d)*N_nodes); %multiply by delta afterwards
ElMatrices.K1Matrix = zeros(CompStruct.Data.DVarNum(ii_d)*N_nodes,CompStruct.Data.DVarNum(ii_d)*N_nodes); %multiply by 1/delta afterwards
ElMatrices.K2Matrix = zeros(CompStruct.Data.DVarNum(ii_d)*N_nodes,CompStruct.Data.DVarNum(ii_d)*N_nodes); %multiply by 1 afterwards
ElMatrices.K3Matrix = zeros(CompStruct.Data.DVarNum(ii_d)*N_nodes,CompStruct.Data.DVarNum(ii_d)*N_nodes); %multiply by delta afterwards

%===============================================================================
% Compute the elements of K_i and M matrices according to the above
% presented formulas
%===============================================================================

   [ElMatrices] = CompStruct.Methods.KM_el_matrix{ii_d}(BasicMatrices, CompStruct,  FEMatrices, ...,
                                                            ii_d, ii_el, ElPhysProps, TriProps);

% mass matrix term - the kinetic energy
ElMatrices.MMatrix = - delta*ElMatrices.NtRhoNMatrix;

% stiffness matrix terms - the potential energy
ElMatrices.K1Matrix = - (1/delta)*ElMatrices.B1tCB1Matrix;
ElMatrices.K2Matrix = - ( ElMatrices.B1tCB2Matrix - ElMatrices.B2tCB1Matrix );
ElMatrices.K3Matrix = - delta*ElMatrices.B2tCB2Matrix;

end
