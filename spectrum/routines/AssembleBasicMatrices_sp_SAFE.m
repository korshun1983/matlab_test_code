function [BasicMatrices] = ...
                    AssembleBasicMatrices_sp_SAFE(CompStruct)

%   Part of the toolbox for solving problems of wave propagation
% in arbitrary anisotropic inhomogeneous waveguides.
% For details see User manual 
% and comments to the main script gen_aniso.m
%
% AssembleBasicMatrices_sp_SAFE M-file      
%      AssembleBasicMatrices_sp_SAFE, by itself, prepares 
%      such matrices as Lx, Ly, Lz,
%      identity matrices, etc. 
%
%   [T.Zharnikov, SMR v0.3_08.2014]
%
% function [BasicMatrices] = AssembleBasicMatrices_sp_SAFE()
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
%   Code for AssembleBasicMatrices_sp_SAFE
%
%###############################################################################

%===============================================================================
% Preparing the matrices
% For reference see Auld book or
% Bartoli_Marzani_DiScalea_Viola_JSoundVibration_v295_p685_2006
%===============================================================================    

BasicMatrices.Lx = [ 1, 0, 0;...
                     0, 0, 0;...
                     0, 0, 0;...
                     0, 0, 0;...
                     0, 0, 1;...
                     0, 1, 0];
   
BasicMatrices.Ly = [ 0, 0, 0;...
                     0, 1, 0;...
                     0, 0, 0;...
                     0, 0, 1;...
                     0, 0, 0;...
                     1, 0, 0];
   
BasicMatrices.Lz = [ 0, 0, 0;...
                     0, 0, 0;...
                     0, 0, 1;...
                     0, 1, 0;...
                     1, 0, 0;...
                     0, 0, 0];

BasicMatrices.Lx_fluid = [ 1; 0; 0];
   
BasicMatrices.Ly_fluid = [ 0; 1; 0];
   
BasicMatrices.Lz_fluid = [ 0; 0; 1];

BasicMatrices.E3 = eye(3);

% compute matrices for the integrals of the L_i functions and their monomials
BasicMatrices.LEdgeIntMatrix9 = CompStruct.Methods.L1L2_int_matrix(CompStruct.Advanced.N_nodes);
BasicMatrices.LIntMatrix9 = CompStruct.Methods.L1L2L3_int_matrix(CompStruct.Advanced.N_nodes);

% N.B.!!! The matrices NL, dNL, and in particular, NNNConv, dNNNConv,
% NNdNConv, dNNdNConv, are sparse (from 0.01% to 2%) and require a lot of
% storage space (up to 100Mb) in full format.
% If necessary, the storage space can be significantly optimized by using
% and converting to the
% sparse storage
% sizeNL = sizeNL;
% sparseNL = sparse(reshape(NL,round(sqrt(prod(size(NL)))),[]));
% fullNL = reshape(sparseNL,sizeNL);

% compute the expansion of the shape functions N_i and their derivatives 
% with respect to the L_i functions (interpolating polynomials)
[BasicMatrices.NLMatrix, BasicMatrices.NodeLCoord ] =  CompStruct.Methods.NL_matrix();
[BasicMatrices.dNLMatrix, BasicMatrices.dNLMatrix_val] =  CompStruct.Methods.dNL_matrices(BasicMatrices);

% compute the expansion of the double convolutions 
% of the shape functions N_i and their derivatives 
% with respect to the L_i functions (interpolating polynomials)
[BasicMatrices] = CompStruct.Methods.ConvolveMatrices(CompStruct,BasicMatrices);

% compute the expansion of the edge shape functions N_i 
% with respect to the L_i functions (interpolating polynomials)
BasicMatrices.NLEdgeMatrix =  CompStruct.Methods.NLEdge_matrix();

% compute the expansion of the double convolutions 
% of the shape functions N_i and their derivatives 
% with respect to the L_i functions (interpolating polynomials)
[BasicMatrices] = CompStruct.Methods.ConvolveEdgeMatrices(CompStruct,BasicMatrices);

end