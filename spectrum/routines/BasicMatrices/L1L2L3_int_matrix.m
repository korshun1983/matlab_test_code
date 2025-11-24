function [LIntMatrix] = L1L2L3_int_matrix(N_degree)
     
%   Part of the toolbox for solving problems of wave propagation
% in arbitrary anisotropic inhomogeneous waveguides.
% For details see User manual 
% and comments to the main script gen_aniso.m
%
%   L1L2L3_int_matrix function computes the integrals of the form
%  int int {over triangle} L1^{a}L2^{b}L3^{c} dx dy using the exact integration 
%  formula from the book of Zienkiewicz (8.38):
%  (a!b!c!/(a+b+c+2)!)*2*delta, where delta is the total area of the triangle
%  This function computes the matrix of coefficients 2*(a!b!c!/(a+b+c+2)!)
% NB!!! This implementation is specific for FE computations
%
%   [T.Zharnikov, SMR v0.3_08.2014]
%
% function [LIntMatrix] = L1L2L3_int_matrix()
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
%   Code for L1L2L3_int_matrix
%
%###############################################################################
%===============================================================================
% Initialization
%===============================================================================

% allocate memory for LIntMatrix (matrix of coefficients to compute 
% intregrals of the form int int {over triangle} L1^{a}L2^{b}L3^{c} )
LIntMatrix = zeros(N_degree,N_degree,N_degree);

%===============================================================================
% Compute the elements of the matrix using the exact formula from the book 
% of Zienkiewicz (8.38): 2*(a!b!c!/(a+b+c+2)!)
%===============================================================================

for aa = 1:N_degree
    for bb = aa:N_degree
        for cc = bb:N_degree
            LIntMatrix(aa,bb,cc) = 2*factorial(aa - 1)*factorial(bb - 1)/prod([ ( (cc  - 1) + 1 ) : ((aa  - 1) + (bb  - 1) + (cc  - 1) + 2)]);
            LIntMatrix(aa,cc,bb) = 2*factorial(aa - 1)*factorial(bb - 1)/prod([ ( (cc  - 1) + 1 ) : ((aa  - 1) + (bb  - 1) + (cc  - 1) + 2)]);
            LIntMatrix(bb,aa,cc) = 2*factorial(aa - 1)*factorial(bb - 1)/prod([ ( (cc  - 1) + 1 ) : ((aa  - 1) + (bb  - 1) + (cc  - 1) + 2)]);
            LIntMatrix(bb,cc,aa) = 2*factorial(aa - 1)*factorial(bb - 1)/prod([ ( (cc  - 1) + 1 ) : ((aa  - 1) + (bb  - 1) + (cc  - 1) + 2)]);
            LIntMatrix(cc,aa,bb) = 2*factorial(aa - 1)*factorial(bb - 1)/prod([ ( (cc  - 1) + 1 ) : ((aa  - 1) + (bb  - 1) + (cc  - 1) + 2)]);
            LIntMatrix(cc,bb,aa) = 2*factorial(aa - 1)*factorial(bb - 1)/prod([ ( (cc  - 1) + 1 ) : ((aa  - 1) + (bb  - 1) + (cc  - 1) + 2)]);
        end
    end
end

end