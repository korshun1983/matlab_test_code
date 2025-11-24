function [BasicMatrices] = ConvolveEdgeMatrices(CompStruct,BasicMatrices)
     
%   Part of the toolbox for solving problems of wave propagation
% in arbitrary anisotropic inhomogeneous waveguides.
% For details see User manual 
% and comments to the main script gen_aniso.m
%
%   ConvolveMatrices function computes the matrices, which indicate the expansion
%  and the coefficients of the expansion of the double convolutions of the
%  types conv(N,conv(N,N)), conv(dN,conv(N,N)), conv(N,conv(N,dN)), conv(dN,conv(N,dN)),
%  into the intepolating functions L_j (L1, L2, L3).
%  These convolutions enter into the expressions for the elements 
%  of the integrands for K_i and M matrices' elements.
%  For HTTI solid.
%
%  Current implementation is for cubic elements.
% NB!!! This implementation is specific for FE computations
%
%   [T.Zharnikov, SMR v0.3_08.2014]
%
% function [NNNConvMatrix,dNNNConvMatrix,NNdNConvMatrix,dNNdNConvMatrix]...
%             = ConvolveMatrices(BasicMatrices)
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
%   Code for ConvolveMatrices
%
%###############################################################################
%===============================================================================
% Initialization
%===============================================================================
fname = strcat(CompStruct.Config.root_path,'output\ConvEdgeMatrices.mat');
%fname='';

if (exist(fname,'file') ~= 0)
    load(fname);
    
    BasicMatrices.NENENEConvMatrixInt = NENENEConvMatrixInt;
%     NNNConvMatrix_size = [ N_nodes,N_nodes,N_nodes,10,10,10 ];
%     dNNNConvMatrix_size = [ 3,N_nodes,N_nodes,N_nodes,10,10,10 ];
%     NNdNConvMatrix_size = [ 3,N_nodes,N_nodes,N_nodes,10,10,10 ];
%     dNNdNConvMatrix_size = [ 3,3,N_nodes,N_nodes,N_nodes,10,10,10 ];
    
%     NNNConvMatrixInt = reshape(full(NNNConvMatrixInt_sparse),NNNConvMatrix_size);
%     dNNNConvMatrixInt = reshape(full(dNNNConvMatrixInt_sparse),dNNNConvMatrix_size);
%     NNdNConvMatrixInt = reshape(full(NNdNConvMatrixInt_sparse),NNdNConvMatrix_size);
%     dNNdNConvMatrixInt = reshape(full(dNNdNConvMatrixInt_sparse),dNNdNConvMatrix_size);
    
else

    % Identify number of nodes
    NEdge_nodes = CompStruct.Advanced.NEdge_nodes;
    
    % Allocate memory for the matrices.
    % The implementation is for cubic elements, hence expansion is up to the
    % 3rd power in L_i
    
    NENENEConvMatrix = zeros(NEdge_nodes,NEdge_nodes,NEdge_nodes,10,10);
%     dNNNConvMatrix = zeros(3,N_nodes,N_nodes,N_nodes,10,10,10);
%     NNdNConvMatrix = zeros(3,N_nodes,N_nodes,N_nodes,10,10,10);
%     dNNdNConvMatrix = zeros(3,3,N_nodes,N_nodes,N_nodes,10,10,10);
    
    LIntMatrix = BasicMatrices.LEdgeIntMatrix9(1:10,1:10);
    LIntMatrixNENENE = zeros(size(NENENEConvMatrix));
%     LIntMatrixNNdN = zeros(size(NNdNConvMatrix));
%     LIntMatrixdNNN = zeros(size(NNdNConvMatrix));
%     LIntMatrixdNNdN = zeros(size(dNNdNConvMatrix));
    
    %===============================================================================
    % Compute the elements of all of the convolution matrices
    % and prepare large LIntMatrices for the integration
    %===============================================================================
    
    % process all of the elements of the
    % NNNConvMatrix, dNNNConvMatrix, NNdNConvMatrix, dNNdNConvMatrix matrices
    'conv';
    %tic
    for bb = 1:NEdge_nodes
        for cc = 1:NEdge_nodes
            NEb = squeeze(BasicMatrices.NLEdgeMatrix(bb,:,:));
            NEc = squeeze(BasicMatrices.NLEdgeMatrix(cc,:,:));
            convn_NEb_NEc = convn(NEb,NEc);
            
            for aa = 1:NEdge_nodes
                NEa = squeeze(BasicMatrices.NLEdgeMatrix(aa,:,:));
                NENENEConvMatrix(aa,bb,cc,:,:) = ...
                    convn(NEa,convn_NEb_NEc);
                
                LIntMatrixNENENE(aa,bb,cc,:,:) = LIntMatrix;
                
%                 for ii = 1:3
%                     dNa_i = squeeze(BasicMatrices.dNLMatrix(ii,aa,:,:,:));
%                     dNNNConvMatrix(ii,aa,bb,cc,:,:,:) = ...
%                         convn(dNa_i,convn_Nb_Nc);
%                     
%                     LIntMatrixdNNN(ii,aa,bb,cc,:,:,:) = LIntMatrix;
%                 end;
            end
            
%             for jj = 1:3
%                 dNc_j = squeeze(BasicMatrices.dNLMatrix(jj,cc,:,:,:));
%                 convn_Nb_dNc_j = convn(Nb,dNc_j);
%                 
%                 for aa = 1:N_nodes
%                     Na = squeeze(BasicMatrices.NLMatrix(aa,:,:,:));
%                     NNdNConvMatrix(jj,aa,bb,cc,:,:,:) = ...
%                         convn(Na,convn_Nb_dNc_j);
%                     
%                     LIntMatrixNNdN(jj,aa,bb,cc,:,:,:) = LIntMatrix;
%                     
%                     for ii = 1:3
%                         dNa_i = squeeze(BasicMatrices.dNLMatrix(ii,aa,:,:,:));
%                         dNNdNConvMatrix(ii,jj,aa,bb,cc,:,:,:) = ...
%                             convn(dNa_i,convn_Nb_dNc_j);
%                         
%                         LIntMatrixdNNdN(ii,jj,aa,bb,cc,:,:,:) = LIntMatrix;
%                     end;
%                 end;
%             end;
            
        end
    end
    %toc;
    %===============================================================================
    % Compute the integrals over the standard triangle of the monoms of the form
    % N_a*N_b*N_c, dN_a*N_b*N_c, N_a*N_b*dN_c, dN_a*N_b*dN_c
    % by multiplying elementwise convolution matrices NNN, etc. on the large
    % LIntMatrices and summing over the L1^i*L2^j*L3^k contributions
    %===============================================================================
    'conv_int';
    %tic;
    NENENEConvMatrix_unsummed = NENENEConvMatrix.*LIntMatrixNENENE;
    BasicMatrices.NENENEConvMatrixInt = squeeze(sum(sum(NENENEConvMatrix_unsummed,5),4));
%     dNNNConvMatrix_unsummed = dNNNConvMatrix.*LIntMatrixdNNN;
%     dNNNConvMatrixInt = squeeze(sum(sum(sum(dNNNConvMatrix_unsummed,7),6),5));
%     NNdNConvMatrix_unsummed = NNdNConvMatrix.*LIntMatrixNNdN;
%     NNdNConvMatrixInt = squeeze(sum(sum(sum(NNdNConvMatrix_unsummed,7),6),5));
%     dNNdNConvMatrix_unsummed = dNNdNConvMatrix.*LIntMatrixdNNdN;
%     dNNdNConvMatrixInt = squeeze(sum(sum(sum(dNNdNConvMatrix_unsummed,8),7),6));
    
    % Save precomputed NNN integral matrices to the external file
%     NNNConvMatrix_s1 = N_nodes*N_nodes*N_nodes;
%     dNNNConvMatrix_s1 = 3*N_nodes*N_nodes*N_nodes;
%     NNdNConvMatrix_s1 = 3*N_nodes*N_nodes*N_nodes;
%     dNNdNConvMatrix_s1 = 3*3*N_nodes*N_nodes*N_nodes;
%     NNNConvMatrix_s2 = 10*10*10;
%     dNNNConvMatrix_s2 = 10*10*10;
%     NNdNConvMatrix_s2 = 10*10*10;
%     dNNdNConvMatrix_s2 = 10*10*10;

%     NNNConvMatrixInt_sparse = sparse(reshape(NNNConvMatrixInt,NNNConvMatrix_s1,NNNConvMatrix_s2));
%     dNNNConvMatrixInt_sparse = sparse(reshape(dNNNConvMatrixInt,dNNNConvMatrix_s1,dNNNConvMatrix_s2));
%     NNdNConvMatrixInt_sparse = sparse(reshape(NNdNConvMatrixInt,NNdNConvMatrix_s1,NNdNConvMatrix_s2));
%     dNNdNConvMatrixInt_sparse = sparse(reshape(dNNdNConvMatrixInt,dNNdNConvMatrix_s1,dNNdNConvMatrix_s2));
    save(fname,'-struct','BasicMatrices','NENENEConvMatrixInt');%,'dNNNConvMatrixInt','NNdNConvMatrixInt','dNNdNConvMatrixInt');
    %toc;

end;

end
