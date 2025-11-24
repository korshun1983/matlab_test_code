function [Results] = St4_ComputeSolution_sp_SAFE(CompStruct,BasicMatrices,FEMatrices,FullMatrices)

%   Part of the toolbox for solving problems of wave propagation
% in arbitrary anisotropic inhomogeneous waveguides.
% For details see User manual 
% and comments to the main script gen_aniso.m
%
% St4_ComputeSolution_sp_SM M-file      
%      St4_ComputeSolution_sp_SM, by itself, 
% organizes and carries out computation of the problem solution.
% (the waveguide spectrum using spectral method).
% It assumes that all basic blocks are available.
%
% NB!!! This implementation is specific for spectrum computation 
% with the spectral method.
%
%   [T.Zharnikov, D.Syresin, SMR v0.12_12.2012]
%
% function [Results] = St4_ComputeSolution_sp_SM(CompStruct,BasicMatrices)
%
%  Inputs - 
%
%       CompStuct - structure containing the parameters of the model 
%                   in the form, which is convenient and ready for
%                   computations
%
%       BasicMatrices - structure containing the basic matrices (differentiation, 
%                   matrix representation blocks, etc.);
%
%  Outputs -
%
%       Results - structure containing the results of the computation. 
%                   It means eigenvalues and eigenvectors. The matrix 
%                   representations of the operators are included as well.
%
%  M-files required-
%
% Last Modified by Timur Zharnikov SMR v0.1_01.2012

%###############################################################################
%
%   Code for St4_ComputeSolution_sp_SM
%
%###############################################################################
%===============================================================================
% Initialization
%===============================================================================

%Results = struct([]);

% % Symmetrize the matices, if necessary - see Magliula, e.g.
% [M_full,M2_full,K1_full,K2_full,K3_full,M_Source] = symmetrization_new( M_full,M2_full,K1_full,K2_full,K3_full,M_Source,nodes);
% % Forse the symmetry of the matrices, again, if necessary - Denis Syresin
% % introduced it. Clear explanation why it is necessary and should be done
% % is missing though.
% K1_full=(triu(K1_full)+(triu(K1_full))'-diag(diag(K1_full)));

%===============================================================================
% Cycle for computing points on dispersion curve according to k or omega grids
%===============================================================================

    
    ii_f=CompStruct.if_grid;
    %==============================================================================
    % Computing and assembling mass and stiffness matrices for FE
    % computations and inserting boundary and interface conditions 
    %===============================================================================
    
    switch CompStruct.Config.OuterBC
        case 'fixed'
            
            % Assembling spectral method matrices

            % Inserting boundary and interface conditions into the spectral method matrices

        case 'adjust'
            
            % Adjust radial extent of the outer layer
            switch CompStruct.Model.LayerType{CompStruct.Data.N_layers}
                case 'HTTI'
                    CompStruct.Asymp = CompStruct.Methods.ComputeAsymptotes(CompStruct);
                    % Adjust the range of phase speeds to study, if necessary
                    if isfield(CompStruct.Asymp,'V_SH')
                        CompStruct.Advanced.V_ref = CompStruct.Asymp.V_SH;
                    else
                        CompStruct.Advanced.V_ref = ( CompStruct.Advanced.V_min + CompStruct.Advanced.V_max )/2;                            
                    end;
                    CompStruct.Advanced.lam_ref = CompStruct.Advanced.V_ref/CompStruct.Advanced.f_ref;
                    k_ref = 2*pi*CompStruct.Advanced.f_ref/CompStruct.Advanced.V_ref;
                otherwise
            end;
            
            switch CompStruct.Config.EigenVar
                case 'k'
                    %Compute spectral method operator for particular k
                    f_val = BasicMatrices.f_grid(ii_f);
                    CompStruct.Model.LayerR(CompStruct.Data.N_layers + 1) = ...
                            CompStruct.Model.LayerR(CompStruct.Data.N_layers) + ...
                            CompStruct.Advanced.RadialExtensionFactor*CompStruct.Advanced.lam_ref*CompStruct.Advanced.f_ref/f_val;
                case 'omega'
                otherwise
            end;
            
            % Adjust blocks of the spectral matrices for the outer layer
            ii_l = CompStruct.Data.N_layers;
            [BasicMatrices.DDR1{ii_l},BasicMatrices.DDR2{ii_l},...
                BasicMatrices.ZM{ii_l},BasicMatrices.IdM{ii_l},...
                BasicMatrices.DiagR1{ii_l},BasicMatrices.DiagR2{ii_l},...
                BasicMatrices.RMatrix{ii_l}] = CompStruct.Methods.AssembleBasicMatrices(CompStruct,ii_l);
            
            [BasicMatrices.PhysProp{ii_l}] = CompStruct.Methods.PreparePhysProp{ii_l}(CompStruct,ii_l);
            
            [BasicMatrices.L_parts{ii_l},...
                BasicMatrices.S_parts{ii_l},...
                BasicMatrices.D_parts{ii_l},...
                BasicMatrices.M_parts{ii_l}] = CompStruct.Methods.AssembleMatricesParts_sp_SM{ii_l}(BasicMatrices,CompStruct,ii_l);
            

            % Assembling spectral method matrices

            % Inserting boundary and interface conditions into the spectral method matrices

        case 'PML'
            
            % Adjust radial extent of the outer layer
            switch CompStruct.Model.LayerType{CompStruct.Data.N_layers}
                case 'HTTI'
                    CompStruct.Asymp = CompStruct.Methods.ComputeAsymptotes(CompStruct);
                    % Adjust the range of phase speeds to study, if necessary
                    if isfield(CompStruct.Asymp,'V_SH')
                        CompStruct.Advanced.V_ref = CompStruct.Asymp.V_SH;
                    else
                        CompStruct.Advanced.V_ref = ( CompStruct.Advanced.V_min + CompStruct.Advanced.V_max )/2;                            
                    end;
                    CompStruct.Advanced.lam_ref = CompStruct.Advanced.V_ref/CompStruct.Advanced.f_ref;
                    k_ref = 2*pi*CompStruct.Advanced.f_ref/CompStruct.Advanced.V_ref;
                otherwise
            end;
            
            switch CompStruct.Config.EigenVar
                case 'k'
                    %Compute spectral method operator for particular k
                    f_val = BasicMatrices.f_grid(ii_f);
                    CompStruct.Model.LayerR(CompStruct.Data.N_layers) = ...
                            CompStruct.Model.LayerR(CompStruct.Data.N_layers - 1) + 3*CompStruct.Advanced.lam_ref*CompStruct.Advanced.f_ref/f_val;
                    CompStruct.Model.LayerR(CompStruct.Data.N_layers + 1) = ...
                            CompStruct.Model.LayerR(CompStruct.Data.N_layers) + 3*CompStruct.Advanced.lam_ref*CompStruct.Advanced.f_ref/f_val;
                case 'omega'
                otherwise
            end;
            
            % Adjust blocks of the spectral matrices for the preouter layer
            omega_val = 2*pi*BasicMatrices.f_grid(ii_f); %????????????????????????

            ii_l = CompStruct.Data.N_layers - 1;
            [BasicMatrices.DDR1{ii_l},BasicMatrices.DDR2{ii_l},...
                BasicMatrices.ZM{ii_l},BasicMatrices.IdM{ii_l},...
                BasicMatrices.DiagR1{ii_l},BasicMatrices.DiagR2{ii_l},...
                BasicMatrices.RMatrix{ii_l}] = CompStruct.Methods.AssembleBasicMatrices(CompStruct,ii_l);
            
            [BasicMatrices.PhysProp{ii_l}] = CompStruct.Methods.PreparePhysProp{ii_l}(CompStruct,ii_l);
            
            [BasicMatrices.L_parts{ii_l},...
                BasicMatrices.S_parts{ii_l},...
                BasicMatrices.D_parts{ii_l},...
                BasicMatrices.M_parts{ii_l}] = CompStruct.Methods.AssembleMatricesParts_sp_SM{ii_l}(BasicMatrices,CompStruct,ii_l);


            % Adjust blocks of the spectral matrices for the outer layer
            CompStruct.Advanced.PML_factor = 3/2*CompStruct.Advanced.V_ref*4.0...
                /( CompStruct.Model.LayerR(CompStruct.Data.N_layers + 1) - CompStruct.Model.LayerR(CompStruct.Data.N_layers) );
            omega_val = 2*pi*BasicMatrices.f_grid(ii_f); %????????????????????????

            ii_l = CompStruct.Data.N_layers;
            [BasicMatrices.DDR1{ii_l},BasicMatrices.DDR2{ii_l},...
                BasicMatrices.ZM{ii_l},BasicMatrices.IdM{ii_l},...
                BasicMatrices.DiagR1{ii_l},BasicMatrices.DiagR2{ii_l},...
                BasicMatrices.RMatrix{ii_l}] = CompStruct.Methods.AssembleBasicMatrices(CompStruct,ii_l);
            
            [BasicMatrices.PhysProp{ii_l}] = CompStruct.Methods.PreparePhysProp{ii_l}(CompStruct,ii_l);
            
            [BasicMatrices.L_parts{ii_l},...
                BasicMatrices.S_parts{ii_l},...
                BasicMatrices.D_parts{ii_l},...
                BasicMatrices.M_parts{ii_l}] = CompStruct.Methods.AssembleMatricesPartsPML_sp_SM{ii_l}(omega_val,BasicMatrices,CompStruct,ii_l);
            

            % Assembling spectral method matrices

            % Inserting boundary and interface conditions into the spectral method matrices

        case 'ABC'
            
            % Adjust radial extent of the outer layer
            switch CompStruct.Model.LayerType{CompStruct.Data.N_layers}
                case 'HTTI'
                    CompStruct.Asymp = CompStruct.Methods.ComputeAsymptotes(CompStruct);
                    % Adjust the range of phase speeds to study, if necessary
                    if isfield(CompStruct.Asymp,'V_SH')
                        CompStruct.Advanced.V_ref = CompStruct.Asymp.V_SH;
                    else
                        CompStruct.Advanced.V_ref = ( CompStruct.Advanced.V_min + CompStruct.Advanced.V_max )/2;                            
                    end;
                    CompStruct.Advanced.lam_ref = CompStruct.Advanced.V_ref/CompStruct.Advanced.f_ref;
                    k_ref = 2*pi*CompStruct.Advanced.f_ref/CompStruct.Advanced.V_ref;
                otherwise
            end;
            
            switch CompStruct.Config.EigenVar
                case 'k'
                    %Compute spectral method operator for particular k
                    f_val = BasicMatrices.f_grid(ii_f);
                    CompStruct.Model.LayerR(CompStruct.Data.N_layers) = ...
                            CompStruct.Model.LayerR(CompStruct.Data.N_layers - 1) + 5*CompStruct.Advanced.lam_ref*CompStruct.Advanced.f_ref/f_val;
                    CompStruct.Model.LayerR(CompStruct.Data.N_layers + 1) = ...
                            CompStruct.Model.LayerR(CompStruct.Data.N_layers) + 5*CompStruct.Advanced.lam_ref*CompStruct.Advanced.f_ref/f_val;
                case 'omega'
                otherwise
            end;
            
            % Adjust blocks of the spectral matrices for the preouter layer
            omega_val = 2*pi*BasicMatrices.f_grid(ii_f); %????????????????????????

            ii_l = CompStruct.Data.N_layers - 1;
            [BasicMatrices.DDR1{ii_l},BasicMatrices.DDR2{ii_l},...
                BasicMatrices.ZM{ii_l},BasicMatrices.IdM{ii_l},...
                BasicMatrices.DiagR1{ii_l},BasicMatrices.DiagR2{ii_l},...
                BasicMatrices.RMatrix{ii_l}] = CompStruct.Methods.AssembleBasicMatrices(CompStruct,ii_l);
            
            [BasicMatrices.PhysProp{ii_l}] = CompStruct.Methods.PreparePhysProp{ii_l}(CompStruct,ii_l);
            
            [BasicMatrices.L_parts{ii_l},...
                BasicMatrices.S_parts{ii_l},...
                BasicMatrices.D_parts{ii_l},...
                BasicMatrices.M_parts{ii_l}] = CompStruct.Methods.AssembleMatricesParts_sp_SM{ii_l}(BasicMatrices,CompStruct,ii_l);

            % Adjust blocks of the spectral matrices for the outer layer
            omega_val = 2*pi*BasicMatrices.f_grid(ii_f); %????????????????????????

            ii_l = CompStruct.Data.N_layers;
            [BasicMatrices.DDR1{ii_l},BasicMatrices.DDR2{ii_l},...
                BasicMatrices.ZM{ii_l},BasicMatrices.IdM{ii_l},...
                BasicMatrices.DiagR1{ii_l},BasicMatrices.DiagR2{ii_l},...
                BasicMatrices.RMatrix{ii_l}] = CompStruct.Methods.AssembleBasicMatrices(CompStruct,ii_l);
            
            [BasicMatrices.PhysProp{ii_l}] = CompStruct.Methods.PreparePhysProp{ii_l}(CompStruct,ii_l);
            
            [BasicMatrices.L_parts{ii_l},...
                BasicMatrices.S_parts{ii_l},...
                BasicMatrices.D_parts{ii_l},...
                BasicMatrices.M_parts{ii_l}] = CompStruct.Methods.AssembleMatricesPartsABC_sp_SM{ii_l}(BasicMatrices,CompStruct,ii_l);
            

            % Assembling spectral method matrices

            % Inserting boundary and interface conditions into the spectral method matrices

        case 'TTBC'
            % to be added, if possible to implement
        otherwise
    end;

    %==============================================================================
    % Formulating generalized eigenvalue problem
    %===============================================================================
    switch CompStruct.Config.EigenVar 
        case 'k'
            
            % The problem statement for the SAFE approach is:
            % K1 + 1i*k*K2 + k^2*K3 - omega^2*M - 1i*omega*P
            % It is equivalent to the polynomial evp
            % M*lambda^2 + C*lambda + K
            %
            % To linearize this problem as:
            % A - lambda*B = 0
            % NB!!! To use eigs, B should positive definite or Hermitian
            % Two possible linearizations of this evp are as follows:
            % (1) [ 0 N ; -K -C ] - lambda*[ N 0 ; 0 M ]
            % (2) [ -K 0 ; 0 N ] - lambda*[ C M ; N 0 ]
            % where N can be any nonsingular n × n matrix
            %
            % Taking k as the eigenvalue and omega as the parameter, one gets:
            % M = K3 ; C = 1i*K2 ; K = K1 - M*omega^2 - 1i*P*omega
            % 
            % Given the problem formulation, the linearization (1) is more
            % convenient. 
            % The simplest choice of N is either the identity matrix,
            % or (-K)
            
            Msize = size(FullMatrices.K1Matrix);
            ZMFull = sparse(Msize(1),Msize(2));
            IdMFull = speye(Msize(1),Msize(2));%normest(FullMatrices.K3Matrix)*speye(Msize(1),Msize(2));%

            % Compute spectral method operator for particular omega
            omega_val = 2*pi*BasicMatrices.f_grid(ii_f);
            omega_var_1000=omega_val*1.e3;
            switch CompStruct.Config.SpeedUp
                case 'yes'
                case 'no'
%                     Results.SMOpFull = sparse(Results.SMOp{1}(:,:) + k_val*Results.SMOp{2}(:,:) + k_val^2*Results.SMOp{3}(:,:));
%                     Results.MMatrixFull = sparse(Results.MMatrix);

                    % N = IdMFull
                    % alternatively 
                    % N = -K'

                    FullMatrices.KMatrix = FullMatrices.K1Matrix ...
                        - omega_var_1000^2*FullMatrices.MMatrix + omega_var_1000*1i*FullMatrices.PMatrix;
                    FullMatrices.AMatrix  = [ ZMFull ,                       IdMFull ; ... 
                                              ( - FullMatrices.KMatrix ),    ( - 1i*FullMatrices.K2Matrix ) ];
                    FullMatrices.BMatrix = [ IdMFull,    ZMFull ; ...
                                             ZMFull,     FullMatrices.K3Matrix ];
%  A=[Z K3norm; 
%      -L -K2_full];  
%  B=[K3norm Z;
%      Z K3_full];   
% L=(K1_full-(omega_ar(i))^2*M_full+omega_ar(i)*M2_full);
                otherwise
            end;
        case 'omega'
            % For SAFE the computations are done currently with 'k' as the
            % eigenvalue, because it gives more flexibility - possibility
            % to account for the leaky modes, frequency dependent
            % attenuation, etc.
            % So, 'omega' option is neither considered nor implemented.
            % 
            %Compute spectral method operator for particular k
%             k_val = BasicMatrices.k_grid(ii_f);
%             switch CompStruct.Config.SpeedUp
%                 case 'yes'
%                     Results.SMOpFull = sparse(Results.SMOp{1}(:,:) + k_val*Results.SMOp{2}(:,:) + k_val^2*Results.SMOp{3}(:,:));
%                     Results.MMatrixFull = sparse(Results.MMatrix);
%                 case 'no'
%                     Results.SMOpFull = Results.SMOp{1}(:,:) + k_val*Results.SMOp{2}(:,:) + k_val^2*Results.SMOp{3}(:,:);
%                     Results.MMatrixFull = Results.MMatrix;
%                 otherwise
%             end;
        otherwise
    end;
    
    %===============================================================================
    % Finding the spectrum of generalized eigenvalue problem
    % NB! Since SAFE works with the sparse matrices (otherwise the
    % computation time is prohibitively expensive),
    % the procedure always calls 'eigs' and there is no need for speed up
    %===============================================================================
    % Compute the wave vector value, from which to start searching
    % for the eigenvalues
    k_start_val = omega_val/CompStruct.Advanced.EigSearchStart;

%    [Results.REig_vecs,Results.Eig_vals,Results.LEig_vecs] = eigs(FullMatrices.AMatrix,FullMatrices.BMatrix,...
    [Results.REig_vecs,Results.REig_vals] = eigs(FullMatrices.AMatrix,FullMatrices.BMatrix,...
        CompStruct.Advanced.num_eig_max,k_start_val,CompStruct.Advanced.EigsOptions);
    
    % we make additionalsorting of eigs program result due to possible mistake of program in sorting
    [~,ind_sort]=sort(abs(diag(Results.REig_vals)),'descend');
    Results.REig_vals=Results.REig_vals(ind_sort,ind_sort);
    Results.REig_vecs=Results.REig_vecs(:,ind_sort);

    %% Find left eigvectors
    %[Results.REig_Lvecs,Results.REig_Lvals] = eigs(FullMatrices.AMatrix',FullMatrices.BMatrix',...
    %    CompStruct.Advanced.num_eig_max,k_start_val,CompStruct.Advanced.EigsOptions);
    %% we make additionalsorting of eigs program result due to possible mistake of program in sorting
    %[~,ind_sort]=sort(abs(diag(Results.REig_Lvals)),'descend');
    %Results.REig_Lvals=Results.REig_Lvals(ind_sort,ind_sort);
    %Results.REig_Lvecs=Results.REig_Lvecs(:,ind_sort);

    
%     FullMatrices.AMatrix_H  = [ ZMFull ,                       IdMFull ; ...
%         ( - FullMatrices.KMatrix )',    ( - 1i*FullMatrices.K2Matrix )' ];
%     FullMatrices.BMatrix_H = [ IdMFull,    ZMFull ; ...
%         ZMFull,     (FullMatrices.K3Matrix)' ];
%     
%     [Results.LEig_vecs,Results.LEig_vals] = eigs(FullMatrices.AMatrix_H,FullMatrices.BMatrix_H,...
%         CompStruct.Advanced.num_eig_max,k_start_val,CompStruct.Advanced.EigsOptions);

%     toc
    
    %===============================================================================
    % Assembling results
    %===============================================================================
%     Results.Rho_vec = [];
%     for ii_l = 1:CompStruct.Data.N_layers
%         Results.Rho_vec = vertcat(Results.Rho_vec,BasicMatrices.PhysProp{ii_l}.rho_dist);
%     end;
%     Results.R_vec = vertcat(BasicMatrices.RMatrix{:});
%     Results.DiffR = [0; diff(Results.R_vec)];
    
    switch CompStruct.Config.SaveData % safe is perpormed in main program
        case 'yes'
            Results.omega_val = omega_val;
        case 'no'
            Results.Momega_val{ii_f} = omega_val;
            Results.MREig_vecs{ii_f} = Results.REig_vecs;
            Results.MREig_vals{ii_f} = Results.REig_vals;
%             Results.MLEig_vecs{ii_f} = Results.LEig_vecs;
%             Results.MLEig_vals{ii_f} = Results.LEig_vals;

%            Results.FullMatrices{ii_f} = FullMatrices;

%             Results.MSOp{ii_f} = Results.SOp;
%             Results.MDOp{ii_f} = Results.DOp;
%             Results.MR_vec{ii_f} = Results.R_vec;
%             Results.MRho_vec{ii_f} = Results.Rho_vec;
%             Results.MDiffR{ii_f} = Results.DiffR;
%            Results.k_grid{ii_f} = BasicMatrices.k_grid(ii_f);
    end
    
%---------------------------------------------------------------------------------

end
