function [ElMatrices] = KM_el_matrix_HTTI_PML(BasicMatrices, CompStruct, FEMatrices, ... 
                                             ii_d, ii_el, ElPhysProps, TriProps)

% NB!!!!!!! WE HAVE PML LAYER
     
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

% Allocate memory for the matrices.
% The implementation is for cubic elements, hence expansion is up to the 3rd power in L_i
B1Matrix = zeros(6,CompStruct.Data.DVarNum(ii_d)*N_nodes,4,4,4);
B2Matrix = zeros(6,CompStruct.Data.DVarNum(ii_d)*N_nodes,4,4,4);
NMatrix = zeros(3,CompStruct.Data.DVarNum(ii_d)*N_nodes,4,4,4);

ElMatrices.B1tCB1Matrix = zeros(CompStruct.Data.DVarNum(ii_d)*N_nodes,CompStruct.Data.DVarNum(ii_d)*N_nodes);
ElMatrices.B1tCB2Matrix = zeros(CompStruct.Data.DVarNum(ii_d)*N_nodes,CompStruct.Data.DVarNum(ii_d)*N_nodes);
ElMatrices.B2tCB1Matrix = zeros(CompStruct.Data.DVarNum(ii_d)*N_nodes,CompStruct.Data.DVarNum(ii_d)*N_nodes);
ElMatrices.B2tCB2Matrix = zeros(CompStruct.Data.DVarNum(ii_d)*N_nodes,CompStruct.Data.DVarNum(ii_d)*N_nodes);

NtRhoNMatrix = zeros(CompStruct.Data.DVarNum(ii_d)*N_nodes,CompStruct.Data.DVarNum(ii_d)*N_nodes);

%===============================================================================
% Compute the elements of K_i and M matrices according to the above presented formulas
%===============================================================================
% Prepare various properties
Rho = ElPhysProps.RhoVec;
CijMatrix = ElPhysProps.CijMatrix;
dxL = TriProps.dxL;
dyL = TriProps.dyL;
Lx = BasicMatrices.Lx;
Ly = BasicMatrices.Ly;
Lz = BasicMatrices.Lz;

OnesMatrix = ones(N_nodes);
OnesVarNum = ones(CompStruct.Data.DVarNum(ii_d));
% OnesVarMatrix = ones(CompStruct.Data.DVarNum(ii_d));
IdMatrix = eye(CompStruct.Data.DVarNum(ii_d));
Msize = CompStruct.Data.DVarNum(ii_d)*N_nodes;

% process all of the elements of the stiffness and mass matrices 
% expansion coefficients (N_nodes x N_nodes x ... x ... x ... )
Rho_large = zeros([Msize N_nodes Msize]);
% Rho_large = zeros(size(BasicMatrices.NNNConvMatrixInt));
LxTCijLx_dxL_dxL_large = zeros([3 3 Msize N_nodes Msize]);
% LxTCijLx_dxL_dxL_large = zeros(size(BasicMatrices.dNNdNConvMatrixInt));
LxTCijLz_dxL_large = zeros([3 Msize N_nodes Msize]);
% LxTCijLz_dxL_large = zeros(size(BasicMatrices.dNNNConvMatrixInt));
LzTCijLx_dxL_large = zeros([3 Msize N_nodes Msize]);
% LzTCijLx_dxL_large = zeros(size(BasicMatrices.NNdNConvMatrixInt));
LzTCijLz_large = zeros([Msize N_nodes Msize]);
% LzTCijLz_large = zeros(size(BasicMatrices.NNNConvMatrixInt));
% dNNdNConvMatrixIntLarge = zeros([3 3 Msize N_nodes Msize]);
% dNNNConvMatrixIntLarge = zeros([3 Msize N_nodes Msize]);
% NNdNConvMatrixIntLarge = zeros([3 Msize N_nodes Msize]);
% NNNConvMatrixIntLarge = zeros([Msize N_nodes Msize]);
% %NodeNumbers = linspace(1,N_nodes,N_nodes);
% for ii_k = 1:CompStruct.Advanced.N_nodes
%     NNNConvMatrixIntLarge(:,ii_k,:) = kron(squeeze(BasicMatrices.NNNConvMatrixInt(:,ii_k,:)),OnesVarNum);
%     for ii_c = 1:3
%         NNdNConvMatrixIntLarge(ii_c,:,ii_k,:) = ...
%             kron(squeeze(BasicMatrices.NNdNConvMatrixInt(ii_c,:,ii_k,:)),OnesVarNum);
%         dNNNConvMatrixIntLarge(ii_c,:,ii_k,:) = ...
%             kron(squeeze(BasicMatrices.dNNNConvMatrixInt(ii_c,:,ii_k,:)),OnesVarNum);
%         for ii_c2 = 1:3
%             dNNdNConvMatrixIntLarge(ii_c,ii_c2,:,ii_k,:) = ...
%                 kron(squeeze(BasicMatrices.dNNdNConvMatrixInt(ii_c,ii_c2,:,ii_k,:)),OnesVarNum);
%         end;
%     end;
% end;

%     NNNConvMatrixIntLarge( ( ii + CompStruct.Data.DVarNum(ii_d)*( NodeNumbers - 1 )),:,...
%         ( ii + CompStruct.Data.DVarNum(ii_d)*( NodeNumbers - 1 )) ) = ...
%         BasicMatrices.NNNConvMatrixInt;
%     NNdNConvMatrixIntLarge(:, ( ii + CompStruct.Data.DVarNum(ii_d)*( NodeNumbers - 1 )),:,...
%         ( ii + CompStruct.Data.DVarNum(ii_d)*( NodeNumbers - 1 )) ) = ...
%         BasicMatrices.NNdNConvMatrixInt;
%     dNNNConvMatrixIntLarge(:, ( ii + CompStruct.Data.DVarNum(ii_d)*( NodeNumbers - 1 )),:,...
%         ( ii + CompStruct.Data.DVarNum(ii_d)*( NodeNumbers - 1 )) ) = ...
%         BasicMatrices.dNNNConvMatrixInt;
%     dNNdNConvMatrixIntLarge(:,:, ( ii + CompStruct.Data.DVarNum(ii_d)*( NodeNumbers - 1 )),:,...
%         ( ii + CompStruct.Data.DVarNum(ii_d)*( NodeNumbers - 1 )) ) = ...
%         BasicMatrices.dNNdNConvMatrixInt;
% 

TriNodes = FEMatrices.DElements{ii_d}(1:10,ii_el);

% PML layer, its thickness is always in SH wavelength
% gamma(r)=1+sigma(r)/iw
% sigma(r)=H(r-r_zv)*alpha*(r-r_zv)^2
% H is Heaviside step function, r_zv is the inner radius of PML layer
% alpha is the constant: 1) =8*2pi*f_peak/L_pml^2, L_pml is the length of PML layer

f_peak=CompStruct.f_grid(CompStruct.if_grid).*CompStruct.Misc.F_conv; % Hz
w_peak=2.*pi*f_peak;
cone=complex(0.,1.);

r_zv=CompStruct.Model.DomainRx(end-1);
L_ext_layer=CompStruct.Model.DomainRx(end)-r_zv;

for kk = 1:N_nodes

    Lx = BasicMatrices.Lx;
    Ly = BasicMatrices.Ly;
    var_df_xy=1.;
    
    if ii_d==CompStruct.Data.N_domain   

       xy_node= FEMatrices.MeshNodes(1:2,TriNodes(kk));
      
       if CompStruct.Model.PML_method==1 % Circle PML

          r_node=sqrt(xy_node(1)^2+xy_node(2)^2);

          if strcmpi(CompStruct.Model.AddDomainLoc,'ext')
             r_zv=min(CompStruct.Model.DomainRx(end-1),CompStruct.Model.DomainRy(end-1));
             Lr_ext_layer=min(CompStruct.Model.DomainRx(end),CompStruct.Model.DomainRy(end))-r_zv;
          end

          if strcmpi(CompStruct.Model.AddDomainLoc,'int')
             Lr_ext_layer=CompStruct.Model.AddDomainL_m;
             r_zv=min(CompStruct.Model.DomainRx(end),CompStruct.Model.DomainRy(end))-Lr_ext_layer;
          end

          if r_node>r_zv
             ivar=CompStruct.Model.PML_degree;

             sigma_r=((r_node-r_zv)/Lr_ext_layer)^ivar;
             gamma_r=cone*CompStruct.Model.PML_factor*sigma_r;
             r_tilde=r_node-cone*(CompStruct.Model.PML_factor/(ivar+1))*((r_node-r_zv)/Lr_ext_layer)^(ivar+1);

             var_xy=xy_node(1)^2/(gamma_r*r_node^2)+xy_node(2)^2/(r_tilde*r_node);
             var_yx=xy_node(2)^2/(gamma_r*r_node^2)+xy_node(1)^2/(r_tilde*r_node);
             var_all=(1./(gamma_r*r_node^2)-1./(r_tilde*r_node))*xy_node(1)*xy_node(2);

             Lx=var_xy.*BasicMatrices.Lx+var_all.*BasicMatrices.Ly;
             Ly=var_all.*BasicMatrices.Lx+var_yx.*BasicMatrices.Ly;

             var_df_xy=gamma_r*r_tilde/r_node;
          end
       end

       if CompStruct.Model.PML_method==2  % Rectange PML
          ivar=CompStruct.Model.PML_degree;
           
          if strcmpi(CompStruct.Model.AddDomainLoc,'ext')
             x_zv=CompStruct.Model.DomainRx(end-1);
             Lx_ext_layer=CompStruct.Model.DomainRx(end)-x_zv;
             y_zv=CompStruct.Model.DomainRy(end-1);
             Ly_ext_layer=CompStruct.Model.DomainRy(end)-y_zv;
          end

          if strcmpi(CompStruct.Model.AddDomainLoc,'int')
             Lx_ext_layer=CompStruct.Model.AddDomainL_m;
             x_zv=CompStruct.Model.DomainRx(end)-Lx_ext_layer;
             Ly_ext_layer=CompStruct.Model.AddDomainL_m;
             y_zv=CompStruct.Model.DomainRx(end)-Ly_ext_layer;
          end

          gamma_x=1.;
          if abs(xy_node(1))>x_zv
             sigma_x=( ( abs(xy_node(1))-x_zv)/Lx_ext_layer )^ivar;
             gamma_x=1.-cone*CompStruct.Model.PML_factor*sigma_x;

             LL_xx=1./gamma_x;
             Lx=LL_xx.*BasicMatrices.Lx;
          end
           
          gamma_y=1.;
          if abs(xy_node(2))>y_zv
             sigma_y=( ( abs(xy_node(2))-y_zv)/Ly_ext_layer )^ivar;
             gamma_y=1.-cone*CompStruct.Model.PML_factor*sigma_y;

             LL_yy=1./gamma_y;
             Ly=LL_yy.*BasicMatrices.Ly;
          end
          
          var_df_xy=gamma_x*gamma_y;
       end

    end
    
%     !!!!!!!!!!!!!!!!!!!
    LxTCijLx_kk = conj(Lx')*squeeze(CijMatrix(kk,:,:))*Lx;
    LxTCijLy_kk = conj(Lx')*squeeze(CijMatrix(kk,:,:))*Ly;
    LyTCijLx_kk = conj(Ly')*squeeze(CijMatrix(kk,:,:))*Lx;
    LyTCijLy_kk = conj(Ly')*squeeze(CijMatrix(kk,:,:))*Ly;
    LxTCijLz_kk = conj(Lx')*squeeze(CijMatrix(kk,:,:))*Lz;
    LzTCijLx_kk = conj(Lz')*squeeze(CijMatrix(kk,:,:))*Lx;
    LyTCijLz_kk = conj(Ly')*squeeze(CijMatrix(kk,:,:))*Lz;
    LzTCijLy_kk = conj(Lz')*squeeze(CijMatrix(kk,:,:))*Ly;
    LzTCijLz_kk = conj(Lz')*squeeze(CijMatrix(kk,:,:))*Lz;
 
%     LxTCijLx_kk = Lx'*squeeze(CijMatrix(kk,:,:))*Lx;
%     LxTCijLy_kk = Lx'*squeeze(CijMatrix(kk,:,:))*Ly;
%     LyTCijLx_kk = Ly'*squeeze(CijMatrix(kk,:,:))*Lx;
%     LyTCijLy_kk = Ly'*squeeze(CijMatrix(kk,:,:))*Ly;
%     LxTCijLz_kk = Lx'*squeeze(CijMatrix(kk,:,:))*Lz;
%     LzTCijLx_kk = Lz'*squeeze(CijMatrix(kk,:,:))*Lx;
%     LyTCijLz_kk = Ly'*squeeze(CijMatrix(kk,:,:))*Lz;
%     LzTCijLy_kk = Lz'*squeeze(CijMatrix(kk,:,:))*Ly;
%     LzTCijLz_kk = Lz'*squeeze(CijMatrix(kk,:,:))*Lz;
    
    LxTCijLx_kk=LxTCijLx_kk*var_df_xy;
    LxTCijLy_kk=LxTCijLy_kk*var_df_xy;
    LyTCijLx_kk=LyTCijLx_kk*var_df_xy;
    LyTCijLy_kk=LyTCijLy_kk*var_df_xy;
    LxTCijLz_kk=LxTCijLz_kk*var_df_xy;
    LzTCijLx_kk=LzTCijLx_kk*var_df_xy;
    LyTCijLz_kk=LyTCijLz_kk*var_df_xy;
    LzTCijLy_kk=LzTCijLy_kk*var_df_xy;
    LzTCijLz_kk=LzTCijLz_kk*var_df_xy;
     
    Rho_large(:,kk,:) = reshape(kron(OnesMatrix,Rho(kk)*IdMatrix),Msize,1,Msize);
    Rho_large(:,kk,:) = Rho_large(:,kk,:)*var_df_xy;

    LzTCijLz_large(:,kk,:) = kron(OnesMatrix,LzTCijLz_kk);
    for ii = 1:3
        LxTCijLz_dxL_large(ii,:,kk,:) = kron(OnesMatrix,...
            ( LxTCijLz_kk*dxL(ii) + ...
              LyTCijLz_kk*dyL(ii) ));
        LzTCijLx_dxL_large(ii,:,kk,:) = kron(OnesMatrix,...
            ( LzTCijLx_kk*dxL(ii) + ...
              LzTCijLy_kk*dyL(ii) ));
        for jj = 1:3
            LxTCijLx_dxL_dxL_large(ii,jj,:,kk,:) = kron(OnesMatrix,...
                ( LxTCijLx_kk*dxL(ii)*dxL(jj) + ... 
                  LxTCijLy_kk*dxL(ii)*dyL(jj) + ...
                  LyTCijLx_kk*dyL(ii)*dxL(jj) + ...
                  LyTCijLy_kk*dyL(ii)*dyL(jj) ));
        end
    end
end

NtRhoNMatrix_unsummed = Rho_large.*BasicMatrices.NNNConvMatrixIntLarge{ii_d};
ElMatrices.NtRhoNMatrix = squeeze(sum(NtRhoNMatrix_unsummed,2));

B1tCB1Matrix_unsummed = LxTCijLx_dxL_dxL_large.*BasicMatrices.dNNdNConvMatrixIntLarge{ii_d};
ElMatrices.B1tCB1Matrix = squeeze(sum(sum(sum(B1tCB1Matrix_unsummed,1),2),4));
B1tCB2Matrix_unsummed = LxTCijLz_dxL_large.*BasicMatrices.dNNNConvMatrixIntLarge{ii_d};
ElMatrices.B1tCB2Matrix = squeeze(sum(sum(B1tCB2Matrix_unsummed,1),3));
B2tCB1Matrix_unsummed = LzTCijLx_dxL_large.*BasicMatrices.NNdNConvMatrixIntLarge{ii_d};
ElMatrices.B2tCB1Matrix = squeeze(sum(sum(B2tCB1Matrix_unsummed,1),3));
B2tCB2Matrix_unsummed = LzTCijLz_large.*BasicMatrices.NNNConvMatrixIntLarge{ii_d};
ElMatrices.B2tCB2Matrix = squeeze(sum(B2tCB2Matrix_unsummed,2));

% % process all of the elements of the 2d B1 and B2 matrices (6 x 3*N_nodes)
% for nN = 1:N_nodes
%     for ii = 1:6
%         for jj = 1:3
%             B1Matrix(ii, ( CompStruct.Data.DVarNum(ii_d)*(nN - 1) + jj ),1:3,1:3,1:3) = ...
%                 BasicMatrices.Lx(ii,jj)*dxNLMatrix(nN,:,:,:) + ...
%                 BasicMatrices.Ly(ii,jj)*dyNLMatrix(nN,:,:,:);
%             B2Matrix(ii, ( CompStruct.Data.DVarNum(ii_d)*(nN - 1) + jj ),:,:,:) = ...
%                 BasicMatrices.Lz(ii,jj)*BasicMatrices.NLMatrix(nN,:,:,:);
%         end;
%     end;
% end;
% 
% % process all of the elements of the 2d N matrix (3 x 3*N_nodes)
% for nN = 1:N_nodes
%     for ii = 1:3
%         for jj = 1:3
%             NMatrix(ii, ( CompStruct.Data.DVarNum(ii_d)*(nN - 1) + jj ),:,:,:) = ...
%                 BasicMatrices.E3(ii,jj)*BasicMatrices.NLMatrix(nN,:,:,:);
%         end;
%     end;
% end;
% 
% % process all of the elements of the 2d 
% % NtRhoNMatrix, etc. matrices (3*N_nodes x 3*N_nodes)
% for aa = 1:(CompStruct.Data.DVarNum(ii_d)*N_nodes)
%     for bb = 1:(CompStruct.Data.DVarNum(ii_d)*N_nodes)
%         
%         % mass matrix term - the kinetic energy
%         for ii = 1:3
%             Rho = squeeze(ElPhysProps.RhoMatrix(:,:,:));
%             Nib = squeeze(NMatrix(ii,bb,:,:,:));
%             Nt_ai = conj(squeeze(NMatrix(ii,aa,:,:,:)));
%             ElMatrices.NtRhoNMatrix(aa,bb,:,:,:) = ...
%                 squeeze(NtRhoNMatrix(aa,bb,:,:,:)) + convn(Nt_ai,convn(Rho,Nib));
%         end;
%         
%         % stiffness matrix terms - the potential energy
%         for ii = 1:6
%             for jj = 1:6
%                 Cij = squeeze(ElPhysProps.CijMatrix(ii,jj,:,:,:));
%                 B1jb = squeeze(B1Matrix(jj,bb,:,:,:));
%                 B2jb = squeeze(B2Matrix(jj,bb,:,:,:));
%                 B1t_ai = conj(squeeze(B1Matrix(ii,aa,:,:,:)));
%                 B2t_ai = conj(squeeze(B2Matrix(ii,aa,:,:,:)));
%                 ElMatrices.B1tCB1Matrix(aa,bb,:,:,:) = ...
%                     squeeze(B1tCB1Matrix(aa,bb,:,:,:)) + convn(B1t_ai,convn(Cij,B1jb));
%                 ElMatrices.B1tCB2Matrix(aa,bb,:,:,:) = ...
%                     squeeze(B1tCB2Matrix(aa,bb,:,:,:)) + convn(B1t_ai,convn(Cij,B2jb));
%                 ElMatrices.B2tCB1Matrix(aa,bb,:,:,:) = ...
%                     squeeze(B2tCB1Matrix(aa,bb,:,:,:)) + convn(B2t_ai,convn(Cij,B1jb));
%                 ElMatrices.B2tCB2Matrix(aa,bb,:,:,:) = ...
%                     squeeze(B2tCB2Matrix(aa,bb,:,:,:)) + convn(B2t_ai,convn(Cij,B2jb));
%             end;
%         end;
%         
%     end;
% end;

end
