function [FullMatrices] = ...
    InsertIC_sp_SAFE(MMatrix,SMOp,SOp,DOp,BasicMatrices,CompStruct,ii_l1,ii_l2,ii_BC)
      
%   Part of the toolbox for solving problems of wave propagation
% in arbitrary anisotropic inhomogeneous waveguides.
% For details see User manual 
% and comments to the main script gen_aniso.m
%
%InsertIC_sp_SM M-file      
%      InsertIC_sp_SM, by itself, inserts rows
%       corresponding to the interface conditions
%       into the spectral method matrix.
%       It assumes that all necessary blocks are available.
%
% NB!!! This implementation is specific for spectrum computation 
% with the spectral method.
%
%   [T.Zharnikov, D.Syresin, SMR v0.12_12.2012]
%
%function [MMatrix,SMOp] = ...
%    InsertIC_sp_SM(MMatrix,SMOp,SOp,DOp,BasicMatrices,CompStruct,ii_l1,ii_l2,ii_BC)
%    
%  Inputs - 
%
%       MMatrix - matrix entering into the left-hand side of the 
%               generalized eigenvalue problem
%
%       SMOp - matrix representation of the governing operator 
%
%       SOp - matrix representation of the stress tensor operator
%
%       DOp - matrix representation of the displacement vector operator
%
%       BasicMatrices - structure containing the basic matrices (differentiation, etc.);
%    
%       CompStuct - structure containing the parameters of the model;
%
%       ii_l1 - number of layer 1;
%
%       ii_l2 - number of layer 2;
%
%       ii_BC - type of the boundary conditions between layer 1 and layer 2
%       (not used now, determined automatically from other user defined source of information)
%
%  Outputs -
%
%       MMatrix - matrix entering into the left-hand side of the 
%               generalized eigenvalue problem after incorporating interface
%               conditions
%
%       SMOp - matrix representation of the governing operator after
%               incorporating interface conditions
%
%  M-files required-
%
% Last Modified by Timur Zharnikov SMR v0.1_01.2012

%###############################################################################
%
%   Code for InsertIC_sp_SM
%
%###############################################################################

%===============================================================================
% Initialization, assembling Computational structure with various parameters
%===============================================================================

%==============================================================================
% Inserting interface conditions into the spectral method matrices
%===============================================================================

%===============================================================================
% Setting conditions at interfaces
%===============================================================================

switch CompStruct.Model.LayerType{ii_D1} % type of the first layer in the pair
    case 'fluid' % fluid
        switch CompStruct.Model.LayerType{ii_D2} % type of the second layer in the pair
            
            case 'fluid' % fluid %fluid-fluid contact
                    %for ideal fluid layers there is just one variable - scalar potential
                    ii_vL1 = 1; %phi
                    for ii_kd = 1:3
                        % fluid - fluid condition for u_r continuity
                        % fluid - fluid condition for s_rr continuity
                    end;
                end;
                
            case 'HTTI' % solid %fluid-solid contact
                
%                 D1Pos = [CompStruct.Data.DVarNum{ii_int}(1):CompStruct.Data.DVarNum{ii_int}(2)];
%                 D2Pos = [CompStruct.Data.DVarNum{ii_int + 1}(1):CompStruct.Data.DVarNum{ii_int + 1}(2)];
                [ D1VecRepRow, D1VecRepCol, D1VecRepV ] = find(BMatrixD1{ii_int});
                [ D2VecRepRow, D2VecRepCol, D2VecRepV ] = find(BMatrixD2{ii_int});
                D1VecRepRow = D1VecRepRow + CompStruct.Data.DVarNum{ii_D2}(1) - 1;
                D1VecRepCol = D1VecRepCol + CompStruct.Data.DVarNum{ii_D1}(1) - 1;
                D2VecRepRow = D2VecRepRow + CompStruct.Data.DVarNum{ii_D1}(1) - 1;
                D2VecRepCol = D2VecRepCol + CompStruct.Data.DVarNum{ii_D2}(1) - 1;
                FullMatrices.PMatrix = FullMatrices.PMatrix + ...
                    sparse(D1VecRepRow, D1VecRepCol, D1VecRepV, ...
                        CompStruct.Data.DVarNum{N_domain}(2),CompStruct.Data.DVarNum{N_domain}(2))
                    + ...
                    sparse(D2VecRepRow, D2VecRepCol, D2VecRepV, ...
                        CompStruct.Data.DVarNum{N_domain}(2),CompStruct.Data.DVarNum{N_domain}(2));
                
            otherwise
        end;
        
    case 'HTTI' % solid
        switch CompStruct.Model.LayerType{ii_l2} % type of the second layer in the pair
            
            case 'fluid' % fluid %solid-fluid contact
                for ii_n = 1:CompStruct.Data.N_harmonics
                    %for ideal fluid layers there is just one variable - scalar potential
                    ii_vL1 = 1; %s_rr
                    CondPosL1 = BasicMatrices.Pos{ii_l1}{ii_n}{ii_vL1}(2); %L1 s_rr terminal radial position, variable 1 - s_rr
                    ii_vL2 = 1; %phi
                    CondPosL2 = BasicMatrices.Pos{ii_l2}{ii_n}{ii_vL2}(1); %L2 initial radial position
                    MMatrix(CondPosL1,CondPosL1) = 0;
                    MMatrix(CondPosL2,CondPosL2) = 0;
                    for ii_kd = 1:3
                        % solid - fluid condition for s_rr continuity
                        SMOp{ii_kd}(CondPosL1,:) = SOp{ii_kd}(CondPosL1,:) - SOp{ii_kd}(CondPosL2,:);
                        % solid - fluid condition for u_r continuity
                        SMOp{ii_kd}(CondPosL2,:) = DOp.ur{ii_kd}(CondPosL1,:) - DOp.ur{ii_kd}(CondPosL2,:);
                    end;
                    ii_vL1 = 2; %s_rth
                    CondPosL1 = BasicMatrices.Pos{ii_l1}{ii_n}{ii_vL1}(2); %L1 s_rth terminal radial position, variable 1 - s_rth
                    MMatrix(CondPosL1,CondPosL1) = 0;
                    for ii_kd = 1:3
                        % fluid - solid condition for s_rth = 0
                        SMOp{ii_kd}(CondPosL1,:) = SOp{ii_kd}(CondPosL1,:);
                    end;
                    ii_vL1 = 3; %s_rz
                    CondPosL1 = BasicMatrices.Pos{ii_l1}{ii_n}{ii_vL1}(2); %L1 s_rz terminal radial position, variable 1 - s_rz
                    MMatrix(CondPosL1,CondPosL1) = 0;
                    for ii_kd = 1:3
                        % fluid - solid condition for s_rz = 0
                        SMOp{ii_kd}(CondPosL1,:) = SOp{ii_kd}(CondPosL1,:);
                    end;
                end;

            case 'HTTI' % solid %solid-solid contact
                switch CompStruct.Model.BCType{ii_BC}

                    case 'SSstiff' % stiff solid-solid contact
                        for ii_n = 1:CompStruct.Data.N_harmonics
                            ii_vL1 = 1; %s_rr
                            CondPosL1 = BasicMatrices.Pos{ii_l1}{ii_n}{ii_vL1}(2); %L1 s_rr terminal radial position, variable 1 - s_rr
                            ii_vL2 = 1; %s_rr
                            CondPosL2 = BasicMatrices.Pos{ii_l2}{ii_n}{ii_vL2}(1); %L2 s_rr initial radial position, variable 1 - s_rr
                            MMatrix(CondPosL1,CondPosL1) = 0;
                            MMatrix(CondPosL2,CondPosL2) = 0;
                            for ii_kd = 1:3
                                % solid - solid condition for s_rr continuity
                                SMOp{ii_kd}(CondPosL1,:) = SOp{ii_kd}(CondPosL1,:) - SOp{ii_kd}(CondPosL2,:);
                                % solid - solid condition for u_r continuity
                                SMOp{ii_kd}(CondPosL2,:) = DOp.ur{ii_kd}(CondPosL1,:) - DOp.ur{ii_kd}(CondPosL2,:);
                            end;
                            ii_vL1 = 2; %s_rth
                            CondPosL1 = BasicMatrices.Pos{ii_l1}{ii_n}{ii_vL1}(2); %L1 s_rth terminal radial position, variable 1 - s_rth
                            ii_vL2 = 2; %s_rth
                            CondPosL2 = BasicMatrices.Pos{ii_l2}{ii_n}{ii_vL2}(1); %L2 s_rth initial radial position, variable 1 - s_rth
                            MMatrix(CondPosL1,CondPosL1) = 0;
                            MMatrix(CondPosL2,CondPosL2) = 0;
                            for ii_kd = 1:3
                                % solid - solid condition for s_rth continuity
                                SMOp{ii_kd}(CondPosL1,:) = SOp{ii_kd}(CondPosL1,:) - SOp{ii_kd}(CondPosL2,:);
                                % solid - solid condition for u_th continuity
                                SMOp{ii_kd}(CondPosL2,:) = DOp.uth{ii_kd}(CondPosL1,:) - DOp.uth{ii_kd}(CondPosL2,:);
                            end;
                            ii_vL1 = 3; %s_rz
                            CondPosL1 = BasicMatrices.Pos{ii_l1}{ii_n}{ii_vL1}(2); %L1 s_rz terminal radial position, variable 1 - s_rz
                            ii_vL2 = 3; %s_rz
                            CondPosL2 = BasicMatrices.Pos{ii_l2}{ii_n}{ii_vL2}(1); %L2 s_rz initial radial position, variable 1 - s_rz
                            MMatrix(CondPosL1,CondPosL1) = 0;
                            MMatrix(CondPosL2,CondPosL2) = 0;
                            for ii_kd = 1:3
                                % solid - solid condition for s_rz continuity
                                SMOp{ii_kd}(CondPosL1,:) = SOp{ii_kd}(CondPosL1,:) - SOp{ii_kd}(CondPosL2,:);
                                % solid - solid condition for u_z continuity
                                SMOp{ii_kd}(CondPosL2,:) = DOp.uz{ii_kd}(CondPosL1,:) - DOp.uz{ii_kd}(CondPosL2,:);
                            end;
                        end;

                    case 'SSslip' % slip solid-solid contact
                        for ii_n = 1:CompStruct.Data.N_harmonics
                            ii_vL1 = 1; %s_rr
                            CondPosL1 = BasicMatrices.Pos{ii_l1}{ii_n}{ii_vL1}(2); %L1 s_rr terminal radial position, variable 1 - s_rr
                            ii_vL2 = 1; %s_rr
                            CondPosL2 = BasicMatrices.Pos{ii_l2}{ii_n}{ii_vL2}(1); %L2 s_rr initial radial position, variable 1 - s_rr
                            MMatrix(CondPosL1,CondPosL1) = 0;
                            MMatrix(CondPosL2,CondPosL2) = 0;
                            for ii_kd = 1:3
                                % solid - solid condition for s_rr continuity
                                SMOp{ii_kd}(CondPosL1,:) = SOp{ii_kd}(CondPosL1,:) - SOp{ii_kd}(CondPosL2,:);
                                % solid - solid condition for u_r continuity
                                SMOp{ii_kd}(CondPosL2,:) = DOp.ur{ii_kd}(CondPosL1,:) - DOp.ur{ii_kd}(CondPosL2,:);
                            end;
                            ii_vL1 = 2; %s_rth
                            CondPosL1 = BasicMatrices.Pos{ii_l1}{ii_n}{ii_vL1}(2); %L1 s_rth terminal radial position, variable 1 - s_rth
                            ii_vL2 = 2; %s_rth
                            CondPosL2 = BasicMatrices.Pos{ii_l2}{ii_n}{ii_vL2}(1); %L2 s_rth initial radial position, variable 1 - s_rth
                            MMatrix(CondPosL1,CondPosL1) = 0;
                            MMatrix(CondPosL2,CondPosL2) = 0;
                            for ii_kd = 1:3
                                % solid - solid condition for s_rth = 0 L1
                                SMOp{ii_kd}(CondPosL1,:) = SOp{ii_kd}(CondPosL1,:);
                                % solid - solid condition for s_rth = 0 L2
                                SMOp{ii_kd}(CondPosL2,:) = SOp{ii_kd}(CondPosL2,:);
                            end;
                            ii_vL1 = 3; %s_rz
                            CondPosL1 = BasicMatrices.Pos{ii_l1}{ii_n}{ii_vL1}(2); %L1 s_rz terminal radial position, variable 1 - s_rz
                            ii_vL2 = 3; %s_rz
                            CondPosL2 = BasicMatrices.Pos{ii_l2}{ii_n}{ii_vL2}(1); %L2 s_rz initial radial position, variable 1 - s_rz
                            MMatrix(CondPosL1,CondPosL1) = 0;
                            MMatrix(CondPosL2,CondPosL2) = 0;
                            for ii_kd = 1:3
                                % solid - solid condition for s_rz = 0 L1
                                SMOp{ii_kd}(CondPosL1,:) = SOp{ii_kd}(CondPosL1,:);
                                % solid - solid condition for s_rz = 0 L2
                                SMOp{ii_kd}(CondPosL2,:) = SOp{ii_kd}(CondPosL2,:);
                            end;
                        end;
                    otherwise
                        %dd
                end;
            otherwise
        end;
    otherwise
end;


end