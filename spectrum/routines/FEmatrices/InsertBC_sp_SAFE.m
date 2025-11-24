function [MMatrix,SMOp] = ...
    InsertBC_sp_SM(MMatrix,SMOp,SOp,DOp,BasicMatrices,CompStruct)
      
%   Part of the toolbox for solving problems of wave propagation
% in arbitrary anisotropic inhomogeneous waveguides.
% For details see User manual 
% and comments to the main script gen_aniso.m
%
%InsertBC_sp_SM M-file      
%      InsertBC_sp_SM, by itself, inserts rows
%       corresponding to the boundary conditions
%       into the spectral method matrix.
%       It assumes that all necessary blocks are available.
%
% NB!!! This implementation is specific for spectrum computation 
% with the spectral method.
%
%   [T.Zharnikov, D.Syresin, SMR v0.12_12.2012]
%
%function [MMatrix,SMOp] = ...
%    InsertBC_sp_SM(MMatrix,SMOp,SOp,DOp,BasicMatrices,CompStruct)
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
%  Outputs -
%
%       MMatrix - matrix entering into the left-hand side of the 
%               generalized eigenvalue problem after incorporating boundary
%               conditions
%
%       SMOp - matrix representation of the governing operator after incorporating boundary
%               conditions
%
%  M-files required-
%
% Last Modified by Timur Zharnikov SMR v0.1_01.2012

%###############################################################################
%
%   Code for InsertBC_sp_SM
%
%###############################################################################

%===============================================================================
% Initialization, assembling Computational structure with various parameters
%===============================================================================

%==============================================================================
% Inserting boundary conditions into the spectral method matrices
%===============================================================================

%===============================================================================
% Setting inner boundary conditions 
%===============================================================================

ii_l = 1;
switch CompStruct.Model.BCType{ii_l} % types of BCs: 0 - free, 1 - rigid, 2 - natural (or stiff), 3 - slip
    case 'free' % free BC
        switch CompStruct.Model.LayerType{ii_l} % types of layers: 0 - fluid, 1 - solid
            case 'fluid' % fluid
                for ii_n = 1:CompStruct.Data.N_harmonics
                    ii_vL = 1;
                    CondPos = BasicMatrices.Pos{ii_l}{ii_n}{ii_vL}(1); %writing conditions at the initial radial position
                    MMatrix(CondPos,CondPos) = 0;
                    for ii_kd = 1:3
                        SMOp{ii_kd}(CondPos,:) = SOp{ii_kd}(CondPos,:);
                    end;
                end;
%                 for ii_m = 1:CompStruct.LayerPosFactor(1)*CompStruct.N_modes
%                     row_pos = ( ii_m - 1 )*CompStruct.LayerNr(1) + 1;
%                     switch CompStruct.modes(ii_m)
%                         case 0
%                             SMOp_parts(row_pos,:,:) = Displacement_parts(row_pos,:,:);
%                         otherwise
%                             SMOp_parts(row_pos,:,:) = Stress_parts(row_pos,:,:);
%                     end;
%                     MMatrix(row_pos,row_pos) = 0;
%                 end;
            case 'HTTI' % solid
                for ii_n = 1:CompStruct.Data.N_harmonics
                    for ii_v = 1:CompStruct.Data.LayerVarNum(ii_l)
                        CondPos = BasicMatrices.Pos{ii_l}{ii_n}{ii_v}(1); %writing conditions at the initial radial position
                        MMatrix(CondPos,CondPos) = 0;
                        for ii_kd = 1:3
                            SMOp{ii_kd}(CondPos,:) = SOp{ii_kd}(CondPos,:);
                        end;
                    end;
                end;
            otherwise
        end;
    case 'rigid' % rigid BC
        switch CompStruct.Model.LayerType{1}
            case 'fluid' % fluid
                for ii_n = 1:CompStruct.Data.N_harmonics
                    for ii_v = 1:CompStruct.Data.LayerVarNum(ii_l)
                        CondPos = BasicMatrices.Pos{ii_l}{ii_n}{ii_v}(1); %writing conditions at the initial radial position
                        MMatrix(CondPos,CondPos) = 0;
                        for ii_kd = 1:3
                            SMOp{ii_kd}(CondPos,:) = DOp.ur{ii_kd}(CondPos,:);
                        end;
                    end;
                end;
            case 'HTTI' % solid
                for ii_n = 1:CompStruct.Data.N_harmonics
                    ii_vL = 1;
                    CondPos = BasicMatrices.Pos{ii_l}{ii_n}{ii_vL}(1); %writing conditions at the terminal radial position
                    MMatrix(CondPos,CondPos) = 0;
                    for ii_kd = 1:3
                        SMOp{ii_kd}(CondPos,:) = DOp.ur{ii_kd}(CondPos,:);
                    end;
                    ii_vL = 2;
                    CondPos = BasicMatrices.Pos{ii_l}{ii_n}{ii_vL}(1); %writing conditions at the terminal radial position
                    MMatrix(CondPos,CondPos) = 0;
                    for ii_kd = 1:3
                        SMOp{ii_kd}(CondPos,:) = DOp.uth{ii_kd}(CondPos,:);
                    end;
                    ii_vL = 3;
                    CondPos = BasicMatrices.Pos{ii_l}{ii_n}{ii_vL}(1); %writing conditions at the terminal radial position
                    MMatrix(CondPos,CondPos) = 0;
                    for ii_kd = 1:3
                        SMOp{ii_kd}(CondPos,:) = DOp.uz{ii_kd}(CondPos,:);
                    end;
                end;
            otherwise
        end;
    otherwise
end;

%===============================================================================
% Setting outer boundary conditions 
%===============================================================================

ii_l = CompStruct.Data.N_layers;
switch CompStruct.Model.BCType{ii_l + 1} % types of BCs: 0 - free, 1 - rigid, 2 - natural (or stiff), 3 - slip
    case 'free' % free BC
        switch CompStruct.Model.LayerType{ii_l} % types of layers: 0 - fluid, 1 - solid
            case 'fluid' % fluid
                for ii_n = 1:CompStruct.Data.N_harmonics
                    for ii_v = 1:CompStruct.Data.LayerVarNum(ii_l)
                        CondPos = BasicMatrices.Pos{ii_l}{ii_n}{ii_v}(2); %writing conditions at the terminal radial position
                        MMatrix(CondPos,CondPos) = 0;
                        for ii_kd = 1:3
                            SMOp{ii_kd}(CondPos,:) = SOp{ii_kd}(CondPos,:);
                        end;
                    end;
                end;
            case 'HTTI' % solid
                for ii_n = 1:CompStruct.Data.N_harmonics
                    for ii_v = 1:CompStruct.Data.LayerVarNum(ii_l)
                        CondPos = BasicMatrices.Pos{ii_l}{ii_n}{ii_v}(2); %writing conditions at the terminal radial position
                        MMatrix(CondPos,CondPos) = 0;
                        for ii_kd = 1:3
                            SMOp{ii_kd}(CondPos,:) = SOp{ii_kd}(CondPos,:);
                        end;
                    end;
                end;
            otherwise
        end;
    case 'rigid' % rigid BC
        switch CompStruct.Model.LayerType{ii_l}
            case 'fluid' % fluid
                for ii_n = 1:CompStruct.Data.N_harmonics
                    for ii_v = 1:CompStruct.Data.LayerVarNum(ii_l)
                        CondPos = BasicMatrices.Pos{ii_l}{ii_n}{ii_v}(2); %writing conditions at the terminal radial position
                        MMatrix(CondPos,CondPos) = 0;
                        for ii_kd = 1:3
                            SMOp{ii_kd}(CondPos,:) = DOp.ur{ii_kd}(CondPos,:);
                        end;
                    end;
                end;
            case 'HTTI' % solid
                for ii_n = 1:CompStruct.Data.N_harmonics
                    ii_vL = 1;
                    CondPos = BasicMatrices.Pos{ii_l}{ii_n}{ii_vL}(2); %writing conditions at the terminal radial position
                    MMatrix(CondPos,CondPos) = 0;
                    for ii_kd = 1:3
                        SMOp{ii_kd}(CondPos,:) = DOp.ur{ii_kd}(CondPos,:);
                    end;
                    ii_vL = 2;
                    CondPos = BasicMatrices.Pos{ii_l}{ii_n}{ii_vL}(2); %writing conditions at the terminal radial position
                    MMatrix(CondPos,CondPos) = 0;
                    for ii_kd = 1:3
                        SMOp{ii_kd}(CondPos,:) = DOp.uth{ii_kd}(CondPos,:);
                    end;
                    ii_vL = 3;
                    CondPos = BasicMatrices.Pos{ii_l}{ii_n}{ii_vL}(2); %writing conditions at the terminal radial position
                    MMatrix(CondPos,CondPos) = 0;
                    for ii_kd = 1:3
                        SMOp{ii_kd}(CondPos,:) = DOp.uz{ii_kd}(CondPos,:);
                    end;
                end;
            otherwise
        end;
    otherwise
end;

end