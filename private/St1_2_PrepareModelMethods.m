function [InputParam] = St1_2_PrepareModelMethods(InputParam)
      
%   Part of the toolbox for solving problems of wave propagation
% in arbitrary anisotropic inhomogeneous waveguides.
% For details see User manual 
% and comments to the main script gen_aniso.m
%
%   St1_2_PrepareModelMethods.m M-file      
%      St1_2_PrepareModelMethods.m assigns 
% the methods, which will be used to solve
% the particular problem.
% This selection is based on the model and the problem, 
% selected as the input in the problem configuration.
%
% NB! Typically this script will not be modified by the user.
%
%   [T.Zharnikov, D.Syresin, SMR v0.12_12.2012]
%
% function [InputParam] = St1_2_PrepareModelMethods(InputParam)
%
%  Inputs - 
%       InputParam - structure containing the input parameters:
%
%  Outputs -
%       InputParam - structure containing the input parameters, 
%               which is updated with Methods structure 
%               with the definitions of the methods to use
%
%  M-files required-
%
% Last Modified by Timur Zharnikov SMR v0.1_01.2012

%###############################################################################
%
%   Code for St1_2_PrepareModelMethods
%
%###############################################################################
%===============================================================================
% Initialization
%===============================================================================

%===============================================================================
% Assigning the methods, which will be used for further steps
%===============================================================================

InputParam.Methods.St1_3_SetModelUser = @St1_3_SetModelUser_sp_SAFE;

cd('..');
root_path = cd('routines\');
root_path = strcat(root_path,'\');
InputParam.Config.root_path = root_path;

InputParam.Methods.chebdif = @chebdif;

InputParam.Methods.em_tensor_VTI = @em_tensor_VTI;
InputParam.Methods.rot_c_ij = @rot_c_ij;
InputParam.Methods.rot_matrix = @rot_matrix;

%InputParam.Methods.ComputeAsymptotes = @ComputeAsymptotes;
InputParam.Methods.V_phase_VTI_exact_RPH = @V_phase_VTI_exact_RPH;

%InputParam.Methods.chebdif = @chebdif;

cd('Mesh2d v24\');
% InputParam.Methods.BH_mesh = @BH_mesh;
InputParam.Methods.MeshFaces = @meshfaces;
cd('..');

cd(root_path);

solver_path = '';

switch InputParam.Config.ProblemType 
% Choice of the problem to be solved 
%   types of problems: 
% - 'spectrum' - find the spectrum of the waveguide;
% - 'source' - compute the problem with the source;
% - 'ExFun' - compute the excitation function;
% - etc.
    case 'spectrum'
        
        solver_path = strcat(solver_path,InputParam.Config.ProblemType,'\');
        switch InputParam.Config.NumMethod 
            % Choice of the method
            %   types of methods:
            % - 'SM' - spectral method;
            % - 'SAFE' - semianalytical FEM;
            % - 'Riccati' - matrix Riccati approach;
            % - etc. (e.g., route search)
            case 'SM'
                % "this case is onmitted for simplification

            case 'SAFE'
                cd(InputParam.Config.root_path);
                cd('routines\');
                InputParam.Methods.ComputeAsymptotes = @ComputeAsymptotesSAFE;
                cd(root_path);

                solver_path = strcat(solver_path,'\');
                
%                 switch InputParam.Config.EigenVar
%                     %Choice of the eigenvariable
%                     %   types of variables: 
%                     % - 'omega' - frequency;
%                     % - 'k' - wavevector (e.g., necessary for viscoelasticity);
%                     % - etc.
%                     case 'omega'
%                         ;
%                     case 'k'
%                         solver_path = strcat(solver_path,'k','\');
%                         
%                         switch InputParam.Config.Symmetry
%                             %Choice of the symmetry used for the solution
%                             %   types of variables:
%                             % - 'none' - no symmetry, compute full series of basis functions;
%                             % - 'plane0' - plane of mirror symmetry at theta=0;
%                             % - etc.
%                             case 'none'
%                                solver_path = strcat(solver_path,'nosym','\');
                InputParam.Config.solver_path = strcat(root_path,solver_path);

                cd(InputParam.Config.solver_path);

                InputParam.Methods.St1_4_SetModelAdvanced = @St1_4_SetModelAdvanced_sp_SAFE;

                InputParam.Methods.St2_PrepareModel = @St2_PrepareModel_sp_SAFE;
                InputParam.Methods.St2_1_PrepareModelParams = @St2_1_PrepareModelParams_sp_SAFE;
                InputParam.Methods.St2_2_PrepareModelMethods = @St2_2_PrepareModelMethods_sp_SAFE;

                InputParam.Methods.St3_PrepareBasicMatrices   = @St3_PrepareBasicMatrices_sp_SAFE;
                InputParam.Methods.St4_ComputeSolution = @St4_ComputeSolution_sp_SAFE;

                cd(root_path);
 %                           otherwise
 %                       end;
 %                   otherwise
 %                   end
            case 'Riccati'
            otherwise
        end
    case 'source'
    case 'ExFun'
    otherwise
end

end
