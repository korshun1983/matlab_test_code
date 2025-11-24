function [Asymptotes] =  ComputeAsymptotesSAFE(CompStruct)

%   Part of the toolbox for solving problems of wave propagation
% in arbitrary anisotropic inhomogeneous waveguides.
% For details see User manual 
% and comments to the main script gen_aniso.m
%
% ComputeAsymptotes M-file      
%      ComputeAsymptotes, by itself, computes various
% asymptotes of dispersion curves, like V_mud, Stoneley,
% anisotropic Stoneley (provided F.Karpfinger and R.Prioul permission),
% low-frequency asymptotes of dipole normal modes, etc.
%
%   [T.Zharnikov, D.Syresin, SMR v0.12_12.2012]
%
% function [Asymptotes] =  ComputeAsymptotes(CompStruct)
%
%  Inputs - 
%       CompStruct - structure containing the information about model;
%
%  Outputs -
%       Asymptotes - structure containing the asymptotes of dispersion
%       curves
%
%  M-files required-
%
% Last Modified by Timur Zharnikov SMR v0.1_01.2012

%###############################################################################
%
%   Code for ComputeAsymptotes
%
%###############################################################################
%===============================================================================
% Initialization
%===============================================================================

% Compute asymptotes, if necessary
% Compute asymptote for mud layer, if it is present in the model
switch CompStruct.Model.mud_domain
    case 0
    otherwise
        %extract mud properties
%        Mud_properties = cell2mat(CompStruct.Model.LayerParam(CompStruct.Model.mud_layer));
        Mud_properties = CompStruct.Model.DomainParam{CompStruct.Model.mud_domain};
        Rho_mud = Mud_properties(1);
        Lambda_mud = Mud_properties(2);
        % compute V_mud
        Asymptotes.V_mud = sqrt(Lambda_mud/Rho_mud);
end;

% Compute asymptotes for the outer formation, if it is TTI
switch CompStruct.Model.DomainType{CompStruct.Data.N_domain}
    case 'HTTI'
        % extract outer formation parameters
        Formation_properties = CompStruct.Model.DomainParam{CompStruct.Data.N_domain};
        Rho = 1e+3*Formation_properties(1);
        C_main = Formation_properties(2:6);
        Theta = Formation_properties(7);
        % compute Christoffel equation solution according to 
        % the Rock Physics Handbook
        [V_qP, V_qSV, V_SH] = CompStruct.Methods.V_phase_VTI_exact_RPH(Rho,C_main,Theta); 
        Asymptotes.V_qP  = 1e-3*V_qP;
        Asymptotes.V_qSV = 1e-3*V_qSV;
        Asymptotes.V_SH  = 1e-3*V_SH;
        % compute the Stoneley wave speed (low-frequency asymptote)
        % for VTI homogeneous formation
        if (CompStruct.Model.mud_domain ~= 0 )
            Asymptotes.V_St  = 1/sqrt( Rho_mud*( 1/Lambda_mud + 1/C_main(5) ) );
        end;
    otherwise
end;

end