function [V_qP, V_qSV, V_SH] = V_phase_VTI_exact_RPH(rho,c_main,theta)

%   This routine computes velocities of the three bulk plane waves
% for VTI medium
%
% V_phase_VTI_exact_RPH M-file      
%      V_phase_VTI_exact_RPH, by itself, solves Kelvin-Chritoffel equations for 
%      VTI media for given density, elastic moduli tensor, 
%      frequency and wave vector direction, which is defined by angle theta.
%      It calculates phase velocities of the bulk plane waves 
%      for bulk homogeneous medium using exact solution for VTI from
%      Rock physics handbook by Mavko et.al.
%      NB! For homogeneous anisotropic media, there is no frequency
%      dependence of V_phase
%
%      [Timur Zharnikov SMR v0.1 2011]
%
%function [V_qP, V_qSV, V_SH] = V_phase_VTI_exact_RPH(rho,c_main,theta)
%
%  Inputs -
%       rho - medium density. In kg/m^3.
%       c_main - 1 by 5 array (vector) of VTI parameters, which
%               are conviniently chosen as [c11 c13 c33 c44 c66]. 
%               In GPa.
%       theta - wave propagation direction inclination to VTI symmetry
%               axis. In radians
%
%  Outputs -
%       V_phase = [V_qP, V_qSV, V_SH] - 1x3 vector of phase velocities 
%                   of qP, qSV, and SH waves. In km/s
%       V_group - 1x3 vector of group velocity in km/s. Not available at
%       the moment
%
%  M-files required-
%
% Last Modified by Timur Zharnikov SMR v0.1 2011.10.28

%###############################################################################
%
%   Code for KC_eq_solve
%
%###############################################################################
%===============================================================================
% Initialization
%===============================================================================

%-------------------------------------------------------------------------------
% Calculating exact solution of Kelvin-Christoffel equation for VTI medium 
%   for plane waves at angle \theta inclination to VTI axis 
%   using formulas by Mavko et. al. Rock physics handbook
%-------------------------------------------------------------------------------

M_th = ( ( c_main(1) - c_main(4) )*sin(theta).^2 - ( c_main(3) - c_main(4) )*cos(theta).^2 ).^2 ...
       + ( c_main(2) + c_main(4) )^2*sin(2*theta).^2;

A_V_qP = 1/2*( c_main(1)*sin(theta).^2  + c_main(3)*cos(theta).^2 + c_main(4) + sqrt(M_th) );
A_V_qSV = 1/2*( c_main(1)*sin(theta).^2  + c_main(3)*cos(theta).^2 + c_main(4) - sqrt(M_th) );
A_V_SH = c_main(5)*sin(theta).^2 + c_main(4)*cos(theta).^2;

V_qP = sqrt(1e+9*A_V_qP/rho); %1e+9 factor just translates c_main from GPa into Pa as required by the formula
V_qSV = sqrt(1e+9*A_V_qSV/rho); %1e+9 factor just translates c_main from GPa into Pa as required by the formula
V_SH = sqrt(1e+9*A_V_SH/rho); %1e+9 factor just translates c_main from GPa into Pa as required by the formula

%V_phase = [V_qP, V_qSV, V_SH];

end