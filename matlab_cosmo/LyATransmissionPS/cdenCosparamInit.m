%
%***********************************************************
%*  [om_m,om_v,om_bh2,h,an,sigma8] = cdenCosparamInit; *
%***********************************************************
%***********************************************************
% Initiates cosmological parameters.
%
% ARGUMENTS
%
% COMPATIBILITY: Matlab(?), Octave
%
% AUTHOR: Avery Meiksin
%
% HISTORY:
%  10 10 17 Creation date.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [om_m,om_v,om_bh2,h,an,sigma8] = cdenCosparamInit;
global omega;
%disp('Planck 2014');
%om_m = 0.308;
%om_v = 1 - om_m;
%om_bh2 = 0.02214;
%h = 0.678;
%an = 0.9608;
%sigma8 = 0.826;
%
%disp('Planck 2016');
%om_m = 0.315;
%om_v = 1 - om_m;
%om_bh2 = 0.02222;
%h = 0.6731;
%an = 0.9655;
%sigma8 = 0.829;
%
disp('Planck 2018');
om_m = 0.3111;
om_v = 1 - om_m;
om_bh2 = 0.02242;
h = 0.6766;
an = 0.9665;
sigma8 = 0.8102;
omega = [om_m, om_v, 0., om_bh2];
