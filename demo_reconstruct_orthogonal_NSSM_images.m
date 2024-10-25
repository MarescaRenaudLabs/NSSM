% Example script demonstrating how to reconstruct Sound Sheet Imaging data.
%
% Date:     2024-10-22
% Author:   B. Heiles, R. Waasdorp
%
% =========================================================================

clear
close all
clc

% Load RFData of representative planes Figure 2E:
%   wells containing two strains of Escherichia coli bacteria, wild type E.
%   coli and GV-expressing E. coli.
load('data\RFData_planes_figure_2E.mat');

% Load RFData of representative planes Figure 3B:
%   Ultrasound imaging of mARG expression in orthopic tumors, 8 days old
%   tumor in mouse.
load('data\RFData_planes_figure_3B.mat');

% all parameters required for reconstruction are defined in struct P

%% Reconstruct RF to linear and nonlinear images

% Split the RFData in the linear and nonlinear part.
[RF_nonlinear, RF_linear] = nssm.recon.splitRF(RFData, P.Nz_RF);

% split data in array 1 and array 2
RF_nonlinear_arr1 = RF_nonlinear(:, 1:end / 2);
RF_nonlinear_arr2 = RF_nonlinear(:, (end / 2 + 1):end);
RF_linear_arr1 = RF_linear(:, 1:end / 2);
RF_linear_arr2 = RF_linear(:, (end / 2 + 1):end);

% construct Lookup table with traveltimes
LUT = nssm.recon.computeLUT(P);

% reconstruct linear and nonlinear data for both arrays
BFIQ_lin_arr1 = nssm.recon.reconDAS(P, LUT, RF_linear_arr1);
BFIQ_lin_arr2 = nssm.recon.reconDAS(P, LUT, RF_linear_arr2);

BFIQ_nonlin_arr1 = nssm.recon.reconDAS(P, LUT, RF_nonlinear_arr1);
BFIQ_nonlin_arr2 = nssm.recon.reconDAS(P, LUT, RF_nonlinear_arr2);

%% Plot Linear and Nonlinear images side by side for both arrays

log_compress = @(x) 20 * log10(abs(x) / max(abs(x(:))));
DynamicRange = [-40 00]; % dB

% array 1
f = figure(1); clf;
f.Position = [680 100 900 420];

subplot(121)
imagesc(P.XRecon * 1e3, P.ZRecon * 1e3, log_compress(BFIQ_lin_arr1), DynamicRange)
daspect([1 1 1])
colormap bone
colorbar
title('SSM - RCA Array 1')
xlabel('x (mm)')
ylabel('z (mm)')

subplot(122)
imagesc(P.XRecon * 1e3, P.ZRecon * 1e3, log_compress(BFIQ_nonlin_arr1), DynamicRange)
daspect([1 1 1])
colormap bone
colorbar
title('NSSM - RCA Array 1')
xlabel('x (mm)')
ylabel('z (mm)')

% array 2
f = figure(2); clf;
f.Position = [680 520 900 420];

subplot(121)
imagesc(P.XRecon * 1e3, P.ZRecon * 1e3, log_compress(BFIQ_lin_arr2), DynamicRange)
daspect([1 1 1])
colormap bone
colorbar
title('SSM - RCA Array 2')
xlabel('y (mm)')
ylabel('z (mm)')

subplot(122)
imagesc(P.XRecon * 1e3, P.ZRecon * 1e3, log_compress(BFIQ_nonlin_arr2), DynamicRange)
daspect([1 1 1])
colormap bone
colorbar
title('NSSM - RCA Array 2')
xlabel('y (mm)')
ylabel('z (mm)')

return
