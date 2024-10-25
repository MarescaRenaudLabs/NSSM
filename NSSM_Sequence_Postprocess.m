%% NSSM_Sequence_Postprocess
% Run this script after running the NSSM_Sequence or after loading data
% acquired with the sequence. The script will reconstruct all sound sheets,
% and compound both arrays into a 3D volume.
%
% Date:     2024-10-22
% Author:   B. Heiles, R. Waasdorp
%
% =========================================================================

% Rearrange RF Data
RFData = RcvData{1}(:, Trans.Connector);

% split in SSM and NSSM data
[RF_NSSM, RF_SSM] = nssm.recon.splitRF(RFData, P.Nz_RF);

%% Init Beamform Config (BFC)
BFC.Angles = P.Angles; % angles in degrees
BFC.HalfAperture = repelem(P.HalfAperture, numel(P.Bisectors)); % half aperture in number of elements
% correct HalfApertures for whole/half bisectors. Whole bisectors have a
% inactive element in the center, so the HalfAperture is equal the number of
% active elements.
% For half bisectors, the aperture is contiguous, therefore the
% HalfAperture becomes the configured HalfAperture-0.5.
BFC.HalfAperture = BFC.HalfAperture - rem(P.Bisectors, 1);

BFC.Fs = Receive(1).decimSampleRate * 1e6;
BFC.SpeedOfSound = P.SpeedOfSound;
BFC.DemodFrequency = Receive(1).demodFrequency * 1e6;
BFC.lambda = P.SpeedOfSound / P.TransmitFrequency;
BFC.Nz_RF = P.Nz_RF;
BFC.DemodMode = 'IQ';
BFC.fNumber = 1;
BFC.TimeToPeak = TW(1).peak / (Trans.frequency * 1e6);

BFC.XElements = Trans.ElementPos(1:Trans.numelements / 2, 1) * 1e-3;
BFC.ZElements = Trans.ElementPos(1:Trans.numelements / 2, 3) * 1e-3;

BFC.BisectorLocation = interp1(1:Trans.numelements / 2, BFC.XElements, P.Bisectors);

BFC.XRecon = BFC.XElements(1):Trans.spacingMm / 1000/2:BFC.XElements(end);
BFC.ZRecon = P.ImagingDepth(1):BFC.lambda / 2:P.ImagingDepth(2);

% compute traveltime lookuptable
LUT = nssm.recon.computeLUT(BFC);
M_lookup = nssm.recon.buildDASmtx(BFC, LUT);

%% reconstruct per array and mode

% reconstruct with DAS (slow)
% BFIQ_SSM_arr1 = nssm.recon.reconDAS(BFC, LUT, RF_SSM(:,1:end/2));
% BFIQ_SSM_arr2 = nssm.recon.reconDAS(BFC, LUT, RF_SSM(:,(end/2+1):end));
% BFIQ_NSSM_arr1 = nssm.recon.reconDAS(BFC, LUT, RF_NSSM(:,1:end/2));
% BFIQ_NSSM_arr2 = nssm.recon.reconDAS(BFC, LUT, RF_NSSM(:,(end/2+1):end));

% reconstruct with sparse matrix multiplication (faster)
BFIQ_SSM_arr1 = nssm.recon.reconDASmtx(BFC, M_lookup, RF_SSM(:, 1:end / 2));
BFIQ_SSM_arr2 = nssm.recon.reconDASmtx(BFC, M_lookup, RF_SSM(:, (end / 2 + 1):end));
BFIQ_NSSM_arr1 = nssm.recon.reconDASmtx(BFC, M_lookup, RF_NSSM(:, 1:end / 2));
BFIQ_NSSM_arr2 = nssm.recon.reconDASmtx(BFC, M_lookup, RF_NSSM(:, (end / 2 + 1):end));

%% Compound the arrays
[~, idx] = min(abs(BFC.XRecon.' - BFC.BisectorLocation));
Nx = numel(BFC.XRecon);
Nz = numel(BFC.ZRecon);

BFIQ_SSM = zeros(Nz, Nx, Nx);
BFIQ_NSSM = zeros(Nz, Nx, Nx);

% SSM
BFIQ_SSM(:, :, idx) = BFIQ_SSM_arr1;
BFIQ_SSM(:, idx, :) = BFIQ_SSM(:, idx, :) + permute(BFIQ_SSM_arr2, [1 3 2]);

% NSSM
BFIQ_NSSM(:, :, idx) = BFIQ_NSSM_arr1;
BFIQ_NSSM(:, idx, :) = BFIQ_NSSM(:, idx, :) + permute(BFIQ_NSSM_arr2, [1 3 2]);
