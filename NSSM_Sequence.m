%% NSSM_Sequence
% This example script shows how to implement a basic NSSM sequence on a
% Verasonics Vantage 256. It requires an RCA transducer defined in the
% computeTrans function from Verasonics.
%
% Define the sequence using the following parameters in a struct P:
%   ProbeName           Set the RCA Probe name that is known to VSX
%                       computeTrans function.
%   TransmitFrequency   Set the transmission frequency in Hz
%   DutyCycle           Set the dutycycle of the transmit waveform
%   HalfCycle           Set the number of halfcycles of the transmit
%                       waveform.
%   ImageVoltage        The voltage to use during imaging. Start low to
%                       avoid collapsing of GVs.
%   Angles              Specify the transmit angles in degrees as a row
%                       vector. The angles will be coherently compounded in
%                       the final image.
%   HalfAperture        Specify the half aperture size, in number of
%                       elements.
%   Bisectors           Specify the bisectors to use. Bisectors are defined
%                       in element indices around which the aperture is
%                       centered. Make sure to select bisectors for which
%                       there are enough elements on both sides to fit the
%                       requested HalfAperture.
%                       Bisectors specified as full number will have an inactive
%                       element in the center, whereas bisectors specified
%                       as index+0.5 have a contiguous aperture.
%                       You can specify multiple bisectors, but only one
%                       bisector per array will be reconstructed in the live
%                       preview.
%   ActiveArrays        Set the active Arrays. Set to 1:2 (i.e. [1,2]) to
%                       transmit X waves on both arrays, sequentially. If
%                       you only need one array, specify the required array
%                       number as 1 or 2.
%   SpeedOfSound        Speed of sound in the medium in m/s.
%   ImagingDepth        Image acquisition and reconstruction depth in
%                       meters.
%
% The sequence will transmit X, / and \ waves for each angle and bisector,
% and show a live preview for one selected bisector.
% After acquisition, the full volume can be reconstructed by closing the
% VSX window, and subsequently run the script NSSM_Sequence_Postprocess.m
%
% Sequence tested with Matlab R2021b and Vantage-4.8.6-2302201700.
%
% Date:     2024-10-22
% Author:   B. Heiles, R. Waasdorp
%
% =========================================================================

clear all % clear all to clear persistent variables in beamformer function
close all;
clc

P.ProbeName = 'RC15gV'; % Name of probe that is known to VSX computeTrans

% Transmit settings
P.TransmitFrequency = 13.6e6; % Transmit frequency in Hz
P.DutyCycle = 0.67; % Dutycycle of the pulse
P.HalfCycle = 2; % number of half cycles to transmit
P.ImageVoltage = 10; % Imaging voltage

P.Angles = [7:1:21];
P.HalfAperture = 20; % half aperture size in number of elements
P.Bisectors = P.HalfAperture + 1:0.5:80 - P.HalfAperture; % step size 1 pitch scanning, 0.5 for half pitch scanning
P.ActiveArrays = 1:2; % only 1 or 2 or 1:2 for volumetric imaging

P.SpeedOfSound = 1480; % speed of sound in m/s
P.ImagingDepth = [1 10] * 1e-3; % Image reconstruction depth in meters

%% Define Transducer
Trans.name = P.ProbeName;
Trans.units = 'mm';
Trans = computeTrans(Trans);

%% Input Validation
UserInputValidation(P, Trans);

P.NumTXRX = numel(P.ActiveArrays) * numel(P.Bisectors) * numel(P.Angles) * 3;
P.maxAcqLength = sqrt(P.ImagingDepth(2) ^ 2 + (Trans.elementLength * 1e-3) ^ 2); % [m]

% Choose demodfrequency close to UserParam.TransmitFrequency
D = [15.625, 12.5, 10.4167, 8.9286];
[~, idx] = min(abs(D * 1e6 - P.TransmitFrequency));
P.DemodFrequency = D(idx) * 1e6;
P.maxAcqLengthWvl = P.maxAcqLength * (P.DemodFrequency / P.SpeedOfSound);
P.maxAcqLengthWvlVSX = P.maxAcqLength * (Trans.frequency * 1e6 / P.SpeedOfSound);
P.Nz_RF = ceil(P.maxAcqLengthWvl * 2 * 4/128) * 128; % fit acquisition length to a multiple of 128 samples

%% Specify Resource Parameters Attributes
Resource.Parameters.numTransmit = Trans.numelements; % number of transmit channels.
Resource.Parameters.numRcvChannels = Trans.numelements; % number of receive channels.
Resource.Parameters.speedOfSound = P.SpeedOfSound; % Speed of sound is in m/s. This will be used to calculate the wavelength in the rest of the script
Resource.Parameters.simulateMode = 0; % Enable simulation with 1/0
Resource.VDAS.dmaTimeout = 2e3; % Timeout for VDAS to wait for data to be received in ms

% allocate RcvBuffer
Resource.RcvBuffer.datatype = 'int16'; % so far only supports int16
Resource.RcvBuffer.rowsPerFrame = numel(P.Bisectors) * numel(P.Angles) * 3 * P.Nz_RF;
Resource.RcvBuffer.colsPerFrame = Resource.Parameters.numRcvChannels; % Will receive on all channels
Resource.RcvBuffer.numFrames = 1; % This is set to the number of transmits for the RCA

%% Create Transmit structures
% Define your Waveform
TW.type = 'parametric'; % Define the waveform type
TW.Parameters = [P.TransmitFrequency .* 1e-6, P.DutyCycle, P.HalfCycle, 1];

% create TX structure
TX = repmat(struct('waveform', 1, ...
    'Apod', zeros(1, Trans.numelements), ...
    'Delay', zeros(1, Trans.numelements)), ...
    1, numel(P.ActiveArrays) * numel(P.Bisectors) * numel(P.Angles) * 3);

% loop over all active arrays, bisectors and angles and define the transmits
for k_array = 1:numel(P.ActiveArrays)
    idx_array = P.ActiveArrays(k_array);
    offset_array = (idx_array > 1) * Trans.numelements / 2;

    for k_bis = 1:numel(P.Bisectors)
        idx_bis = P.Bisectors(k_bis) + 0.5 + offset_array;

        half_aperture = P.HalfAperture;
        idx_aperture_left = floor(idx_bis) - half_aperture + (0:half_aperture - 1);
        idx_aperture_right = ceil(idx_bis) + (0:half_aperture - 1);

        aperture_apod = ones(1, numel(idx_aperture_left));

        for k_angle = 1:numel(P.Angles)
            angle_rad = deg2rad(P.Angles(k_angle));

            delay_left = Trans.ElementPos(1:floor(half_aperture), 1) .* 1e-3 .* sin(angle_rad) ./ P.SpeedOfSound;
            delay_left = (delay_left - min(delay_left)) * Trans.frequency * 1e6; % convert to wvl
            delay_right = flip(delay_left);

            % index variable
            k_tx = 3 * (k_angle - 1) + 3 * numel(P.Angles) * (k_bis - 1) + 3 * numel(P.Angles) * numel(P.Bisectors) * (k_array - 1);

            % 1st transmit is X
            TX(k_tx + 1).Delay(idx_aperture_left) = delay_left;
            TX(k_tx + 1).Delay(idx_aperture_right) = delay_right;
            TX(k_tx + 1).Apod(idx_aperture_left) = aperture_apod;
            TX(k_tx + 1).Apod(idx_aperture_right) = aperture_apod;

            % 2nd transmit is /
            TX(k_tx + 2).Delay(idx_aperture_left) = delay_left;
            TX(k_tx + 2).Apod(idx_aperture_left) = aperture_apod;

            % 2nd transmit is /
            TX(k_tx + 3).Delay(idx_aperture_right) = delay_right;
            TX(k_tx + 3).Apod(idx_aperture_right) = aperture_apod;

        end
    end
end

% set transmit voltage
TPC(1).hv = P.ImageVoltage;

%% Define Receive structures
% TGC
TGC.CntrlPts = 1023 .* ones(1, 8);
TGC.rangeMax = P.ImagingDepth(2) ./ (P.SpeedOfSound / (Trans.frequency * 1e6)); % in wavelengths at Trans.frequency
TGC.Waveform = computeTGCWaveform(TGC);

% setup Receives
Receive = repmat(struct('mode', 0, ...
    'Apod', zeros(1, Trans.numelements), ...
    'startDepth', 0, ...
    'endDepth', P.maxAcqLengthWvlVSX, ...
    'TGC', 1, ...
    'bufnum', 1, ...
    'framenum', 1, ...
    'acqNum', [], ...
    'demodFrequency', P.DemodFrequency / 1e6, ...
    'sampleMode', 'NS200BW'), 1, P.NumTXRX);

% default receive apodization
ApodReceive = [ones(1, Trans.numelements / 2), zeros(1, Trans.numelements / 2)];

for k_frame = 1:Resource.RcvBuffer.numFrames
    for k_array = 1:numel(P.ActiveArrays)
        idx_array = P.ActiveArrays(k_array);

        accum_mode = 0;
        apod_receive = ApodReceive;
        if idx_array == 1; apod_receive = flip(ApodReceive); end
        if idx_array == 2 && k_array == 2; accum_mode = 1; end

        for k_bis = 1:numel(P.Bisectors)

            for k_angle = 1:numel(P.Angles)

                % indexing variable
                offset_acqnum = 3 * (k_angle - 1) + 3 * numel(P.Angles) * (k_bis - 1);
                k_rx = offset_acqnum + 3 * numel(P.Angles) * numel(P.Bisectors) * (k_array - 1) + ...
                    numel(P.ActiveArrays) * numel(P.Bisectors) * numel(P.Angles) * 3 * (k_frame - 1);

                % Receive X
                Receive(k_rx + 1).mode = accum_mode;
                Receive(k_rx + 1).Apod = apod_receive;
                Receive(k_rx + 1).acqNum = 1 + offset_acqnum;
                Receive(k_rx + 1).framenum = k_frame;

                % Receive /
                Receive(k_rx + 2).mode = accum_mode;
                Receive(k_rx + 2).Apod = apod_receive;
                Receive(k_rx + 2).acqNum = 2 + offset_acqnum;
                Receive(k_rx + 2).framenum = k_frame;

                % Receive \
                Receive(k_rx + 3).mode = accum_mode;
                Receive(k_rx + 3).Apod = apod_receive;
                Receive(k_rx + 3).acqNum = 3 + offset_acqnum;
                Receive(k_rx + 3).framenum = k_frame;

            end
        end
    end
end

%% Define External Processing

np = 1;

PN.beamform = np;
Process(np).classname = 'External';
Process(np).method = 'nssm.sequence.beamform';
Process(np).Parameters = {'srcbuffer', 'receive', ...
                              'srcbufnum', 1, ...
                              'srcframenum', -1, ...
                              'dstbuffer', 'none'};
np = np + 1;

PN.plot = np;
Process(np).classname = 'External';
Process(np).method = 'nssm.sequence.plotOrtho';
Process(np).Parameters = {'srcbuffer', 'none', ...
                              'dstbuffer', 'none'};
np = np + 1;

%% Define Sequence Control Objects
PRF = 5e3; % Hz
timeBtwTX = round(1 / PRF * 1e6);

nsc = 1;

SC.timebtwTX = nsc;
SeqControl(nsc).command = 'timeToNextAcq'; % time between acquisitions (each Angle)
SeqControl(nsc).argument = timeBtwTX;
nsc = nsc + 1;

SC.timebtwVolumes = nsc;
SeqControl(nsc).command = 'timeToNextAcq';
SeqControl(nsc).argument = 10 * timeBtwTX;
SeqControl(nsc).condition = 'ignore'; % not time sensitive, dont report
nsc = nsc + 1;

SC.matlab = nsc;
SeqControl(nsc).command = 'returnToMatlab';
nsc = nsc + 1;

SC.jump = nsc;
SeqControl(nsc).command = 'jump'; % Jump back command
SeqControl(nsc).argument = 1;
nsc = nsc + 1;

SC.sync = nsc;
SeqControl(nsc).command = 'sync';
nsc = nsc + 1;

%% Define Events to execute
n = 1;

for k_frame = 1:Resource.RcvBuffer.numFrames
    for k_array = 1:numel(P.ActiveArrays)
        for k_bis = 1:numel(P.Bisectors)
            for k_angle = 1:numel(P.Angles)
                k_tx = 3 * (k_angle - 1) + 3 * numel(P.Angles) * (k_bis - 1) + 3 * numel(P.Angles) * numel(P.Bisectors) * (k_array - 1);
                k_rx = k_tx + numel(P.ActiveArrays) * numel(P.Bisectors) * numel(P.Angles) * 3 * (k_frame - 1);

                Event(n).info = sprintf('F%i Arr%02i Bis%02i Ang%02i - X', k_frame, k_array, k_bis, k_angle);
                Event(n).tx = k_tx + 1;
                Event(n).rcv = k_rx + 1;
                Event(n).recon = 0;
                Event(n).process = 0;
                Event(n).seqControl = SC.timebtwTX;
                n = n + 1;

                Event(n).info = sprintf('F%i Arr%02i Bis%02i Ang%02i - /', k_frame, k_array, k_bis, k_angle);
                Event(n).tx = k_tx + 2;
                Event(n).rcv = k_rx + 2;
                Event(n).recon = 0;
                Event(n).process = 0;
                Event(n).seqControl = SC.timebtwTX;
                n = n + 1;

                Event(n).info = sprintf('F%i Arr%02i Bis%02i Ang%02i - \\', k_frame, k_array, k_bis, k_angle);
                Event(n).tx = k_tx + 3;
                Event(n).rcv = k_rx + 3;
                Event(n).recon = 0;
                Event(n).process = 0;
                Event(n).seqControl = SC.timebtwTX;
                n = n + 1;
            end
        end
    end

    % set different time between TX for last transmit event
    Event(n - 1).seqControl = SC.timebtwVolumes;

    % transfer data
    Event(n).info = 'transfer to host';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'transferToHost';
    nsc = nsc + 1;
    n = n + 1;

    Event(n).info = 'wait for transfer complete';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [nsc, SC.matlab];
    SeqControl(nsc).command = 'waitForTransferComplete';
    SeqControl(nsc).argument = nsc - 1;
    nsc = nsc + 1;
    n = n + 1;

    Event(n).info = 'Reconstruct';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = PN.beamform;
    Event(n).seqControl = 0;
    n = n + 1;

    Event(n).info = 'Plot';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = PN.plot;
    Event(n).seqControl = 0;
    n = n + 1;

    % syncronization events for hardware and software sequencers
    Event(n).info = 'mark transfer processed';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [nsc, SC.matlab];
    SeqControl(nsc).command = 'markTransferProcessed';
    SeqControl(nsc).argument = nsc - 2;
    nsc = nsc + 1;
    n = n + 1;

    Event(n).info = 'sync hardware and software sequencers';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [SC.sync];
    n = n + 1;

end

% jump back to beginning
Event(n).info = 'jump back to start';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = [SC.matlab, SC.jump];
n = n + 1;

%% Add UI to select Bisector for live reconstruction
P.ActiveBisectorIdx = round(numel(P.Bisectors) / 2);

bis_idx_max = numel(P.Bisectors);
bis_range = bis_idx_max - 1;
n_ui = 1;
UI(n_ui).Control = {'UserB2', ...
                        'Style', 'VsSlider', ...
                        'Label', 'Bisector Idx', ...
                        'SliderMinMaxVal', [1, bis_idx_max, P.ActiveBisectorIdx], ...
                        'SliderStep', [1 / bis_range, 10 / bis_range], ...
                        'ValueFormat', '%.0f', ...
                    };
UI(n_ui).Callback = {'nssm.sequence.BisectorIdxSliderCallback(UIValue)'};

%% Save matfile and run VSX
matfilename = 'NSSM_Sequence.mat';
save(matfilename)
filename = matfilename; % set filename for VSX
VSX

return

% =========================================================================
% Local functions
% =========================================================================
function UserInputValidation(P, Trans)

    if min(P.Bisectors - P.HalfAperture) < 1
        error('Selected Half Aperture Size and Bisectors do not fit. Either choose a smaller HalfAperture or increase the smallest bisector');
    end
    if max(P.Bisectors + P.HalfAperture) > Trans.numelements
        error('Selected Half Aperture Size and Bisectors do not fit. Either choose a smaller HalfAperture or increase the largest bisector');
    end

end
