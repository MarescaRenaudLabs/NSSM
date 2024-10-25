function BFIQ = reconDAS(P, LUT, RF)
    % BFIQ = reconDAS(P, LUT, RF)
    % Reconstruct the RFData using Delay-and-Sum Beamforming with the provided traveltime lookup table (LUT).
    %
    % Input struct P must contain the following fields:
    %   - ZRecon          [Nz] z-coordinates of the reconstruction grid in meters
    %   - XRecon          [Nx] x-coordinates of the reconstruction grid in meters
    %   - XElements       [Nelements] x-coordinates of the transducer elements in meters
    %   - ZElements       [Nelements] z-coordinates of the transducer elements in meters
    %   - Angles          [Nangles] transmit angles in degrees
    %   - HalfAperture    [Nbisectors] half aperture sizes defined in number of elements
    %   - SpeedOfSound    Medium speed of sound in m/s
    %   - TimeToPeak      Time to peak of the imaging waveform in seconds
    %   - fNumber         F-number to be applied to the beamforming
    %   - Fs              Sampling frequency in Hz
    %   - DemodMode       Demodulation mode ('IQ', 'hilbert', 'none')
    %   - DemodFrequency  Demodulation frequency in Hz
    %   - Nz_RF           Number of samples in one TX/RX event
    %                       (Receive(1).endSample - Receive(1).startSample)
    %
    % Input struct LUT must contain the following fields:
    %   - TX      [Nz, Nbisectors, Nangles] array with transmit travel times
    %   - RX      [Nz, Nx, Nelements] array with receive travel times
    %
    % RF data must be an array of size [Nz_RF*NAngles*NBisectors, NChannels_per_array]
    %
    % Date:     2024-10-22
    % Author:   B. Heiles, R. Waasdorp
    %
    % =========================================================================

    Nangles = numel(P.Angles);
    Nbisectors = numel(P.HalfAperture);
    Nelements = numel(P.XElements);

    Nx = numel(P.XRecon);
    Nz = numel(P.ZRecon);
    zax = P.ZRecon(:);
    xax = P.XRecon(:).';

    % angular frequency for upmixing if data is IQ demodulated
    wc = 2 * pi * P.DemodFrequency;

    % time vector RF
    t_rf = (0:P.Nz_RF - 1).' / P.Fs;

    % initialize output
    BFIQ = complex(zeros(Nz, Nx, Nbisectors));

    % demodulate RF
    RFreformat = nssm.utils.reshaper2to4D(RF, P.Nz_RF, Nangles);
    switch P.DemodMode
        case 'IQ'
            RFIQ = nssm.utils.rf2iq(RFreformat, P.Fs, P.DemodFrequency);
        case 'hilbert'
            RFIQ = hilbert(RFreformat);
        case 'none'
            RFIQ = double(RFreformat);
    end

    % beamform per bisector and angle and sum angles coherently
    if Nbisectors > 2; nssm.utils.progressbar_ui(0, Nbisectors); end
    for kb = 1:Nbisectors
        for ka = 1:Nangles
            for ke = 1:Nelements
                tof_img = LUT.RX(:, :, ke) + LUT.TX(:, ka, kb) + P.TimeToPeak;
                tof_img = reshape(tof_img, [], 1);
                bfiq = interp1(t_rf, RFIQ(:, ke, ka, kb), tof_img, 'linear', 0);

                % upmixing if RF is IQ demodulated
                if strcmp(P.DemodMode, 'IQ')
                    bfiq = bfiq .* exp(1i * wc * tof_img);
                end

                bfiq = reshape(bfiq, Nz, Nx);

                % apply f-number
                f_mask = zax ./ (2 .* abs(P.XElements(ke) - xax)) >= P.fNumber;
                BFIQ(:, :, kb) = BFIQ(:, :, kb) + bfiq .* f_mask;
            end
        end

        if Nbisectors > 2; nssm.utils.progressbar_ui(kb, Nbisectors); end
    end

end
