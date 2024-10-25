function BFIQ = reconDASmtx(P, M_lookup, RF)
    % BFIQ = reconDASmtx(P, M_lookup, RF)
    % Reconstruct the RFData using Delay-and-Sum Beamforming with the provided sparse DAS matrix.
    % The sparse DAS matrix is computed using the function nssm.recon.buildDASmtx.
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
    % Input M_lookup must be a scruct containing the following fields:
    %   - HalfApertures             Array containing the unique half apertures
    %   - BisIdxToHalfApertureIdx   [Nbisectors] lookup table to convert a bisector
    %                               index to the index in the unique HalfAperture array.
    %   - Angles                    Array containing the unique angles
    %   - AngleIdxToAngleIdx        [Nangles] lookup table to convert angle
    %                               index to the index in the unique Angle array.
    %   - DASMTX                    [Nangles, Nbisectors] cell array containing
    %                               the sparse DAS matrices for all unique
    %                               combinations of angles and half apertures.
    %
    % RF data must be an array of size [Nz_RF*NAngles*NBisectors, NChannels_per_array]
    %
    % Date:     2024-10-22
    % Author:   B. Heiles, R. Waasdorp
    %
    % =========================================================================

    Nz = numel(P.ZRecon);
    Nx = numel(P.XRecon);
    Nangles = numel(P.Angles);
    Nbisectors = numel(P.HalfAperture);

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
            RFIQ = RFreformat;
    end

    % and beamform
    if Nbisectors > 2; nssm.utils.progressbar_ui(0, Nbisectors); end
    for kb = 1:Nbisectors
        for ka = 1:Nangles
            % lookup DASMTX
            kb_lookup = M_lookup.BisIdxToHalfApertureIdx(kb);
            ka_lookup = M_lookup.AngleIdxToAngleIdx(ka);
            M = M_lookup.DASMTX{ka_lookup, kb_lookup};

            % Beamform
            RFIQtmp = RFIQ(:, :, ka, kb);
            BFIQtmp = reshape(M * RFIQtmp(:), Nz, Nx);
            BFIQ(:, :, kb) = BFIQ(:, :, kb) + BFIQtmp;
        end
        if Nbisectors > 2; nssm.utils.progressbar_ui(kb, Nbisectors); end
    end

end
