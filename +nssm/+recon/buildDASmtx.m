function M_lookup = buildDASmtx(P, LUT)
    % M_lookup = buildDASmtx(P, LUT)
    %
    % Construct the sparse DAS matrices for unique combination of angle and
    % half aperture size, to allow fast beamforming by sparse matrix
    % multiplication, described in detail in (Perrot et al. 2021).
    % Code is inspired by the MUST Toolbox (Garcia 2021).
    %
    % The resulting sparse matrices in M_lookup.DASMTX can be used in the
    % function nssm.recon.reconDASmtx to reconstruct the sound sheet data.
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
    % Output struct M_lookup contains the following fields:
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
    % References:
    %   D. Garcia,
    %       "Make the most of MUST, an open-source MATLAB UltraSound Toolbox,"
    %       2021 IEEE International Ultrasonics Symposium (IUS), 2021, pp. 1-4,
    %       doi: 10.1109/IUS52206.2021.9593605.
    %
    %   V. Perrot, M. Polichetti, F. Varray and D. Garcia,
    %       "So you think you can DAS? A viewpoint on delay-and-sum beamforming,"
    %       Ultrasonics, 111, 2021, p. 106309,
    %       doi: 10.1016/j.ultras.2020.106309.
    %
    % Date:     2024-10-22
    % Author:   B. Heiles, R. Waasdorp
    %
    % =========================================================================

    % get element and reconstruction positions
    xe = P.XElements(:).';
    ze = P.ZElements(:).';
    [x, z] = meshgrid(P.XRecon, P.ZRecon);
    N = numel(x);

    % get number of samples and number of elements
    nl = P.Nz_RF; % number of samples in fast time
    nc = numel(xe); % number of elements

    % check if we have I/Q signals
    isIQ = strcmp(P.DemodMode, 'IQ');

    % Sampling frequency (in Hz)
    Fs = P.Fs;
    Fc = P.DemodFrequency;
    fNumber = P.fNumber;
    t0 = 0; % start time recording

    % Center frequency (in Hz)
    wc = 2 * pi * Fc;

    % Receive distances for fNumber masking
    dxT = x(:) - xe;

    % We only need to compute a DASMTX for a unique combination of angle
    % and halfaperture.
    % Here we find unique combinations HalfApertures and Angles
    [uha, ~, idx_ha] = unique(P.HalfAperture);
    M_lookup.HalfApertures = uha;
    M_lookup.BisIdxToHalfApertureIdx = idx_ha;

    [ua, ~, idx_angles] = unique(P.Angles);
    M_lookup.Angles = ua;
    M_lookup.AngleIdxToAngleIdx = idx_angles;

    % initialize the DASMTX as cell array, for unique combinations
    % of halfAperture and angle.
    Nbisectors = numel(M_lookup.HalfApertures);
    Nangles = numel(M_lookup.Angles);
    M_lookup.DASMTX = cell(Nangles, Nbisectors);

    if Nbisectors > 2; nssm.utils.progressbar_ui(0, Nbisectors); end

    for k_bisector = 1:Nbisectors
        for k_angle = 1:Nangles
            % total travel time for TX/RX combination
            tau_tmp = LUT.TX(:, k_angle, k_bisector) + LUT.RX + P.TimeToPeak;

            % Reshape travel time in matrix of [Nz*Nx, Nelements]
            tau = reshape(tau_tmp, [], nc);

            % Corresponding time indices
            idxt = (tau - t0(:)') * Fs + 1;
            idxt = double(idxt); % in case tau is in single precision

            % Check if the time indices are within the range of the RF data
            I = idxt >= 1 & idxt <= nl - 1;

            % Aperture using the f-number:
            if fNumber > 0
                Iaperture = abs(dxT) <= (z(:) / 2 / fNumber);
                I = I & Iaperture;
            end

            % subscripts to linear indices
            idx = idxt + (0:nc - 1) * nl;
            idx = idx(I);

            % weights (aperture size correction)
            W = ones(N, 1);

            idxf = floor(idx);
            idx = idxf - idx;

            %Build DAS matrix for Linear interpolation
            [i, ~] = find(I);
            j = [idxf; idxf + 1];
            s = [idx + 1; -idx];
            s = s .* repmat(W(i), 2, 1);

            if isIQ % upmixing if data is IQ demodulated
                tau = tau(I);
                s = s .* exp(1i * wc * [tau; tau]);
            end

            % ORIGINAL:
            M = sparse([i; i], j, s, numel(x), nl * nc);

            % and put in DASMTX lookup table
            M_lookup.DASMTX{k_angle, k_bisector} = M;

        end
        if Nbisectors > 2; nssm.utils.progressbar_ui(k_bisector, Nbisectors); end

    end
end
