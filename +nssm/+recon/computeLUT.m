function LUT = computeLUT(P)
    % LUT = computeLUT(P)
    % Computes Lookup table for Sound Sheet Delay-and-Sum Beamforming.
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
    % Outputs struct LUT with the following fields:
    %   - TX      [Nz, Nbisectors, Nangles] array with transmit travel times
    %   - RX      [Nz, Nx, Nelements] array with receive travel times
    %
    % Output LUT struct can be provided to the function nssm.recon.reconDAS
    % to reconstruct the sound sheet RFData.
    %
    % Date:     2024-10-22
    % Author:   B. Heiles, R. Waasdorp
    %
    % =========================================================================

    % extract parameters
    Nx = numel(P.XRecon);
    Nz = numel(P.ZRecon);
    xr = P.XRecon(:).';
    zr = P.ZRecon(:);
    xe = P.XElements;
    ze = P.ZElements;
    pitch = xe(2) - xe(1); % pitch in m
    c0 = P.SpeedOfSound;

    Nelem = numel(xe);
    Nangles = numel(P.Angles);
    Nbisectors = numel(P.HalfAperture);

    % initialize LUT arrays
    LUT.TX = zeros(Nz, Nangles, Nbisectors);
    LUT.RX = zeros(Nz, Nx, Nelem);

    [X, Z] = meshgrid(xr, zr);

    % reshape variables to 3D arrays to make LUTs using implicit expansion
    xe = reshape(xe, 1, 1, []);
    ze = reshape(ze, 1, 1, []);
    angles_rad = reshape(deg2rad(P.Angles), 1, []); % convert angles to radians
    half_aperture_m = reshape(P.HalfAperture, 1, 1, []) * pitch; % convert half aperture in number of elements to meters
    cx = c0 ./ cos(angles_rad); % determine supersonic speed

    % transmit (TX) lookup table
    LUT.TX = (zr + half_aperture_m .* tan(angles_rad)) ./ cx;

    % receive (RX) lookup table
    LUT.RX = 1.0 / c0 * sqrt((X - xe) .^ 2 + (Z - ze) .^ 2);

end
