function [RF_nonlinear, RF_linear] = splitRF(RF, Nz_RF)
    % [RF_nonlinear, RF_linear] = splitRF(RF, Nz_RF)
    % Split the received RF data into linear and nonlinear components.
    %
    % The nonlinear component is obtained by subtracting the two half
    % amplitude transmits (/ and \) from the X transmit.
    %
    % The linear component consists of only the X transmit.
    %
    % Inputs:
    %   RF      RFData, size [3 * NAngles * NBisectors * Nz_RF, NChannels]
    %           The data must be ordered as follows:
    %                [angle1, bisector1, X
    %                 angle1, bisector1, /
    %                 angle1, bisector1, \
    %                 angle2, bisector1, X
    %                 angle2, bisector1, /
    %                 angle2, bisector1, \
    %                 ...
    %                 angleN, bisectorN, \]
    %
    %   Nz_RF   Number of samples in one TX/RX event
    %               (Receive(1).endSample - Receive(1).startSample)
    %
    % Date:     2024-10-22
    % Author:   B. Heiles, R. Waasdorp
    %
    % =========================================================================

    NTXPerAngle = 3;
    NTXTotal = size(RF, 1) / Nz_RF;

    inds = (0:NTXPerAngle:NTXTotal - NTXPerAngle) * Nz_RF;
    inds_X = inds + (1:Nz_RF).' + 0 * Nz_RF;
    inds_L = inds + (1:Nz_RF).' + 1 * Nz_RF;
    inds_R = inds + (1:Nz_RF).' + 2 * Nz_RF;

    RF_linear = RF(inds_X(:), :);
    RF_nonlinear = RF_linear - RF(inds_L(:), :) - RF(inds_R(:), :);

end
