function RFout = reshaper4to3D(RF)
    % RF = reshaper4to3D(RF)
    %
    % date:    25-09-2023
    % author:  R. Waasdorp
    % ==============================================================================

    [s1, s2, s3, s4] = size(RF);
    ts1 = s1 * s3;
    ts2 = s2;
    ts3 = s4;

    RFout = zeros(ts1, ts2, ts3, 'like', RF);
    for k = 1:s4
        RFout(:, :, k) = reshaper3to2D(RF(:, :, :, k));
    end

end
