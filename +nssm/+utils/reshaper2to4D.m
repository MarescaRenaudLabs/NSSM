function RF = reshaper2to4D(RF, lSig, nSig)
    % RF = reshaper2to4D(RF, lSig, nSig)
    % Resulting RF has size [lSig, size(RF,2), nSig, size(RF,1)/(lSig*nSig)]
    %
    % date:    03-04-2023
    % author:  R. Waasdorp
    % ==============================================================================

    [s1, s2] = size(RF);
    RF = reshape(permute(RF, [2 1]), s2, lSig, s1 / lSig, []);
    RF = permute(RF, [2 1 3]);
    RF = reshape(RF, lSig, s2, nSig, []);
end
