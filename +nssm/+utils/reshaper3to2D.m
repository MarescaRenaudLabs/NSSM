function RF = reshaper3to2D(RF)
    % RF = reshaper3to2D(RF)
    %
    % date:    03-04-2023
    % author:  R. Waasdorp
    % ==============================================================================
    RF = reshape(permute(RF, [1 3 2]), [], size(RF, 2), 1);
end
