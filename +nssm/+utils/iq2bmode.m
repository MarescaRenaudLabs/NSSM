function bmode = iq2bmode(IQ, DR)
    % convert IQ data to BMode image

    I = abs(IQ); % real envelope
    % I = 20 * log10(I ./ max(I, [], [1 2])) + DR; % log compress
    I = 20 * log10(I ./ max(I(:))) + DR; % log compress
    bmode = uint8(255 * I / DR); % 8-bit log-compressed image

end
