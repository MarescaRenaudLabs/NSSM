function BisectorIdxSliderCallback(UIValue)
    % BisectorSliderCallback - Callback for the Bisector slider to set the active bisector for live reconstruction.
    fprintf('Setting Bisector Index to %.0f\n', UIValue);

    P = evalin('base', 'P');
    P.ActiveBisectorIdx = round(UIValue);
    assignin('base', 'P', P);

    % if VSX if frozen, update bisectorIdx, beamform and update plot
    if evalin('base', 'freeze')
        % make RcvBuffer available in current workspace
        Control.Command = 'copyBuffers';
        runAcq(Control);

        % beamform and update plot
        nssm.sequence.beamform(RcvData{1});
        nssm.sequence.plotOrtho();
    end

end
