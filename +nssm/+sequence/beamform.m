function BFData = beamform(varargin)
    % Beamform function for use in the NSSM sequence for live reconstruction.

    % persistent variables
    persistent M_lookup BFC Trans P isInitialized

    if isa(varargin{1}, 'char')
        switch varargin{1}
            case 'clear'
                isInitialized = false;
                return;
        end
    else
        RF = varargin{1};
    end

    P = evalin('base', 'P');
    Trans = evalin('base', 'Trans');

    if isempty(isInitialized) || ~isInitialized
        fprintf('Initializing Sound Sheet Beamformer\n')

        TW = evalin('base', 'TW');
        Receive = evalin('base', 'Receive');

        % setup Beamformer Configuration
        BFC.Angles = P.Angles; % angles in degrees
        BFC.HalfAperture = P.HalfAperture; % half aperture in number of elements

        BFC.Fs = Receive(1).decimSampleRate * 1e6;
        BFC.SpeedOfSound = P.SpeedOfSound;
        BFC.DemodFrequency = Receive(1).demodFrequency * 1e6;
        BFC.lambda = P.SpeedOfSound / P.TransmitFrequency;
        BFC.Nz_RF = P.Nz_RF;
        BFC.DemodMode = 'IQ';
        BFC.fNumber = 1;
        BFC.TimeToPeak = TW(1).peak / (Trans.frequency * 1e6);

        BFC.XElements = Trans.ElementPos(1:Trans.numelements / 2, 1) * 1e-3;
        BFC.ZElements = Trans.ElementPos(1:Trans.numelements / 2, 3) * 1e-3;

        BFC.BisectorLocation = interp1(1:Trans.numelements / 2, BFC.XElements, P.Bisectors);

        BFC.XRecon = BFC.XElements(1):Trans.spacingMm / 1000/2:BFC.XElements(end);
        BFC.ZRecon = P.ImagingDepth(1):BFC.lambda / 2:P.ImagingDepth(2);

        % compute traveltime lookuptable
        LUT = nssm.recon.computeLUT(BFC);

        % construct DAS matrix
        M_lookup = nssm.recon.buildDASmtx(BFC, LUT);
        isInitialized = true;
        assignin('base', 'BFC', BFC);
        assignin('base', 'M_lookup', M_lookup);
    end

    % reconstruct for active arrays
    % Reorder RF based on Trans.Connector
    RF = RF(:, Trans.Connector);

    % Select RF Data for one Bisector
    tmp = reshaper2to3D(RF, P.Nz_RF * 3 * numel(P.Angles));
    RF = tmp(:, :, P.ActiveBisectorIdx);

    % Split RF in linear and nonlinear part
    [RF_nonlin, RF_lin] = nssm.recon.splitRF(RF, BFC.Nz_RF);

    % split in array 1 and array 2
    RFData.nonlin{1} = RF_nonlin(:, 1:Trans.numelements / 2);
    RFData.nonlin{2} = RF_nonlin(:, (1 + Trans.numelements / 2):end);

    RFData.lin{1} = RF_lin(:, 1:Trans.numelements / 2);
    RFData.lin{2} = RF_lin(:, (1 + Trans.numelements / 2):end);

    fnames = fieldnames(RFData);
    for kf = 1:numel(fnames)
        for ka = P.ActiveArrays(:).' % make sure its a row vector
            curfield = fnames{kf};
            BFData.(curfield){ka} = nssm.recon.reconDASmtx(BFC, M_lookup, RFData.(curfield){ka});
        end
    end

    % assign output in base workspace for plotting
    if nargout == 0; assignin('base', 'BFData', BFData); end

end
