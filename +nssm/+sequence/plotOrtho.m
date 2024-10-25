function plotOrtho(BFData)

    if ~exist('BFData', 'var'); BFData = evalin('base', 'BFData'); end

    P = evalin('base', 'P');
    BFC = evalin('base', 'BFC');
    NRows = numel(P.ActiveArrays);
    active_bisector = P.Bisectors(P.ActiveBisectorIdx);

    DR_nonlin = [-40 0];
    DR_lin = [-40 0];
    log_compress = @(x) 20 * log10(abs(x) / max(abs(x(:))));

    % persistent variables
    persistent fig s_lin s_nonlin titleh

    % Plot
    if isempty(fig) || ~isvalid(fig)
        % if no existing figure, create new figure

        fig = figure(11);
        fig.Position = [100 300 900 600];
        s_lin = gobjects(2, 1);
        s_nonlin = gobjects(2, 1);
        titleh = gobjects(2, 2);

        % plot for active arrays only
        for ka = P.ActiveArrays(:).'
            subplot(NRows, 2, 2 * (ka - 1) + 1)
            s_lin(ka) = imagesc(BFC.XRecon * 1e3, BFC.ZRecon * 1e3, log_compress(BFData.lin{ka}(:, :, 1)), DR_lin);
            daspect([1 1 1])
            colormap gray
            if ka == 1; xlabel('x (mm)'); end
            if ka == 2; xlabel('y (mm)'); end
            ylabel('z (mm)')
            titleh(ka, 1) = title(sprintf('SSM - Array %i - Bisector %.1f', ka, active_bisector));

            subplot(NRows, 2, 2 * (ka - 1) + 2)
            s_nonlin(ka) = imagesc(BFC.XRecon * 1e3, BFC.ZRecon * 1e3, log_compress(BFData.nonlin{ka}(:, :, 1)), DR_nonlin);
            daspect([1 1 1])
            colormap gray
            if ka == 1; xlabel('x (mm)'); end
            if ka == 2; xlabel('y (mm)'); end
            ylabel('z (mm)')
            titleh(ka, 2) = title(sprintf('NSSM - Array %i - Bisector %.1f', ka, active_bisector));
        end
    else
        % reusing existing figure, update image data
        for ka = P.ActiveArrays(:).'
            s_lin(ka).CData = log_compress(BFData.lin{ka}(:, :, 1));
            s_nonlin(ka).CData = log_compress(BFData.nonlin{ka}(:, :, 1));

            titleh(ka, 1).String = sprintf('SSM - Array %i - Bisector %.1f', ka, active_bisector);
            titleh(ka, 2).String = sprintf('NSSM - Array %i - Bisector %.1f', ka, active_bisector);
        end
    end

end
