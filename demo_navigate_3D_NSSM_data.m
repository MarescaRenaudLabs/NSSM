% Example script to view volumetric reconstructed Sound Sheet data.
%
% Date:     2024-10-22
% Author:   B. Heiles, R. Waasdorp
%
% =========================================================================

clc
clear
close all

% Select a dataset to load

% Load volumetric image data of Figure 2E.
%   wells containing two strains of Escherichia coli bacteria, wild type E.
%   coli and GV-expressing E. coli.
load('data/beamformed_volume_figure_2E.mat')

% Load volumetric image data of Figure 3B.
%   Ultrasound imaging of mARG expression in orthopic tumors, 8 days old
%   tumor in mouse.
% load('data/beamformed_volume_figure_3B.mat')

%% View Array 1
log_compress = @(x) 20 * log10(abs(x) ./ max(abs(x(:))));
img_log_NSSI = log_compress(data.IQ_NSSI);
img_log_SSI = log_compress(data.IQ_SSI);

kbis = 50;

f = figure(1); clf;
f.Position = [100 250 1200 600];
subplot(121)
s_ssi = imagesc(data.XRecon * 1e3, data.ZRecon * 1e3, img_log_SSI(:, :, kbis));
caxis([-60 0])
colorbar
daspect([1 1 1])
colormap bone
t_ssi = title(sprintf('SSM - Array 1 Bisector %i', kbis));

subplot(122);
s_nssi = imagesc(data.XRecon * 1e3, data.ZRecon * 1e3, img_log_NSSI(:, :, kbis));
caxis([-40 0])
colorbar
daspect([1 1 1])
colormap bone
t_nssi = title(sprintf('NSSM - Array 1 Bisector %i', kbis));

%% Animate Array 1
for kbis = 1:size(img_log_NSSI, 3)
    s_ssi.CData = img_log_SSI(:, :, kbis);
    s_nssi.CData = img_log_NSSI(:, :, kbis);

    t_ssi.String = sprintf('SSM - Array 1 Bisector %i', kbis);
    t_nssi.String = sprintf('NSSM - Array 1 Bisector %i', kbis);

    drawnow
    pause(0.1)
end

%% View Array 2
log_compress = @(x) 20 * log10(abs(x) ./ max(abs(x(:))));
img_log_NSSI = log_compress(data.IQ_NSSI);
img_log_SSI = log_compress(data.IQ_SSI);

kbis = 76;

f = figure(2); clf;
f.Position = [100 250 1200 600];
subplot(121)
s_ssi = imagesc(data.XRecon * 1e3, data.ZRecon * 1e3, squeeze(img_log_SSI(:, kbis, :)));
caxis([-60 0])
colorbar
daspect([1 1 1])
colormap bone
t_ssi = title(sprintf('SSM - Array 2 - Bisector %i', kbis));

subplot(122);
s_nssi = imagesc(data.XRecon * 1e3, data.ZRecon * 1e3, squeeze(img_log_NSSI(:, kbis, :)));
caxis([-40 0])
colorbar
daspect([1 1 1])
colormap bone
t_nssi = title(sprintf('NSSM - Array 2 - Bisector %i', kbis));

%% Animate Array 2
for kbis = 1:size(img_log_NSSI, 3)
    s_ssi.CData = squeeze(img_log_SSI(:, kbis, :));
    s_nssi.CData = squeeze(img_log_NSSI(:, kbis, :));

    t_ssi.String = sprintf('SSM - Array 2 - Bisector %i', kbis);
    t_nssi.String = sprintf('NSSM - Array 2 - Bisector %i', kbis);

    drawnow
    pause(0.1)
end
