clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
cd(script_dir);

result_dir = fullfile(script_dir,'..','data');
out_dir = fullfile(tempdir,'AIPE_TC_optional_sensitivity');

if exist(out_dir,'dir') ~= 7
    mkdir(out_dir);
end

data_file = fullfile(result_dir,'parameter_sensitivity_full.csv');
if exist(data_file,'file') ~= 2
    error('Missing finalized sensitivity file: %s',data_file);
end

T = readtable(data_file,'VariableNamingRule','preserve');

required = { ...
    'parameter','value', ...
    'Nu_mean','Jrho_mean','Ju_mean'};

for k = 1:numel(required)
    assert(ismember(required{k},T.Properties.VariableNames), ...
        'Missing required column: %s',required{k});
end

set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');

groups = {'sigma0','sigma1','lambda_c','noise','gain'};

titles = { ...
    '(a) Absolute threshold $\sigma_0$', ...
    '(b) Relative threshold $\sigma_1$', ...
    '(c) Active-rate parameter $\lambda_c$', ...
    '(d) Noise intensity', ...
    '(e) Pinning gain $k$'};

baseline_values = [5e-6, 1e-4, 0.75, 0.2, 26];

fig = figure('Color','w');
set(fig,'Units','centimeters','Position',[2 2 22 14.5]);

tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

for ig = 1:numel(groups)

    nexttile;
    hold on;

    rows = strcmp(string(T.parameter),groups{ig});
    R = T(rows,:);

    assert(~isempty(R), ...
        'Missing finalized sensitivity group: %s',groups{ig});

    [x,ord] = sort(R.value);
    R = R(ord,:);

    [~,ib] = min(abs(x-baseline_values(ig)));

    assert(abs(x(ib)-baseline_values(ig)) < 1e-12, ...
        'Missing baseline for %s',groups{ig});

    Nu = R.Nu_mean ./ R.Nu_mean(ib);
    Jr = R.Jrho_mean ./ R.Jrho_mean(ib);
    Ju = R.Ju_mean ./ R.Ju_mean(ib);

    if ig <= 2
        semilogx(x,Nu,'o-','LineWidth',1.35,'MarkerSize',5);
        semilogx(x,Jr,'s--','LineWidth',1.35,'MarkerSize',5);
        semilogx(x,Ju,'d-.','LineWidth',1.35,'MarkerSize',5);
    else
        plot(x,Nu,'o-','LineWidth',1.35,'MarkerSize',5);
        plot(x,Jr,'s--','LineWidth',1.35,'MarkerSize',5);
        plot(x,Ju,'d-.','LineWidth',1.35,'MarkerSize',5);
    end

    yline(1,':','LineWidth',0.9);

    grid on;
    box on;

    title(titles{ig}, ...
        'Interpreter','latex', ...
        'FontWeight','normal', ...
        'FontSize',11.5);

    switch ig
        case 1
            xlabel('$\sigma_0$','Interpreter','latex');
        case 2
            xlabel('$\sigma_1$','Interpreter','latex');
        case 3
            xlabel('$\lambda_c$','Interpreter','latex');
        case 4
            xlabel('Noise intensity');
        case 5
            xlabel('$k$','Interpreter','latex');
            xticks(x);
    end

    ylabel('Ratio to reference setting');

    ymin = min([Nu;Jr;Ju]);
    ymax = max([Nu;Jr;Ju]);
    pad = 0.08*(ymax-ymin+eps);
    ylim([max(0,ymin-pad), ymax+pad]);

    if ig == 1
        legend( ...
            '$N_u/N_{u,0}$', ...
            '$J_\varrho/J_{\varrho,0}$', ...
            '$J_u/J_{u,0}$', ...
            'Interpreter','latex', ...
            'Location','best', ...
            'FontSize',8.5);
    end

end

% The sixth tile is intentionally unused because finalized Table 8
% contains five one-at-a-time parameter groups.
nexttile;
axis off;

set(findall(fig,'Type','axes'), ...
    'FontSize',10, ...
    'LineWidth',1.0);

exportgraphics( ...
    fig, ...
    fullfile(out_dir,'FigureS3_parameter_sensitivity.pdf'), ...
    'ContentType','vector');

exportgraphics( ...
    fig, ...
    fullfile(out_dir,'FigureS3_parameter_sensitivity.png'), ...
    'Resolution',300);

close(fig);

fprintf(['Created optional finalized sensitivity curves in:\n', ...
         '%s\n'],out_dir);
