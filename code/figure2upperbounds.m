clc;
clear;
close all;
a = 3;
b = 1;
L = 0.25;
deltaH = 0.1;
y0 = 1;
Nplot = 1001;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
figure_dir = fullfile(script_dir, '..', 'figures');
data_dir = fullfile(script_dir, '..', 'data', 'final_AD');
if exist(figure_dir, 'dir') ~= 7
    mkdir(figure_dir);
end
if exist(data_dir, 'dir') ~= 7
    mkdir(data_dir);
end
if a <= b || b <= 0 || L <= 0 || deltaH < 0 || y0 < 0
    error('Parameters must satisfy a>b>0, L>0, deltaH>=0 and y0>=0.');
end
r_left = 0;
r_right = a-b;
for k = 1:100
    r_mid = (r_left+r_right)/2;
    F_mid = r_mid-a+b*exp(r_mid*L);
    if F_mid > 0
        r_right = r_mid;
    else
        r_left = r_mid;
    end
end
r = (r_left+r_right)/2;
root_residual = r-a+b*exp(r*L);
tau = linspace(0,L,Nplot);
homogeneous = y0*exp(-r*tau);
Dsharp = deltaH/(a-b)*(1-exp(-(a-b)*tau));
Dslow = deltaH/r*(1-exp(-r*tau));
Dconstant = deltaH/(a-b)*ones(size(tau));
Bsharp = homogeneous+Dsharp;
Bslow = homogeneous+Dslow;
Bconstant = homogeneous+Dconstant;
if max(Bsharp-Bslow) > 1e-12 || max(Bsharp-Bconstant) > 1e-12
    error('Analytical ordering check failed.');
end
%% Plot settings -- consistent with Figure 1
FIG_POS   = [2 2 47 15];
AXIS_FS   = 16;
LABEL_FS  = 20;
TITLE_FS  = 16;
LEGEND_FS = 17;
AXIS_LW   = 1.0;

c_lemma8 = [0.00 0.28 0.80];      % blue
c_common = [0.88 0.12 0.10];      % red
c_const  = [0.12 0.55 0.20];      % green

fig = figure('Color','w');
set(fig,'Units','centimeters','Position',FIG_POS);

tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

%% ------------------------------------------------------------------------
% (a) Total upper estimates
% -------------------------------------------------------------------------
ax1 = nexttile(tl,1);
hold(ax1,'on');

h_const = plot(ax1,tau,Bconstant,'-.', ...
    'Color',c_const,'LineWidth',2.0);

h_common = plot(ax1,tau,Bslow,'--', ...
    'Color',c_common,'LineWidth',2.4);

h_lemma8 = plot(ax1,tau,Bsharp,'-', ...
    'Color',c_lemma8,'LineWidth',2.2);

hold(ax1,'off');
box(ax1,'on');
grid(ax1,'off');
set(ax1,'FontName','Times New Roman','FontSize',AXIS_FS,'LineWidth',AXIS_LW);
xlim(ax1,[0 L]);

xlabel(ax1,'$\tau=t-\iota$','Interpreter','latex','FontSize',LABEL_FS);
ylabel(ax1,'Upper estimate','Interpreter','latex','FontSize',LABEL_FS);
title(ax1,'(a) Total upper estimates','FontWeight','normal','FontSize',TITLE_FS);

lg1 = legend(ax1,[h_lemma8,h_common,h_const], ...
    {'Lemma 8 (this paper)', ...
     'Common-rate estimate', ...
     'Constant-term estimate'}, ...
    'Location','northeast','FontSize',LEGEND_FS,'Box','on');
set(lg1,'Units','normalized');
legend_position = get(lg1,'Position');
legend_position(2) = 0.70;
set(lg1,'Position',legend_position);

%% Zoomed inset for the close Lemma-8/common-rate curves
% Use the terminal part of the finite interval where their separation is
% largest, while still showing the local trend of both curves.
idx_zoom1 = tau >= 0.235 & tau <= 0.25;

% -------------------------------------------------------------------------
% Inset for panel (a): move further left, use a tighter interval
% -------------------------------------------------------------------------
ax_in1 = axes('Position',[0.095 0.205 0.145 0.255]);

idx_in1 = (tau >= 0.195 & tau <= 0.222);

hold(ax_in1,'on');

plot(ax_in1, ...
    tau(idx_in1), ...
    Bsharp(idx_in1), ...
    '-', ...
    'Color',[0.00 0.28 0.80], ...
    'LineWidth',1.8);

plot(ax_in1, ...
    tau(idx_in1), ...
    Bslow(idx_in1), ...
    '--', ...
    'Color',[0.90 0.15 0.10], ...
    'LineWidth',1.8);

hold(ax_in1,'off');

box(ax_in1,'on');
grid(ax_in1,'off');

xlim(ax_in1,[0.21 0.23]);

ylim(ax_in1,[0.73 0.74]);
yticks(ax_in1,[0.73 0.735 0.74]);   % 可选：让刻度更清楚

set(ax_in1, ...
    'FontName','Times New Roman', ...
    'FontSize',9, ...
    'LineWidth',0.8);

%% ------------------------------------------------------------------------
% (b) Perturbation contributions
% -------------------------------------------------------------------------
ax2 = nexttile(tl,2);
hold(ax2,'on');

hD_const = plot(ax2,tau,Dconstant,'-.', ...
    'Color',c_const,'LineWidth',2.0);

hD_common = plot(ax2,tau,Dslow,'--', ...
    'Color',c_common,'LineWidth',2.4);

hD_lemma8 = plot(ax2,tau,Dsharp,'-', ...
    'Color',c_lemma8,'LineWidth',2.2);

hold(ax2,'off');
box(ax2,'on');
grid(ax2,'off');
set(ax2,'FontName','Times New Roman','FontSize',AXIS_FS,'LineWidth',AXIS_LW);
xlim(ax2,[0 L]);
ylim(ax2,[0 1.12*deltaH/(a-b)]);

xlabel(ax2,'$\tau=t-\iota$','Interpreter','latex','FontSize',LABEL_FS);
ylabel(ax2,'Perturbation contribution','Interpreter','latex','FontSize',LABEL_FS);
title(ax2,'(b) Perturbation contributions','FontWeight','normal','FontSize',TITLE_FS);

lg2 = legend(ax2,[hD_lemma8,hD_common,hD_const], ...
    {'Lemma 8 (this paper)', ...
     'Common-rate estimate', ...
     'Constant-term estimate'}, ...
    'Location','northeast','FontSize',LEGEND_FS,'Box','on');
set(lg2,'Units','normalized');
legend_position = get(lg2,'Position');
legend_position(2) = 0.70;
set(lg2,'Position',legend_position);

%% Zoomed inset for the close Lemma-8/common-rate disturbance terms
idx_zoom2 = tau >= 0.20 & tau <= 0.25;

% -------------------------------------------------------------------------
% Inset for panel (b): move lower, use a tighter interval
% -------------------------------------------------------------------------
ax_in2 = axes('Position',[0.675 0.375 0.145 0.255]);

idx_in2 = (tau >= 0.185 & tau <= 0.245);

hold(ax_in2,'on');

plot(ax_in2, ...
    tau(idx_in2), ...
    Dsharp(idx_in2), ...
    '-', ...
    'Color',[0.00 0.28 0.80], ...
    'LineWidth',1.8);

plot(ax_in2, ...
    tau(idx_in2), ...
    Dslow(idx_in2), ...
    '--', ...
    'Color',[0.90 0.15 0.10], ...
    'LineWidth',1.8);

hold(ax_in2,'off');

box(ax_in2,'on');
grid(ax_in2,'off');

xlim(ax_in2,[0.185 0.245]);

y2_min = min([Dsharp(idx_in2), Dslow(idx_in2)],[],'all');
y2_max = max([Dsharp(idx_in2), Dslow(idx_in2)],[],'all');
margin2 = 0.05*(y2_max-y2_min);

ylim(ax_in2,[y2_min-margin2, y2_max+margin2]);

set(ax_in2, ...
    'FontName','Times New Roman', ...
    'FontSize',9, ...
    'LineWidth',0.8);

%% Hide axes toolbars
axes_handles = findall(fig,'Type','axes');
for iaxes = 1:length(axes_handles)
    try
        toolbar_handle = get(axes_handles(iaxes),'Toolbar');
        set(toolbar_handle,'Visible','off');
    catch
    end
end

set(fig,'PaperPositionMode','auto');
exportgraphics(fig,fullfile(figure_dir,'bound_comparison.pdf'),'ContentType','vector');
data_out = [tau(:),Bsharp(:),Bslow(:),Bconstant(:),Dsharp(:),Dslow(:),Dconstant(:)];
dlmwrite(fullfile(data_dir,'bound_comparison.csv'),data_out,'delimiter',',','precision',16);
fid = fopen(fullfile(data_dir,'bound_comparison_parameters.txt'),'w');
if fid < 0
    error('Cannot write the parameter record.');
end
fprintf(fid,'Finite-interval comparison; 0 <= tau <= L.\n');
fprintf(fid,'a = %.16g\nb = %.16g\nL = %.16g\ndelta_H = %.16g\ny0 = %.16g\n',a,b,L,deltaH,y0);
fprintf(fid,'r = %.16g\nroot residual = %.16g\n',r,root_residual);
fprintf(fid,'CSV columns: tau,Bsharp,Bslow,Bconstant,Dsharp,Dslow,Dconstant\n');
fprintf(fid,'All bounds share the same parameters. No experimental data are used in this figure.\n');
fclose(fid);
fprintf('Halanay bound figure completed: r = %.10f, residual = %.3e\n',r,root_residual);
fprintf('Output: %s\n',fullfile(figure_dir,'bound_comparison.pdf'));
