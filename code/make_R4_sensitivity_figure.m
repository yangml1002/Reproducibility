clear; clc; close all;
script_dir=fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir=pwd; end
cd(script_dir);
result_dir=fullfile(script_dir,'..','data','reviewer4_remaining');
out_dir=fullfile(script_dir,'..','figures');
if exist(out_dir,'dir')~=7, mkdir(out_dir); end
f=fullfile(result_dir,'r4_sensitivity_summary.csv');
if exist(f,'file')~=2, error('Missing %s',f); end
T=readtable(f,'VariableNamingRule','preserve');
T.Properties.VariableNames{2}='case_name';
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');
fig=figure('Color','w');
set(fig,'Units','centimeters','Position',[2 2 22 14.5]);
tl=tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
groups={ ...
    'M15_sigma0','M15_sigma1','M15_lambda_c', ...
    'M15_noise','M13_sparse_pin_count','M15_gain'};
titles={ ...
    '(a) Absolute threshold $\sigma_0$', ...
    '(b) Relative threshold $\sigma_1$', ...
    '(c) Active-rate parameter $\lambda_c$', ...
    '(d) Noise intensity', ...
    '(e) Pinned nodes $p$ (sparse)', ...
    '(f) Pinning gain $k$'};
baseline_values=[5e-6 1e-4 .75 .2 4 26];
for ig=1:6
    nexttile; hold on;
    R=T(strcmp(string(T.group),groups{ig}),:);
    [x,ord]=sort(R.parameter_value);
    R=R(ord,:);
    [~,ib]=min(abs(x-baseline_values(ig)));
    Nu=R.updates_mean/R.updates_mean(ib);
    J=R.J_mean/R.J_mean(ib);
    E=R.E_mean/R.E_mean(ib);
    if ig<=2
        semilogx(x,Nu,'o-','LineWidth',1.35,'MarkerSize',5,'MarkerFaceColor','w');
        semilogx(x,J,'s--','LineWidth',1.35,'MarkerSize',5,'MarkerFaceColor','w');
        semilogx(x,E,'d-.','LineWidth',1.35,'MarkerSize',5,'MarkerFaceColor','w');
    else
        plot(x,Nu,'o-','LineWidth',1.35,'MarkerSize',5,'MarkerFaceColor','w');
        plot(x,J,'s--','LineWidth',1.35,'MarkerSize',5,'MarkerFaceColor','w');
        plot(x,E,'d-.','LineWidth',1.35,'MarkerSize',5,'MarkerFaceColor','w');
    end
    yline(1,':','LineWidth',.9);
    grid on; box on;
    title(titles{ig},'Interpreter','latex','FontWeight','normal','FontSize',11.5);
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
            xlabel('$p$','Interpreter','latex'); xticks(x);
        case 6
            xlabel('$k$','Interpreter','latex'); xticks(x);
    end
    ylabel('Ratio to reference setting');
    ymin=min([Nu;J;E]); ymax=max([Nu;J;E]);
    pad=.08*(ymax-ymin+eps);
    ylim([max(0,ymin-pad) ymax+pad]);
    if ig==1
        legend('$N_u/N_{u,0}$','$J_\varrho/J_{\varrho,0}$','$J_u/J_{u,0}$', ...
            'Interpreter','latex','Location','best','FontSize',8.5);
    end
end
set(findall(fig,'Type','axes'),'FontSize',10,'LineWidth',1.0);
exportgraphics(fig,fullfile(out_dir,'FigureS3_parameter_sensitivity.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(out_dir,'FigureS3_parameter_sensitivity.png'),'Resolution',300);
close(fig);
fprintf('Created FigureS3_parameter_sensitivity.pdf/png in:\n%s\n',out_dir);
