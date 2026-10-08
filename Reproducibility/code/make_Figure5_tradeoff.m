clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
cd(script_dir);
result_dir = fullfile(script_dir,'..','data');
out_dir = fullfile(script_dir,'..','figures');
if exist(result_dir,'dir') ~= 7
    error('Cannot find ../data.');
end
if exist(out_dir,'dir') ~= 7
    mkdir(out_dir);
end
f = fullfile(result_dir,'final_main_summary.csv');
if exist(f,'file') ~= 2
    error('Missing %s',f);
end
T = readtable(f,'VariableNamingRule','preserve');
T.Properties.VariableNames{1} = 'case_name';
kap = [.1 .3 1 3];
cp_cases = { ...
    'CP_grid_x0.1', ...
    'CP_grid_x0.3', ...
    'CP_same_threshold', ...
    'CP_grid_x3.0'};
ai_cases = { ...
    'AIPE_grid_x0.1', ...
    'AIPE_grid_x0.3', ...
    'AIPE_aperiodic_A0.4', ...
    'AIPE_grid_x3.0'};
cpN=zeros(1,4); cpJ=cpN; cpJlo=cpN; cpJhi=cpN; cpE=cpN;
aiN=zeros(1,4); aiJ=aiN; aiJlo=aiN; aiJhi=aiN; aiE=aiN;
for k=1:4
    r = strcmp(string(T.case_name),cp_cases{k});
    cpN(k)   = T.mean(r & strcmp(string(T.metric),'updates'));
    cpJ(k)   = T.mean(r & strcmp(string(T.metric),'J_error'));
    cpJlo(k) = T.ci_low(r & strcmp(string(T.metric),'J_error'));
    cpJhi(k) = T.ci_high(r & strcmp(string(T.metric),'J_error'));
    cpE(k)   = T.mean(r & strcmp(string(T.metric),'E_input'));
    r = strcmp(string(T.case_name),ai_cases{k});
    aiN(k)   = T.mean(r & strcmp(string(T.metric),'updates'));
    aiJ(k)   = T.mean(r & strcmp(string(T.metric),'J_error'));
    aiJlo(k) = T.ci_low(r & strcmp(string(T.metric),'J_error'));
    aiJhi(k) = T.ci_high(r & strcmp(string(T.metric),'J_error'));
    aiE(k)   = T.mean(r & strcmp(string(T.metric),'E_input'));
end
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');
fig = figure('Color','w');
set(fig,'Units','centimeters','Position',[3 3 29 8.2]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
idx_k01 = 1;
idx_k3  = 4;
nexttile; hold on;
errorbar(cpN,cpJ,cpJ-cpJlo,cpJhi-cpJ, ...
    'o-','LineWidth',1.45,'MarkerSize',6, ...
    'MarkerFaceColor','w','CapSize',7);
errorbar(aiN,aiJ,aiJ-aiJlo,aiJhi-aiJ, ...
    's--','LineWidth',1.45,'MarkerSize',6, ...
    'MarkerFaceColor','w','CapSize',7);
yline(.05,':','LineWidth',1.0);
text(cpN(idx_k01)+10,cpJ(idx_k01)-1.4e-4,'$\xi=0.1$', ...
    'Interpreter','latex','FontSize',9);
text(cpN(idx_k3)-10,cpJ(idx_k3)-1.4e-4,'$\xi=3$', ...
    'Interpreter','latex','FontSize',9,'HorizontalAlignment','right');
text(aiN(idx_k01)+10,aiJ(idx_k01)+1.4e-4,'$\xi=0.1$', ...
    'Interpreter','latex','FontSize',9);
text(aiN(idx_k3)-10,aiJ(idx_k3)+1.4e-4,'$\xi=3$', ...
    'Interpreter','latex','FontSize',9,'HorizontalAlignment','right');
box on; grid off;
xlabel('Mean global controller updates $N_u$','Interpreter','latex');
ylabel('$J_\varrho$','Interpreter','latex');
title('(a) Error--update trade-off','Interpreter','latex','FontWeight','normal');
ylim([.0415 .0505]);
legend('CP-PETC','Aperiodic AIPE-TC','$J_\varrho=0.05$', ...
    'Interpreter','latex','Location','north','FontSize',9);
nexttile; hold on;
errorbar(cpE,cpJ,cpJ-cpJlo,cpJhi-cpJ, ...
    'o-','LineWidth',1.45,'MarkerSize',6, ...
    'MarkerFaceColor','w','CapSize',7);
errorbar(aiE,aiJ,aiJ-aiJlo,aiJhi-aiJ, ...
    's--','LineWidth',1.45,'MarkerSize',6, ...
    'MarkerFaceColor','w','CapSize',7);
yline(.05,':','LineWidth',1.0);
text(cpE(idx_k01)-.003,cpJ(idx_k01)-1.5e-4,'$\xi=0.1$', ...
    'Interpreter','latex','FontSize',9,'HorizontalAlignment','right');
text(cpE(idx_k3)+.003,cpJ(idx_k3)-1.5e-4,'$\xi=3$', ...
    'Interpreter','latex','FontSize',9);
text(aiE(idx_k01)-.003,aiJ(idx_k01)+1.5e-4,'$\xi=0.1$', ...
    'Interpreter','latex','FontSize',9,'HorizontalAlignment','right');
text(aiE(idx_k3)+.003,aiJ(idx_k3)+1.5e-4,'$\xi=3$', ...
    'Interpreter','latex','FontSize',9);
box on; grid off;
xlabel('Actuation effort $J_u$','Interpreter','latex');
ylabel('$J_\varrho$','Interpreter','latex');
title('(b) Error--effort trade-off','Interpreter','latex','FontWeight','normal');
ylim([.0415 .0505]);
legend('CP-PETC','Aperiodic AIPE-TC','$J_\varrho=0.05$', ...
    'Interpreter','latex','Location','north','FontSize',9);
nexttile; hold on;
resources = [1.00 1.00; 0.75 22512/30000];
bars = bar(resources','grouped');
bars(1).FaceColor = [0 0.4470 0.7410];
bars(2).FaceColor = [0.8500 0.3250 0.0980];
set(gca,'XTick',1:2,'XTickLabel',{'$R_c$','$N_s/30000$'},'TickLabelInterpreter','latex');
ylabel('Normalized resource','Interpreter','latex');
title('(c) Activation and sensing resources','Interpreter','latex','FontWeight','normal');
ylim([0 1.18]); box on; grid off;
legend('CP-PETC','Aperiodic AIPE-TC','Location','north','FontSize',9);
for ib = 1:2
    for ir = 1:2
        text(bars(ib).XEndPoints(ir),resources(ib,ir)+.025,sprintf('%.4g',resources(ib,ir)), ...
            'HorizontalAlignment','center','FontSize',9);
    end
end
set(findall(fig,'Type','axes'),'FontSize',11,'LineWidth',1.0);
% Historical manuscript filename retained for compatibility;
% the current figure contains three panels.
exportgraphics(fig, ...
    fullfile(out_dir,'threshold_tradeoff_2panel.pdf'), ...
    'ContentType','vector');
close(fig);
fprintf('Created threshold_tradeoff_2panel.pdf\n');
