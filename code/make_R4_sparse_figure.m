clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
cd(script_dir);
result_dir = fullfile(script_dir,'..','data','reviewer4_remaining');
out_dir = fullfile(script_dir,'..','figures');
if exist(result_dir,'dir') ~= 7
    error('Cannot find ../data/reviewer4_remaining.');
end
if exist(out_dir,'dir') ~= 7
    mkdir(out_dir);
end
f1 = fullfile(result_dir,'r4_sparse_summary.csv');
f2 = fullfile(result_dir,'sparse_certified_addon_summary.csv');
if exist(f1,'file') ~= 2, error('Missing %s',f1); end
if exist(f2,'file') ~= 2, error('Missing %s',f2); end
T = readtable(f1,'VariableNamingRule','preserve');
T.Properties.VariableNames{1} = 'case_name';
A = readtable(f2,'VariableNamingRule','preserve');
A.Properties.VariableNames{1} = 'case_name';
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');
fig = figure('Color','w');
set(fig,'Units','centimeters','Position',[2 2 26 8.5]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile; hold on; axis equal; axis off;
M = 10;
theta = linspace(pi/2,pi/2+2*pi,10);
theta(end) = [];
x = zeros(M,1);
y = zeros(M,1);
x(1)=0; y(1)=0;
x(2:10)=1.45*cos(theta(:));
y(2:10)=1.45*sin(theta(:));
for j=2:M
    plot([x(1) x(j)],[y(1) y(j)],'-', ...
        'LineWidth',0.85,'Color',[.65 .65 .65]);
    dx = x(j)-x(1);
    dy = y(j)-y(1);
    quiver(x(1)+.55*dx,y(1)+.55*dy,.13*dx,.13*dy,0, ...
        'Color',[.35 .35 .35], ...
        'LineWidth',.9,'MaxHeadSize',1.0);
    quiver(x(1)+.35*dx,y(1)+.35*dy,-.13*dx,-.13*dy,0, ...
        'Color',[.35 .35 .35], ...
        'LineWidth',.9,'MaxHeadSize',1.0);
end
ring_nodes = [2 10 9 8 7 6 5 4 3 2];
for k=1:(numel(ring_nodes)-1)
    i = ring_nodes(k);
    j = ring_nodes(k+1);
    dx = x(j)-x(i);
    dy = y(j)-y(i);
    plot([x(i) x(j)],[y(i) y(j)],'--', ...
        'LineWidth',1.15,'Color',[.12 .12 .12]);
    quiver(x(i)+.40*dx,y(i)+.40*dy,.22*dx,.22*dy,0, ...
        'Color',[.05 .05 .05], ...
        'LineWidth',1.0,'MaxHeadSize',1.1);
end
scatter(x(2:10),y(2:10),135,'w','filled', ...
    'MarkerEdgeColor','k','LineWidth',1.1);
scatter(x(1),y(1),175,[.82 .82 .82],'filled', ...
    'MarkerEdgeColor','k','LineWidth',1.3);
for i=1:M
    text(x(i),y(i),num2str(i), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontSize',10.5,'FontWeight','bold');
end
xlim([-1.9 1.9]);
ylim([-1.8 1.8]);
title('(a) Sparse directed hub-and-ring graph', ...
    'FontWeight','normal','FontSize',12);
nexttile; hold on;
pvals = [2 4 6];
J=zeros(1,3); Jlo=J; Jhi=J;
E=zeros(1,3); Elo=E; Ehi=E;
Nu=zeros(1,3);
for k=1:3
    cname = sprintf('sparse_degree_p%d',pvals(k));
    r = strcmp(string(T.case_name),cname);
    J(k)=T.J_mean(r);
    Jlo(k)=T.J_ci_low(r);
    Jhi(k)=T.J_ci_high(r);
    E(k)=T.E_mean(r);
    Elo(k)=T.E_ci_low(r);
    Ehi(k)=T.E_ci_high(r);
    Nu(k)=T.updates_mean(r);
end
yyaxis left
errorbar(pvals,J,J-Jlo,Jhi-J,'o-', ...
    'LineWidth',1.45,'MarkerSize',6, ...
    'MarkerFaceColor','w','CapSize',7);
ylabel('$J_\varrho$','Interpreter','latex');
yyaxis right
errorbar(pvals,E,E-Elo,Ehi-E,'s--', ...
    'LineWidth',1.45,'MarkerSize',6, ...
    'MarkerFaceColor','w','CapSize',7);
ylabel('$J_u$','Interpreter','latex');
xlabel('Pinned nodes $p$','Interpreter','latex');
xticks(pvals);
grid on; box on;
title('(b) Degree-first pin-count sensitivity ($k=26$)', ...
    'Interpreter','latex','FontWeight','normal','FontSize',10.5);
yyaxis left
for k=1:3
    text(pvals(k),J(k)+.0023,sprintf('$N_u=%.0f$',Nu(k)), ...
        'Interpreter','latex', ...
        'HorizontalAlignment','center', ...
        'FontSize',8.5);
end
legend('$J_\varrho$','$J_u$', ...
    'Interpreter','latex', ...
    'Location','southwest', ...
    'FontSize',9);
nexttile; hold on;
cases = {'sparse_degree_p6_g30','sparse_leaf_p6_g30'};
cats = {'Degree-first','Leaf-first'};
J=zeros(1,2); Jlo=J; Jhi=J;
E=zeros(1,2); Elo=E; Ehi=E;
Nu=zeros(1,2); eigA=zeros(1,2);
for k=1:2
    r = strcmp(string(A.case_name),cases{k});
    J(k)=A.J_mean(r);
    Jlo(k)=A.J_ci_low(r);
    Jhi(k)=A.J_ci_high(r);
    E(k)=A.E_mean(r);
    Elo(k)=A.E_ci_low(r);
    Ehi(k)=A.E_ci_high(r);
    Nu(k)=A.updates_mean(r);
    eigA(k)=A.active_matrix_maxeig(r);
end
xx = 1:2;
yyaxis left
errorbar(xx,J,J-Jlo,Jhi-J,'o-', ...
    'LineWidth',1.45,'MarkerSize',6.5, ...
    'MarkerFaceColor','w','CapSize',7);
ylabel('$J_\varrho$','Interpreter','latex');
yyaxis right
errorbar(xx,E,E-Elo,Ehi-E,'s--', ...
    'LineWidth',1.45,'MarkerSize',6.5, ...
    'MarkerFaceColor','w','CapSize',7);
ylabel('$J_u$','Interpreter','latex');
xlim([.65 2.35]);
xticks(xx);
xticklabels(cats);
grid on; box on;
title('(c) Certified location comparison ($p=6,\ k=30$)', ...
    'Interpreter','latex','FontWeight','normal','FontSize',10.5);
yyaxis left
text(xx(1),J(1)-0.0018, ...
    sprintf('$N_u=%.0f$',Nu(1)), ...
    'Interpreter','latex', ...
    'HorizontalAlignment','center', ...
    'FontSize',8.5);
text(xx(2),J(2)+0.0018, ...
    sprintf('$N_u=%.0f$',Nu(2)), ...
    'Interpreter','latex', ...
    'HorizontalAlignment','center', ...
    'FontSize',8.5);
text(1.5,0.1415, ...
    sprintf('$\\lambda_{\\max}^{\\rm act}=%.3f\\,/\\,%.3f$', ...
    eigA(1),eigA(2)), ...
    'Interpreter','latex', ...
    'HorizontalAlignment','center', ...
    'FontSize',8.3);
legend('$J_\varrho$','$J_u$', ...
    'Interpreter','latex', ...
    'Location','northeast', ...
    'FontSize',9);
set(findall(fig,'Type','axes'),'FontSize',10.5,'LineWidth',1.0);
exportgraphics(fig, ...
    fullfile(out_dir,'sparse_pinning.pdf'), ...
    'ContentType','vector');
close(fig);
fprintf('Created sparse_pinning.pdf\n');
fprintf('Outer ring direction: 2 -> 10 -> 9 -> ... -> 3 -> 2\n');
