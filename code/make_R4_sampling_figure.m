clear; clc; close all;
script_dir=fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir=pwd; end
cd(script_dir);
result_dir=fullfile(script_dir,'..','data','reviewer4_remaining');
out_dir=fullfile(script_dir,'..','figures');
if exist(out_dir,'dir')~=7, mkdir(out_dir); end
T=readtable(fullfile(result_dir,'sampling_consistent_summary.csv'),'VariableNamingRule','preserve');
S=readtable(fullfile(result_dir,'r4_schedule_checks.csv'),'VariableNamingRule','preserve');
hs=T.hs(:);
J=T.J_mean(:); Jlo=T.J_ci_low(:); Jhi=T.J_ci_high(:);
Nu=T.updates_mean(:); Nulo=T.updates_ci_low(:); Nuhi=T.updates_ci_high(:);
E=T.E_mean(:); Elo=T.E_ci_low(:); Ehi=T.E_ci_high(:);
hcert=T.h_certified_limit(1);
rref=strcmp(string(S.case),'reference_complete');
if ~any(rref)
    rref=strcmp(string(S.group),'reference');
end
h_active_max=max(S.max_active_width(rref));
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');
fig=figure('Color','w');
set(fig,'Units','centimeters','Position',[3 3 21 8.4]);
tl=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
errorbar(hs,J,J-Jlo,Jhi-J,'o-','LineWidth',1.5,'MarkerSize',5.5, ...
    'MarkerFaceColor','w','CapSize',6);
xline(hcert,'--','LineWidth',1.2);
yline(.05,':','LineWidth',1.2);
[~,idx5]=min(abs(hs-5e-3));
plot(hs(idx5),J(idx5),'o','MarkerSize',8,'LineWidth',1.4);
text(hs(idx5)*1.15,J(idx5)+.0032, ...
    sprintf('$%.2f\\,h_{\\rm cert}$',hs(idx5)/hcert), ...
    'Interpreter','latex','FontSize',9);
text(hcert*.86,.083,'$h_{\rm cert}$','Interpreter','latex', ...
    'Rotation',90,'HorizontalAlignment','center','FontSize',9);
text(1.15e-5,.052,'$J_\varrho=0.05$','Interpreter','latex','FontSize',9);
set(gca,'XScale','log');
grid on; box on;
xlabel('Sampling period $h_s$ (s)','Interpreter','latex');
ylabel('$J_\varrho$','Interpreter','latex');
title('(a) Certified limit and tested error','FontWeight','normal');
xlim([min(hs)*.8 max(hs)*1.25]); ylim([.04 .13]);
nexttile; hold on;
yyaxis left
errorbar(hs,Nu,Nu-Nulo,Nuhi-Nu,'o-','LineWidth',1.5,'MarkerSize',5.5, ...
    'MarkerFaceColor','w','CapSize',6);
ylabel('Mean global updates $N_u$','Interpreter','latex');
ylim([0 max(Nuhi)*1.12]);
yyaxis right
errorbar(hs,E,E-Elo,Ehi-E,'s--','LineWidth',1.5,'MarkerSize',5.5, ...
    'MarkerFaceColor','w','CapSize',6);
ylabel('$J_u$','Interpreter','latex');
ylim([0 max(Ehi)*1.15]);
set(gca,'XScale','log');
xline(hcert,'--','LineWidth',1.2);
xline(h_active_max,':','LineWidth',1.1);
grid on; box on;
xlabel('Sampling period $h_s$ (s)','Interpreter','latex');
title('(b) Resource trade-off','FontWeight','normal');
yyaxis left
text(h_active_max*.88,.82*max(Nuhi), ...
    sprintf('$h_s>%.4f$ s: no interior samples',h_active_max), ...
    'Interpreter','latex','Rotation',90,'HorizontalAlignment','right','FontSize',8.5);
legend('$N_u$','$J_u$','Interpreter','latex','Location','northeast','FontSize',9);
set(findall(fig,'Type','axes'),'FontSize',11,'LineWidth',1);
exportgraphics(fig,fullfile(out_dir,'sampling_range.pdf'),'ContentType','vector');
close(fig);
fprintf('Created sampling_range.pdf\n');
fprintf('Maximum active width read from schedule = %.8g s\n',h_active_max);
