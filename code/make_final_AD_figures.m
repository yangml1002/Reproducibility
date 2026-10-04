clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
cd(script_dir);
result_dir = fullfile(script_dir,'..','data','final_AD');
out_dir = fullfile(script_dir,'..','figures');
supp_dir = out_dir;
if exist(result_dir,'dir') ~= 7
    error('Cannot find ../data/final_AD: %s', result_dir);
end
if exist(out_dir,'dir') ~= 7, mkdir(out_dir); end
if exist(supp_dir,'dir') ~= 7, mkdir(supp_dir); end
required = { ...
    'final_main_results.mat', ...
    'final_bias_results.mat', ...
    'final_main_summary.csv', ...
    'final_bias_summary.csv', ...
    'final_bias_certificates.csv', ...
    'final_schedule_boundaries.csv'};
for kk = 1:numel(required)
    f = fullfile(result_dir,required{kk});
    if exist(f,'file') ~= 2
        error('Missing required file: %s',f);
    end
end
fprintf('Final figure input folder: %s\n',result_dir);
fprintf('Final figure output folder: %s\n',out_dir);
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');
set(0,'DefaultAxesFontSize',15);
set(0,'DefaultTextFontSize',15);
set(0,'DefaultLegendFontSize',12);
colors = [ ...
    0.0000 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.4660 0.6740 0.1880;
    0.4940 0.1840 0.5560;
    0.3010 0.7450 0.9330;
    0.6350 0.0780 0.1840;
    0.9290 0.6940 0.1250;
    0.0000 0.2000 0.6000;
    0.5500 0.2700 0.0700;
    0.0000 0.0000 0.0000];
S = load(fullfile(result_dir,'final_main_results.mat'));
B = load(fullfile(result_dir,'final_bias_results.mat'));
if iscell(S.case_names)
    main_names = S.case_names;
else
    main_names = cellstr(S.case_names);
end
if iscell(B.case_names)
    bias_names = B.case_names;
else
    bias_names = cellstr(B.case_names);
end
idx_ai = find(strcmp(main_names,'AIPE_aperiodic_A0.4'),1);
idx_ai_x3 = find(strcmp(main_names,'AIPE_grid_x3.0'),1);
if isempty(idx_ai) || isempty(idx_ai_x3)
    error('Required AIPE cases are missing from final_main_results.mat.');
end
t = S.curve_time(:);
state = S.state_example{idx_ai};
if size(state,2) < 20
    error('state_example for AIPE_aperiodic_A0.4 does not contain 20 state columns.');
end
phi = state(:,1:10);
omega = state(:,11:20);
idxZoom = (t <= 0.10);
fig = figure('Color','w');
set(fig,'Units','centimeters','Position',[2,2,40,15]);
subplot(1,2,1);
hold on;
for i=1:10
    plot(t,phi(:,i),'LineWidth',1.15,'Color',colors(i,:));
end
hold off; box on;
xlabel('Time (s)','Interpreter','latex','FontSize',18);
ylabel('$\varpi_i(t)$','Interpreter','latex','FontSize',20);
title('(a) Phase-like states','FontWeight','normal','FontSize',16);
xlim([0 S.Tsim]);
set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',1.0);
legendText = cell(1,10);
for i=1:10, legendText{i}=['$\varpi_{',num2str(i),'}$']; end
legend(legendText,'Interpreter','latex','NumColumns',2, ...
    'Location','northeast','FontSize',11,'Box','on');
axes('Position',[0.23 0.17 0.13 0.25]);
hold on;
for i=1:10
    plot(t(idxZoom),phi(idxZoom,i),'LineWidth',0.9,'Color',colors(i,:));
end
hold off; box on; xlim([0 0.10]);
phiMin=min(phi(idxZoom,:),[],'all'); phiMax=max(phi(idxZoom,:),[],'all');
phiMargin=.05*(phiMax-phiMin+eps);
ylim([phiMin-phiMargin phiMax+phiMargin]);
set(gca,'FontName','Times New Roman','FontSize',8,'LineWidth',0.8);
subplot(1,2,2);
hold on;
for i=1:10
    plot(t,omega(:,i),'LineWidth',1.15,'Color',colors(i,:));
end
hold off; box on;
xlabel('Time (s)','Interpreter','latex','FontSize',18);
ylabel('$\hat{\varpi}_i(t)$','Interpreter','latex','FontSize',20);
title('(b) Velocity-like states','FontWeight','normal','FontSize',16);
xlim([0 S.Tsim]);
set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',1.0);
legendText = cell(1,10);
for i=1:10, legendText{i}=['$\hat{\varpi}_{',num2str(i),'}$']; end
legend(legendText,'Interpreter','latex','NumColumns',2, ...
    'Location','northeast','FontSize',11,'Box','on');
axes('Position',[0.72 0.17 0.13 0.25]);
hold on;
for i=1:10
    plot(t(idxZoom),omega(idxZoom,i),'LineWidth',0.9,'Color',colors(i,:));
end
hold off; box on; xlim([0 0.10]);
omMin=min(omega(idxZoom,:),[],'all'); omMax=max(omega(idxZoom,:),[],'all');
omMargin=.05*(omMax-omMin+eps);
ylim([omMin-omMargin omMax+omMargin]);
set(gca,'FontName','Times New Roman','FontSize',8,'LineWidth',0.8);
set(fig,'PaperPositionMode','auto');
exportgraphics(fig,fullfile(out_dir,'Figure1_states_aperiodic.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(out_dir,'Figure1_states_aperiodic.png'),'Resolution',300);
close(fig);
hs = 5e-5;
Tplot = 0.70;
event_1 = S.event_example{idx_ai};
event_3 = S.event_example{idx_ai_x3};
event_1 = event_1(event_1 <= Tplot);
event_3 = event_3(event_3 <= Tplot);
if numel(event_1)<2 || numel(event_3)<2
    error('Not enough update events for the inter-update figure.');
end
time_1 = event_1(2:end);
time_3 = event_3(2:end);
ratio_1 = diff(event_1)/hs;
ratio_3 = diff(event_3)/hs;
ymax = max([ratio_1(:);ratio_3(:);1]);
ymax = max(2,1.20*ymax);
fig = figure('Color','w');
set(fig,'Units','centimeters','Position',[3,3,19,8.0]);
subplot(1,2,1); hold on;
for k=1:numel(time_1)
    semilogy([time_1(k) time_1(k)],[1 max(1,ratio_1(k))], ...
        'b-','LineWidth',0.55);
end
semilogy(time_1,ratio_1,'ro','MarkerSize',4.0,'LineWidth',0.9);
hRef1 = semilogy([0 Tplot],[1 1],'r--','LineWidth',1.2);
hold off; box on; grid on;
xlabel('Time (s)','Interpreter','latex');
ylabel('$(t_{k+1}-t_k)/h_s$','Interpreter','latex');
title('(a) $(\sigma_1,\sigma_0)=(10^{-4},5\times10^{-6})$', ...
    'Interpreter','latex','FontWeight','normal','FontSize',13);
xlim([0 Tplot]); ylim([0.8 ymax]);
set(gca,'YScale','log','FontName','Times New Roman','FontSize',12,'LineWidth',1);
legend(hRef1,'One sampling period','Location','northwest', ...
    'Interpreter','latex','FontSize',10,'Box','on');
subplot(1,2,2); hold on;
for k=1:numel(time_3)
    semilogy([time_3(k) time_3(k)],[1 max(1,ratio_3(k))], ...
        'b-','LineWidth',0.55);
end
semilogy(time_3,ratio_3,'ro','MarkerSize',4.0,'LineWidth',0.9);
hRef3 = semilogy([0 Tplot],[1 1],'r--','LineWidth',1.2);
hold off; box on; grid on;
xlabel('Time (s)','Interpreter','latex');
ylabel('$(t_{k+1}-t_k)/h_s$','Interpreter','latex');
title('(b) $(\sigma_1,\sigma_0)=(3\times10^{-4},1.5\times10^{-5})$', ...
    'Interpreter','latex','FontWeight','normal','FontSize',13);
xlim([0 Tplot]); ylim([0.8 ymax]);
set(gca,'YScale','log','FontName','Times New Roman','FontSize',12,'LineWidth',1);
legend(hRef3,'One sampling period','Location','northwest', ...
    'Interpreter','latex','FontSize',10,'Box','on');
set(fig,'PaperPositionMode','auto');
exportgraphics(fig,fullfile(out_dir,'Figure3_interupdate_thresholds.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(out_dir,'Figure3_interupdate_thresholds.png'),'Resolution',300);
close(fig);
Tbs = readtable(fullfile(result_dir,'final_bias_summary.csv'),'VariableNamingRule','preserve');
Tbc = readtable(fullfile(result_dir,'final_bias_certificates.csv'),'VariableNamingRule','preserve');
Tbs.Properties.VariableNames{1} = 'case_name';
Tbc.Properties.VariableNames{1} = 'case_name';
q_main = [0 .02 .05 .10];
case_main = {'bias_q0.00','bias_q0.02','bias_q0.05','bias_q0.10'};
idx_main = zeros(size(q_main));
for j=1:numel(q_main)
    idx_main(j)=find(strcmp(bias_names,case_main{j}),1);
end
fig = figure('Color','w');
set(fig,'Units','centimeters','Position',[4,3,20,8.3]);
subplot(1,2,1); hold on;
y_floor = 1e-6;
y_ceil  = 1e-2;
patch([4 5 5 4],[y_floor y_floor y_ceil y_ceil],[0.92 0.92 0.92], ...
    'EdgeColor','none','FaceAlpha',0.35,'HandleVisibility','off');
styles={'k:','b-','r--','m-.'};
for j=1:numel(idx_main)
    semilogy(B.curve_time,max(B.curves(idx_main(j),:),1e-16),styles{j},'LineWidth',1.8);
end
hold off; box on; grid on;
set(gca,'YScale','log', ...
    'FontName','Times New Roman','FontSize',12,'LineWidth',1);
xlim([0.4 5]);
ylim([y_floor y_ceil]);
xlabel('Time (s)','Interpreter','latex');
ylabel('$\mathrm{E}\|\varrho(t)\|^2$','Interpreter','latex');
title('(a) Post-transient mismatch response','FontWeight','normal');
legend('$q_h=0$','$q_h=0.02$','$q_h=0.05$','$q_h=0.10$', ...
    'Interpreter','latex','Location','southwest','FontSize',10);
text(4.02,1.6e-6,'tail window','FontName','Times New Roman', ...
    'FontSize',9,'Color',[0.30 0.30 0.30]);
q_pos = [.02 .05 .10];
meas = zeros(size(q_pos)); lo=zeros(size(q_pos)); hi=zeros(size(q_pos));
cert = zeros(size(q_pos));
for j=1:numel(q_pos)
    cname=sprintf('bias_q%.2f',q_pos(j));
    rows = strcmp(string(Tbs.case_name),cname) & strcmp(string(Tbs.metric),'tail_error');
    meas(j)=Tbs.mean(rows); lo(j)=Tbs.ci_low(rows); hi(j)=Tbs.ci_high(rows);
    crow = strcmp(string(Tbc.case_name),cname);
    cert(j)=Tbc.R_ultimate(crow);
end
subplot(1,2,2); hold on;
errorbar(q_pos,meas,meas-lo,hi-meas,'bo-','LineWidth',1.4, ...
    'MarkerSize',6,'MarkerFaceColor','w','CapSize',8);
plot(q_pos,cert,'rs--','LineWidth',1.4,'MarkerSize',6,'MarkerFaceColor','w');
hold off; box on; grid on;
set(gca,'YScale','log','FontName','Times New Roman','FontSize',12,'LineWidth',1);
xlabel('Mismatch bound $q_h$','Interpreter','latex');
ylabel('Mean-square error','Interpreter','latex');
title('(b) Measured residual and certified bound','FontWeight','normal');
legend('Tail-window mean [95\% CI]','Certified ultimate bound', ...
    'Interpreter','latex','Location','northwest','FontSize',9);
set(fig,'PaperPositionMode','auto');
exportgraphics(fig,fullfile(out_dir,'Figure4_persistent_mismatch.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(out_dir,'Figure4_persistent_mismatch.png'),'Resolution',300);
close(fig);
q_all = [0 .02 .05 .10 .20];
case_all = {'bias_q0.00','bias_q0.02','bias_q0.05','bias_q0.10','bias_q0.20'};
idx_all=zeros(size(q_all));
for j=1:numel(q_all)
    idx_all(j)=find(strcmp(bias_names,case_all{j}),1);
end
fig=figure('Visible','off','Color','w','Position',[100 100 850 520]); hold on;
styles_all={'k:','b-','r--','m-.','g-'};
for j=1:numel(idx_all)
    semilogy(B.curve_time,max(B.curves(idx_all(j),:),1e-16),styles_all{j},'LineWidth',1.6);
end
hold off; box on; grid on;
xlabel('Time (s)','Interpreter','latex');
ylabel('$\mathrm{E}\|\varrho(t)\|^2$','Interpreter','latex');
legend('$q_h=0$','$q_h=0.02$','$q_h=0.05$','$q_h=0.10$','$q_h=0.20$', ...
    'Interpreter','latex','Location','southwest');
set(gca,'FontName','Times New Roman','FontSize',13);
exportgraphics(fig,fullfile(supp_dir,'FigureS2_mismatch_all_qh.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(supp_dir,'FigureS2_mismatch_all_qh.png'),'Resolution',250);
close(fig);
bounds = dlmread(fullfile(result_dir,'final_schedule_boundaries.csv'),',');
if size(bounds,2)~=3
    error('final_schedule_boundaries.csv should contain [alpha_m,beta_m,alpha_{m+1}].');
end
alpha=bounds(:,1); beta=bounds(:,2); right=bounds(:,3);
active_width=beta-alpha;
rest_width=right-beta;
fig=figure('Visible','off','Color','w','Position',[100 100 900 620]);
subplot(2,1,1); hold on;
for j=1:size(bounds,1)
    if alpha(j)>0.35, break; end
    x1=alpha(j); x2=min(beta(j),0.35);
    if x2>x1
        patch([x1 x2 x2 x1],[0 0 1 1],[.85 .85 .85], ...
            'EdgeColor','none');
    end
end
plot([0 .35],[1.02 1.02],'k-','LineWidth',.01);
hold off; box on;
xlim([0 .35]); ylim([0 1.15]); yticks([0 1]); yticklabels({'rest','active'});
xlabel('Time (s)'); title('(a) Genuinely aperiodic activation schedule','FontWeight','normal');
subplot(2,1,2);
plot(0:numel(active_width)-1,active_width,'bo-','LineWidth',1.2,'MarkerSize',4); hold on;
plot(0:numel(rest_width)-1,rest_width,'rs--','LineWidth',1.2,'MarkerSize',4);
hold off; box on; grid on;
xlabel('Cycle index $m$','Interpreter','latex');
ylabel('Width (s)');
title('(b) Cycle-dependent active and rest widths','FontWeight','normal');
legend('Active width','Rest width','Location','best');
exportgraphics(fig,fullfile(supp_dir,'FigureS1_aperiodic_schedule.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(supp_dir,'FigureS1_aperiodic_schedule.png'),'Resolution',250);
close(fig);
fprintf('\nFINAL FIGURES CREATED:\n');
fprintf('  Figure1_states_aperiodic.pdf/png\n');
fprintf('  Figure3_interupdate_thresholds.pdf/png\n');
fprintf('  Figure4_persistent_mismatch.pdf/png\n');
fprintf('ADDITIONAL REPRODUCIBILITY FIGURES:\n');
fprintf('  FigureS1_aperiodic_schedule.pdf/png\n');
fprintf('  FigureS2_mismatch_all_qh.pdf/png\n');
fprintf('\nNote: Manuscript Fig. 3 uses the separate revised Halanay-bound script.\n');
