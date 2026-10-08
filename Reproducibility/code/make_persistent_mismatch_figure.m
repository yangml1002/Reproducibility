function make_persistent_mismatch_figure(out_dir)
script_dir=fileparts(mfilename('fullpath'));
result_dir=fullfile(script_dir,'..','data');
if nargin<1, out_dir=fullfile(tempdir,'AIPE_TC_reproduction_figures'); end
if ~isfolder(out_dir), mkdir(out_dir); end
B=load(fullfile(result_dir,'final_bias_results.mat'),'case_names','curve_time','curves');
bias_names=cellstr(string(B.case_names));
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');
set(0,'DefaultAxesFontSize',15);
set(0,'DefaultTextFontSize',15);
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
legend('Tail-window mean [95\% CI]','Average-rate certificate $R_{\rm av}$', ...
    'Interpreter','latex','Location','northwest','FontSize',9);
set(fig,'PaperPositionMode','auto');
exportgraphics(fig,fullfile(out_dir,'Figure4_persistent_mismatch.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(out_dir,'Figure4_persistent_mismatch.png'),'Resolution',300);
close(fig);
end
