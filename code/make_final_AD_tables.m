clear; clc;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
cd(script_dir);
result_dir = fullfile(script_dir,'..','data');
out_dir = fullfile(script_dir,'..','tables');
if exist(result_dir,'dir') ~= 7
    error('Cannot find ../data: %s',result_dir);
end
if exist(out_dir,'dir') ~= 7, mkdir(out_dir); end
required={ ...
    'final_main_summary.csv', ...
    'final_main_paired.csv', ...
    'final_main_certificates.csv', ...
    'final_bias_summary.csv', ...
    'final_bias_certificates.csv'};
for kk=1:numel(required)
    f=fullfile(result_dir,required{kk});
    if exist(f,'file')~=2, error('Missing required file: %s',f); end
end
T = readtable(fullfile(result_dir,'final_main_summary.csv'),'VariableNamingRule','preserve');
P = readtable(fullfile(result_dir,'final_main_paired.csv'),'VariableNamingRule','preserve');
C = readtable(fullfile(result_dir,'final_main_certificates.csv'),'VariableNamingRule','preserve');
BT = readtable(fullfile(result_dir,'final_bias_summary.csv'),'VariableNamingRule','preserve');
BC = readtable(fullfile(result_dir,'final_bias_certificates.csv'),'VariableNamingRule','preserve');
T.Properties.VariableNames{1} = 'case_name';
C.Properties.VariableNames{1} = 'case_name';
BT.Properties.VariableNames{1} = 'case_name';
BC.Properties.VariableNames{1} = 'case_name';
fprintf('Final table output folder: %s\n',out_dir);
fid=fopen(fullfile(out_dir,'table_strategy_comparison.tex'),'w');
fprintf(fid,'\\begin{table*}[t]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\caption{Control-strategy comparisons under mechanism-specific and identical thresholds.}\n');
fprintf(fid,'\\label{tab:strategy_final}\n');
fprintf(fid,'\\scriptsize\n');
fprintf(fid,'\\renewcommand{\\arraystretch}{1.18}\n');
fprintf(fid,'\\begin{tabular}{lccccc}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,'Method & $R_c$ & $N_s$ & $N_u$ & $J_\\varrho$ & $J_u$ \\\\\n');
fprintf(fid,'\\midrule\n');
fprintf(fid,'\\multicolumn{6}{l}{\\textit{(a) Mechanism‑specific operating configurations}}\\\\\n');
caseA={'CP_original','AIP_aperiodic_A0.4','AIPE_aperiodic_A0.4'};
labelA={'CP‑PETC','AIP‑SDC, aperiodic','AIPE‑TC, aperiodic'};
for ii=1:numel(caseA)
    c=caseA{ii};
    rAct=strcmp(string(T.case_name),c)&strcmp(string(T.metric),'active_fraction');
    rSen=strcmp(string(T.case_name),c)&strcmp(string(T.metric),'sensing_instants');
    rNu=strcmp(string(T.case_name),c)&strcmp(string(T.metric),'updates');
    rJ=strcmp(string(T.case_name),c)&strcmp(string(T.metric),'J_error');
    rE=strcmp(string(T.case_name),c)&strcmp(string(T.metric),'E_input');
    act=T.mean(rAct); sen=T.mean(rSen);
    nu=T.mean(rNu); nul=T.ci_low(rNu); nuh=T.ci_high(rNu);
    jj=T.mean(rJ); jl=T.ci_low(rJ); jh=T.ci_high(rJ);
    ee=T.mean(rE); el=T.ci_low(rE); eh=T.ci_high(rE);
    if abs(nuh-nul)<1e-10
        nu_cell=sprintf('%.0f',nu);
    else
        nu_cell=sprintf('%.2f [%.2f, %.2f]',nu,nul,nuh);
    end
    fprintf(fid,'%s & %.2f & %.0f & %s & %.6f [%.6f, %.6f] & %.4f [%.4f, %.4f] \\\\\n', ...
        labelA{ii},act,sen,nu_cell,jj,jl,jh,ee,el,eh);
end
fprintf(fid,'\\midrule\n');
fprintf(fid,'\\multicolumn{6}{l}{\\textit{(b) Identical event thresholds $(\\sigma_1,\\sigma_0)=(10^{-4},5\\times10^{-6})$}}\\\\\n');
caseB={'CP_same_threshold','AIPE_aperiodic_A0.4'};
labelB={'CP‑PETC','AIPE‑TC, aperiodic'};
for ii=1:numel(caseB)
    c=caseB{ii};
    rAct=strcmp(string(T.case_name),c)&strcmp(string(T.metric),'active_fraction');
    rSen=strcmp(string(T.case_name),c)&strcmp(string(T.metric),'sensing_instants');
    rNu=strcmp(string(T.case_name),c)&strcmp(string(T.metric),'updates');
    rJ=strcmp(string(T.case_name),c)&strcmp(string(T.metric),'J_error');
    rE=strcmp(string(T.case_name),c)&strcmp(string(T.metric),'E_input');
    act=T.mean(rAct); sen=T.mean(rSen);
    nu=T.mean(rNu); nul=T.ci_low(rNu); nuh=T.ci_high(rNu);
    jj=T.mean(rJ); jl=T.ci_low(rJ); jh=T.ci_high(rJ);
    ee=T.mean(rE); el=T.ci_low(rE); eh=T.ci_high(rE);
    fprintf(fid,'%s & %.2f & %.0f & %.2f [%.2f, %.2f] & %.6f [%.6f, %.6f] & %.4f [%.4f, %.4f] \\\\\n', ...
        labelB{ii},act,sen,nu,nul,nuh,jj,jl,jh,ee,el,eh);
end
fprintf(fid,'\\bottomrule\n');
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\end{table*}\n');
fclose(fid);
fid=fopen(fullfile(out_dir,'table_threshold_grid.tex'),'w');
fprintf(fid,'\\begin{table*}[t]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\caption{Common threshold-grid comparison between CP-PETC and aperiodic AIPE-TC.}\n');
fprintf(fid,'\\label{tab:threshold_grid_final}\n');
fprintf(fid,'\\scriptsize\n');
fprintf(fid,'\\renewcommand{\\arraystretch}{1.15}\n');
fprintf(fid,'\\begin{tabular}{c l c c c}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,'$\\xi$ & Method & $N_u$ & $J_\\varrho$ & $J_u$ \\\\\n');
fprintf(fid,'\\midrule\n');
kap=[.1 .3 1 3];
cp_cases={'CP_grid_x0.1','CP_grid_x0.3','CP_same_threshold','CP_grid_x3.0'};
ai_cases={'AIPE_grid_x0.1','AIPE_grid_x0.3','AIPE_aperiodic_A0.4','AIPE_grid_x3.0'};
for jjk=1:4
    for method=1:2
        if method==1
            c=cp_cases{jjk}; label='CP-PETC';
        else
            c=ai_cases{jjk}; label='AIPE-TC, aperiodic';
        end
        rNu=strcmp(string(T.case_name),c)&strcmp(string(T.metric),'updates');
        rJ=strcmp(string(T.case_name),c)&strcmp(string(T.metric),'J_error');
        rE=strcmp(string(T.case_name),c)&strcmp(string(T.metric),'E_input');
        nu=T.mean(rNu); nul=T.ci_low(rNu); nuh=T.ci_high(rNu);
        er=T.mean(rJ); erl=T.ci_low(rJ); erh=T.ci_high(rJ);
        eu=T.mean(rE); eul=T.ci_low(rE); euh=T.ci_high(rE);
        fprintf(fid,'%.1g & %s & %.2f [%.2f, %.2f] & %.6f [%.6f, %.6f] & %.4f [%.4f, %.4f] \\\\\n', ...
            kap(jjk),label,nu,nul,nuh,er,erl,erh,eu,eul,euh);
    end
    if jjk<4, fprintf(fid,'\\addlinespace[1pt]\n'); end
end
fprintf(fid,'\\bottomrule\n');
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\end{table*}\n');
fclose(fid);
fid=fopen(fullfile(out_dir,'table_persistent_mismatch.tex'),'w');
fprintf(fid,'\\begin{table}[t]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\caption{Persistent-mismatch results and certified ultimate bounds.}\n');
fprintf(fid,'\\label{tab:mismatch_final}\n');
fprintf(fid,'\\small\n');
fprintf(fid,'\\renewcommand{\\arraystretch}{1.15}\n');
fprintf(fid,'\\begin{tabular}{c c c c}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,'$q_{\\rm h}$ & $N_u$ & Tail MSE & $R_{\\rm cert}$ \\\\\n');
fprintf(fid,'\\midrule\n');
qh_main=[0 .02 .05 .10];
for jjq=1:numel(qh_main)
    cname=sprintf('bias_q%.2f',qh_main(jjq));
    rNu=strcmp(string(BT.case_name),cname)&strcmp(string(BT.metric),'updates');
    rTail=strcmp(string(BT.case_name),cname)&strcmp(string(BT.metric),'tail_error');
    rC=strcmp(string(BC.case_name),cname);
    nu=BT.mean(rNu); nul=BT.ci_low(rNu); nuh=BT.ci_high(rNu);
    tail=BT.mean(rTail); taill=BT.ci_low(rTail); tailh=BT.ci_high(rTail);
    rc=BC.R_ultimate(rC);
    if qh_main(jjq)==0
        fprintf(fid,'%.2f & %.2f [%.2f, %.2f] & $\\approx 0$ (numerical) & %.6g \\\\\n', ...
            qh_main(jjq),nu,nul,nuh,rc);
    else
    fprintf(fid,'%.2f & %.2f [%.2f, %.2f] & %s [%s, %s] & %.6g \\\\\n', ...
        qh_main(jjq),nu,nul,nuh,latex_sci(tail),latex_sci(taill),latex_sci(tailh),rc);
    end
end
fprintf(fid,'\\bottomrule\n');
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\end{table}\n');
fclose(fid);
fid=fopen(fullfile(out_dir,'key_paired_differences.csv'),'w');
fprintf(fid,'comparison,metric,mean_difference,ci_low,ci_high,orientation\n');
pairA={'AIP_aperiodic_A0.4','AIPE_aperiodic_A0.4','AIP_minus_AIPE'};
pairB={'CP_original','AIPE_aperiodic_A0.4','CP_original_minus_AIPE'};
pairC={'AIPE_aperiodic_A0.4','CP_same_threshold','AIPE_minus_CP_same'};
pairD={'AIPE_periodic','AIPE_aperiodic_A0.4','AIPE_periodic_minus_aperiodic'};
pairs={pairA,pairB,pairC,pairD};
for pp=1:numel(pairs)
    a=pairs{pp}{1}; b=pairs{pp}{2}; label=pairs{pp}{3};
    for met={'updates','J_error','E_input'}
        m=met{1};
        r=strcmp(string(P.case_a),a)&strcmp(string(P.case_b),b)&strcmp(string(P.metric),m);
        if any(r)
            fprintf(fid,'%s,%s,%.12g,%.12g,%.12g,%s-minus-%s\n', ...
                label,m,P.mean_a_minus_b(r),P.ci_low(r),P.ci_high(r),a,b);
        else
            rr=strcmp(string(P.case_a),b)&strcmp(string(P.case_b),a)&strcmp(string(P.metric),m);
            if any(rr)
                mu=-P.mean_a_minus_b(rr);
                cil=-P.ci_high(rr);
                cih=-P.ci_low(rr);
                fprintf(fid,'%s,%s,%.12g,%.12g,%.12g,%s-minus-%s\n', ...
                    label,m,mu,cil,cih,a,b);
            end
        end
    end
end
fclose(fid);
rCP=strcmp(string(T.case_name),'CP_original')&strcmp(string(T.metric),'updates');
rAI=strcmp(string(T.case_name),'AIPE_aperiodic_A0.4')&strcmp(string(T.metric),'updates');
rAIP=strcmp(string(T.case_name),'AIP_aperiodic_A0.4')&strcmp(string(T.metric),'updates');
rCPsame=strcmp(string(T.case_name),'CP_same_threshold')&strcmp(string(T.metric),'updates');
CP_N=T.mean(rCP); AI_N=T.mean(rAI); AIP_N=T.mean(rAIP); CPsame_N=T.mean(rCPsame);
op_reduction=100*(1-AI_N/CP_N);
aip_reduction=100*(1-AI_N/AIP_N);
same_update_increase=100*(AI_N/CPsame_N-1);
sensing_reduction=100*(1-22512/30000);
qh=[.02 .05 .10]';
tail=zeros(numel(qh),1);
for jj=1:numel(qh)
    cname=sprintf('bias_q%.2f',qh(jj));
    rr=strcmp(string(BT.case_name),cname)&strcmp(string(BT.metric),'tail_error');
    tail(jj)=BT.mean(rr);
end
X=[qh.^2 ones(size(qh))];
coef=X\tail;
pred=X*coef;
R2=1-sum((tail-pred).^2)/sum((tail-mean(tail)).^2);
fid=fopen(fullfile(out_dir,'final_numeric_highlights.txt'),'w');
fprintf(fid,'Final A-D numerical highlights (100 paired Brownian paths)\n');
fprintf(fid,'-------------------------------------------------------\n');
fprintf(fid,'Aperiodic amplitude A_ap = 0.4\n');
fprintf(fid,'Mechanism-specific CP -> AIPE update reduction = %.4f %%\n',op_reduction);
fprintf(fid,'AIP-SDC -> AIPE update reduction under same aperiodic windows = %.4f %%\n',aip_reduction);
fprintf(fid,'Same-threshold AIPE update increase relative to CP = %.4f %%\n',same_update_increase);
fprintf(fid,'Activation/sensing-opportunity reduction = %.4f %%\n',sensing_reduction);
fprintf(fid,'Positive-q_h fit: tail_MSE ~= a*q_h^2+b\n');
fprintf(fid,'a = %.12g, b = %.12g, R^2 = %.9f\n',coef(1),coef(2),R2);
fprintf(fid,'q_h=0.20 is retained only as an additional diagnostic case and is not used in the manuscript fit, figures, or tables.\n');
fprintf(fid,'Main-text mismatch set: q_h = 0, 0.02, 0.05, 0.10\n');
fprintf(fid,'Certified threshold grid in main text: threshold multiplier xi = 0.1, 0.3, 1, 3.\n');
fprintf(fid,'xi=10 is omitted from the principal certified comparison.\n');
fclose(fid);
fprintf('\nFINAL TABLE FILES CREATED:\n');
fprintf('  table_strategy_comparison.tex\n');
fprintf('  table_threshold_grid.tex\n');
fprintf('  table_persistent_mismatch.tex\n');
fprintf('  key_paired_differences.csv\n');
fprintf('  final_numeric_highlights.txt\n');

function s=latex_sci(x)
parts=strsplit(sprintf('%.4e',x),'e');
s=sprintf('$%s\\times10^{%d}$',parts{1},str2double(parts{2}));
end
