clear; clc;
script_dir=fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir=pwd; end
cd(script_dir);
result_dir=fullfile(script_dir,'..','data','reviewer4_remaining');
out_dir=fullfile(script_dir,'..','tables');
if exist(result_dir,'dir')~=7
    error('Cannot find ../data/reviewer4_remaining.');
end
if exist(out_dir,'dir')~=7, mkdir(out_dir); end
Ts=readtable(fullfile(result_dir,'r4_sparse_summary.csv'),'VariableNamingRule','preserve');
Ts.Properties.VariableNames{1}='case_name';
Ta=readtable(fullfile(result_dir,'sparse_certified_addon_summary.csv'),'VariableNamingRule','preserve');
Ta.Properties.VariableNames{1}='case_name';
Th=readtable(fullfile(result_dir,'sampling_consistent_summary.csv'),'VariableNamingRule','preserve');
Tsen=readtable(fullfile(result_dir,'r4_sensitivity_summary.csv'),'VariableNamingRule','preserve');
Tsen.Properties.VariableNames{2}='case_name';
fid=fopen(fullfile(out_dir,'table_R4_sparse_pinning.tex'),'w');
fprintf(fid,'\\begin{table*}[t]\n\\centering\n\\scriptsize\n');
fprintf(fid,'\\caption{Sparse-topology pinning results.}\n');
fprintf(fid,'\\label{tab:r4_sparse_pinning}\n');
fprintf(fid,'\\renewcommand{\\arraystretch}{1.15}\n');
fprintf(fid,'\\begin{tabular}{lccccc}\n\\toprule\n');
fprintf(fid,'Configuration & $N_u$ & $J_\\varrho$ & $J_u$ & $\\lambda_{\\max}^{\\rm act}$ & Theorem conditions \\\\\n\\midrule\n');
fprintf(fid,'\\multicolumn{6}{l}{\\textit{(a) Degree-first pin-count sensitivity, gain $k=26$}}\\\\\n');
for p=[2 4 6]
    cname=sprintf('sparse_degree_p%d',p);
    r=strcmp(string(Ts.case_name),cname);
    cert=Ts.active_ok(r)&&Ts.rest_ok(r)&&Ts.sampling_ok(r);
    fprintf(fid,'$p=%d$ & %.2f [%.2f, %.2f] & %.6f [%.6f, %.6f] & %.4f [%.4f, %.4f] & %.4f & %s \\\\\n', ...
        p,Ts.updates_mean(r),Ts.updates_ci_low(r),Ts.updates_ci_high(r), ...
        Ts.J_mean(r),Ts.J_ci_low(r),Ts.J_ci_high(r), ...
        Ts.E_mean(r),Ts.E_ci_low(r),Ts.E_ci_high(r), ...
        Ts.active_matrix_maxeig(r),yesno(cert));
end
fprintf(fid,'\\midrule\n');
fprintf(fid,'\\multicolumn{6}{l}{\\textit{(b) Theorem-certified location comparison, $p=6$, gain $k=30$}}\\\\\n');
labels={'Degree-first','Leaf-first'};
cases={'sparse_degree_p6_g30','sparse_leaf_p6_g30'};
for k=1:2
    r=strcmp(string(Ta.case_name),cases{k});
    cert=Ta.active_matrix_ok(r)&&Ta.rest_matrix_ok(r)&&Ta.sampling_ok(r);
    fprintf(fid,'%s & %.2f [%.2f, %.2f] & %.6f [%.6f, %.6f] & %.4f [%.4f, %.4f] & %.4f & %s \\\\\n', ...
        labels{k},Ta.updates_mean(r),Ta.updates_ci_low(r),Ta.updates_ci_high(r), ...
        Ta.J_mean(r),Ta.J_ci_low(r),Ta.J_ci_high(r), ...
        Ta.E_mean(r),Ta.E_ci_low(r),Ta.E_ci_high(r), ...
        Ta.active_matrix_maxeig(r),yesno(cert));
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\end{table*}\n');
fclose(fid);
fid=fopen(fullfile(out_dir,'table_R4_sampling.tex'),'w');
fprintf(fid,'\\begin{table*}[t]\n\\centering\n\\scriptsize\n');
fprintf(fid,'\\caption{Certified versus tested sampling periods.}\n');
fprintf(fid,'\\label{tab:r4_sampling}\n');
fprintf(fid,'\\renewcommand{\\arraystretch}{1.12}\n');
fprintf(fid,'\\begin{tabular}{cccccc}\n\\toprule\n');
fprintf(fid,'$h_s$ (s) & $h_s/h_{\\rm cert}$ & Sampling certificate & $N_u$ & $J_\\varrho$ & $J_u$ \\\\\n\\midrule\n');
for i=1:height(Th)
    cert=logical(Th.sampling_ok(i));
    fprintf(fid,'%.1e & %.2f & %s & %.2f [%.2f, %.2f] & %.6f [%.6f, %.6f] & %.4f [%.4f, %.4f] \\\\\n', ...
        Th.hs(i),Th.hs_over_hcert(i),yesno(cert), ...
        Th.updates_mean(i),Th.updates_ci_low(i),Th.updates_ci_high(i), ...
        Th.J_mean(i),Th.J_ci_low(i),Th.J_ci_high(i), ...
        Th.E_mean(i),Th.E_ci_low(i),Th.E_ci_high(i));
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\end{table*}\n');
fclose(fid);
fid=fopen(fullfile(out_dir,'table_R4_sensitivity.tex'),'w');
fprintf(fid,'\\begin{table*}[t]\n\\centering\n\\scriptsize\n');
fprintf(fid,'\\caption{One-at-a-time sensitivity of aperiodic AIPE-TC. The pin-count rows use the sparse degree-first graph; all other rows use the complete-graph reference configuration.}\n');
fprintf(fid,'\\label{tab:r4_sensitivity}\n');
fprintf(fid,'\\renewcommand{\\arraystretch}{1.08}\n');
fprintf(fid,'\\begin{tabular}{llcccc}\n\\toprule\n');
fprintf(fid,'Parameter & Value & $N_u$ & $J_\\varrho$ & $J_u$ & Matrix/sampling conditions \\\\\n\\midrule\n');
groups={'M15_sigma0','M15_sigma1','M15_lambda_c','M15_noise','M13_sparse_pin_count','M15_gain'};
display_names={'$\sigma_0$','$\sigma_1$','$\lambda_c$','Noise','$p$ (sparse)','$k$'};
for g=1:numel(groups)
    R=Tsen(strcmp(string(Tsen.group),groups{g}),:);
    [~,ord]=sort(R.parameter_value);
    R=R(ord,:);
    for i=1:height(R)
        cert=logical(R.active_matrix_ok(i)&&R.rest_matrix_ok(i)&&R.sampling_ok(i));
        fprintf(fid,'%s & %.6g & %.2f [%.2f, %.2f] & %.6f [%.6f, %.6f] & %.4f [%.4f, %.4f] & %s \\\\\n', ...
            display_names{g},R.parameter_value(i), ...
            R.updates_mean(i),R.updates_ci_low(i),R.updates_ci_high(i), ...
            R.J_mean(i),R.J_ci_low(i),R.J_ci_high(i), ...
            R.E_mean(i),R.E_ci_low(i),R.E_ci_high(i),yesno(cert));
    end
    if g<numel(groups), fprintf(fid,'\\addlinespace[2pt]\n'); end
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\end{table*}\n');
fclose(fid);
fprintf('Created candidate LaTeX tables in:\n%s\n',out_dir);
function s=yesno(tf)
if tf
    s='Yes';
else
    s='No';
end
end
