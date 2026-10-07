clear; clc;
script_dir=fileparts(mfilename('fullpath'));
root_dir=fileparts(script_dir); data_dir=fullfile(root_dir,'data'); table_dir=fullfile(root_dir,'tables');
S=readtable(fullfile(data_dir,'parameter_sensitivity_full.csv'),'TextType','string');
names={'$\sigma_0$','$\sigma_1$','$\lambda_c$','Noise intensity','Pinning gain'};
fid=start_table(table_dir,'table_parameter_sensitivity.tex','One-at-a-time parameter-sensitivity results.','tab:parameter_sensitivity','cccccc','Parameter & Value & $N_u$ & $J_\varrho$ & $J_u$ & Theorem 1 conditions');
for i=1:height(S)
 param=char(S.parameter(i));
 idx=find(strcmp({'sigma0','sigma1','lambda_c','noise','gain'},param),1);
 assert(~isempty(idx),'Unknown sensitivity parameter: %s',param);
 if strcmp(param,'sigma0') || strcmp(param,'sigma1')
  value_string=sci_tex(S.value(i));
 else
  value_string=sprintf('$%g$',S.value(i));
 end
 row={names{idx},value_string,stat(S,i,'Nu',2),stat(S,i,'Jrho',6),stat(S,i,'Ju',4),char(S.theorem_conditions(i))};
 put(fid,[strjoin(row,' & ') ' \\']);
end
finish(fid);
S=readtable(fullfile(data_dir,'performance_match_summary.csv'),'TextType','string');
D=readtable(fullfile(data_dir,'performance_match_paired.csv'),'TextType','string');
B=readtable(fullfile(data_dir,'same_threshold_paired.csv'),'TextType','string');
fid=start_table(table_dir,'table_performance_matched.tex','Near-error-matched and paired Monte Carlo comparison.','tab:performance_matched','cccccccc','Method & $\lambda_c$ & $\xi$ & $N_u$ & $J_\varrho$ & $J_u$ & $R_c$ & $N_s$');
for i=1:height(S)
 if i==1, rate='--'; else, rate=sprintf('%.5f',S.lambda_c(i)); end
 row={char(S.method(i)),rate,sprintf('%g',S.xi(i)),stat(S,i,'Nu',2),stat(S,i,'Jrho',6),stat(S,i,'Ju',4),sprintf('%g',S.Rc(i)),sprintf('%g',S.Ns(i))};
 put(fid,[strjoin(row,' & ') ' \\']);
end
put(fid,'\midrule');put(fid,'Metric & \multicolumn{3}{c}{AIPE-TC $-$ CP-PETC} & \multicolumn{4}{c}{95\% paired CI} \\');
paired_rows(fid,D,[2 8 6]);
put(fid,'\midrule');put(fid,'\multicolumn{8}{c}{Same-threshold baseline ($\lambda_c=0.75$ for AIPE-TC), AIPE-TC $-$ CP-PETC}\\');
paired_rows(fid,B,[2 6 6]);
finish(fid,'\par\smallskip The AIPE-TC activation rate is selected from the tested operating points to minimize the symmetric relative gap in the mean synchronization error with respect to the CP-PETC reference.');
S=readtable(fullfile(data_dir,'timestep_refinement_matched_summary.csv'),'TextType','string');
fid=start_table(table_dir,'table_timestep_refinement.tex','Time-step refinement check.','tab:timestep_refinement','ccccc','Method & $\Delta t$ & $N_u$ & $J_\varrho$ & $J_u$');
for i=1:height(S)
 row={char(S.method(i)),sci_tex(S.step(i)),stat(S,i,'Nu',2),stat(S,i,'Jrho',6),stat(S,i,'Ju',4)};
 put(fid,[strjoin(row,' & ') ' \\']);
end
finish(fid);
function fid=start_table(folder,name,caption,label,columns,heading)
 fid=fopen(fullfile(folder,name),'w');assert(fid>=0);
 put(fid,'\begin{table}[H]');put(fid,'\centering');put(fid,'\setlength{\tabcolsep}{3pt}');
 put(fid,['\caption{' caption '}\label{' label '}']);put(fid,'\scriptsize');
 put(fid,['\begin{tabular}{' columns '}']);put(fid,'\toprule');put(fid,[heading ' \\']);put(fid,'\midrule');
end
function finish(fid,note)
 put(fid,'\bottomrule');put(fid,'\end{tabular}');
 if nargin>1, put(fid,note); end
 put(fid,'\end{table}');fclose(fid);
end
function put(fid,s)
 fprintf(fid,'%s\n',s);
end
function s=stat(T,i,p,n)
 f=['%.' num2str(n) 'f [%.' num2str(n) 'f, %.' num2str(n) 'f]'];
 s=sprintf(f,T.([p '_mean'])(i),T.([p '_ci_low'])(i),T.([p '_ci_high'])(i));
end
function paired_rows(fid,T,precision)
 names={'$\Delta N_u$','$\Delta J_\varrho$','$\Delta J_u$'};
 for i=1:3
  f=['%.' num2str(precision(i)) 'f'];
  value=sprintf(f,T.mean_difference(i));
  low=sprintf(f,T.ci_low(i));high=sprintf(f,T.ci_high(i));
  put(fid,[names{i} ' & \multicolumn{3}{c}{' value '} & \multicolumn{4}{c}{[' low ', ' high ']} \\']);
 end
end
function s=sci_tex(x)
 if x==0
  s='$0$';
  return;
 end
 exponent=floor(log10(abs(x)));
 coeff=x/(10^exponent);
 if abs(coeff-1)<1e-12
  s=sprintf('$10^{%d}$',exponent);
 else
  s=sprintf('$%g\\times10^{%d}$',coeff,exponent);
 end
end
