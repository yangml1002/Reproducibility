clear; clc; close all;
script_dir=fileparts(mfilename('fullpath')); stage_dir=fileparts(script_dir); data_dir=fullfile(stage_dir,'data');
Nmc=100; hEM=1e-5; Tsim=1.5; seed=1; Aaper=.4;
base=[3 2 5e-5 1e-4 5e-6 .75 .2 4 26 1 0 .4];
case_names={'AIPE baseline','CP reference'}; case_par=[base;1 1 base(3:end)];
run(fullfile(script_dir,'simulation_engine.m'));
baseline_metrics=all_metrics; baseline_mean=mean(all_metrics(:,1:3,1),1);
expected=[372.28 .0465761381918816 7.1970799158426];
save(fullfile(data_dir,'baseline_reproduction.mat'),'baseline_metrics','baseline_mean','expected');
if any(abs(baseline_mean-expected)>[1e-10 1e-12 1e-10]), error('STOP: formal baseline reproduction failed.'); end
cp_metrics=all_metrics(:,:,2); J_target=mean(cp_metrics(:,2));
fprintf('BASELINE_REPRODUCTION_PASS\n');
groups={'sigma0','sigma1','lambda_c','noise','gain'}; grids={[1e-6 5e-6 2e-5],[2e-5 1e-4 3e-4],[.5 .75 .9],[.1 .2 .4],[18 26 34]}; cols=[5 4 6 7 9];
case_par=repmat(base,15,1); case_names=cell(15,1); parameter=strings(15,1); value=zeros(15,1);
for g=1:5
 for j=1:3
  ii=(g-1)*3+j; parameter(ii)=groups{g}; value(ii)=grids{g}(j); case_par(ii,cols(g))=value(ii); case_names{ii}=sprintf('%s_%g',groups{g},value(ii));
 end
end
full_case_par=case_par;
run(fullfile(script_dir,'simulation_engine.m'));
sensitivity_metrics=all_metrics;
case_par=full_case_par;
schedule_records=cell(15,1);
for ii=1:15, schedule_records{ii}=make_boundaries(case_par(ii,6),hEM,Tsim); end
for ii=[2 5 8 11 14]
 if max(abs(all_metrics(:,1:3,ii)-baseline_metrics(:,1:3,1)),[],'all')>1e-12, error('STOP: sensitivity baseline identity failed.'); end
end
schedule_stats=zeros(15,2);
for ii=1:15, schedule_stats(ii,2)=endpoint_deficit(schedule_records{ii},case_par(ii,6)); end
run(fullfile(script_dir,'certificate_R4_remaining.m'));
all_window_bound=case_par(:,6).*(1-case_par(:,6))*.07;
actual_max_rest_rise=zeros(15,1); cycle_fraction_ok=false(15,1);
for ii=1:15
 B=schedule_records{ii}; a=B(:,2)-B(:,1); periods=B(:,3)-B(:,1);
 cycle_fraction_ok(ii)=all(a+1e-12>=case_par(ii,6)*periods);
 actual_max_rest_rise(ii)=max(case_par(ii,6)*(periods-a));
end
average_control_rate_ok=cycle_fraction_ok & all_window_bound<=.05+1e-12;
S=stats_table(sensitivity_metrics); active_matrix_ok=logical(cert_values(:,10)); rest_matrix_ok=logical(cert_values(:,11)); sampling_ok=logical(cert_values(:,12)); theorem_conditions=repmat("No",15,1); theorem_conditions(average_control_rate_ok & active_matrix_ok & rest_matrix_ok & sampling_ok)="Yes";
S=[table(parameter,value) S table(active_matrix_ok,rest_matrix_ok,sampling_ok,average_control_rate_ok,all_window_bound,actual_max_rest_rise,theorem_conditions)];
writetable(S,fullfile(data_dir,'parameter_sensitivity_full.csv'));
save(fullfile(data_dir,'parameter_sensitivity_workspace.mat'),'sensitivity_metrics','case_par','cert_values','schedule_stats');

function T=stats_table(A)
 n=size(A,3); V=zeros(n,9); for i=1:n, for j=1:3, x=A(:,j,i); ci=paired_ci(x); V(i,(j-1)*3+(1:3))=ci; end, end
 T=array2table(V,'VariableNames',{'Nu_mean','Nu_ci_low','Nu_ci_high','Jrho_mean','Jrho_ci_low','Jrho_ci_high','Ju_mean','Ju_ci_low','Ju_ci_high'});
end
function deficit=endpoint_deficit(B,duty)
 t=unique([0;B(:)]); active=zeros(size(t)); for j=1:size(B,1), active=active+max(0,min(t,B(j,2))-B(j,1)); end
 z=duty*t-active; deficit=max(z-cummin(z));
end


function B=make_boundaries(duty,h,T)
 n=round(T/h); left=0; mc=0; B=[];
 while left<n
 width=max(4,4*round(.05*(1+.4*sin(sqrt(2)*mc))/(4*h)));width=min(width,n-left);right=left+width;
 B(end+1,:)=[left min(left+ceil(duty*width),n) right]*h;left=right;mc=mc+1;
 end
end
