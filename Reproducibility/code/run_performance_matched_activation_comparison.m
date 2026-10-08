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
lambda_all=[]; J_all=[]; metrics_all=cell(0); nominal=[.90 .95 .975 .9875 .99375];
for jj=1:numel(nominal)
 lam=nominal(jj); case_par=base; case_par(6)=lam; case_names={sprintf('AIPE_lambda_%g',lam)};
 run(fullfile(script_dir,'simulation_engine.m'));
 lambda_all(end+1,1)=lam; J_all(end+1,1)=mean(all_metrics(:,2,1)); metrics_all{end+1}=all_metrics(:,:,1);
 % All five prescribed initial rates are tested before refinement.
end
for refine=1:6
 above=find(J_all>J_target); below=find(J_all<=J_target);
 if isempty(above)||isempty(below), break; end
 left=max(lambda_all(above)); right=min(lambda_all(below)); lam=(left+right)/2;
 case_par=base; case_par(6)=lam; case_names={sprintf('AIPE_refine_%g',lam)};
 run(fullfile(script_dir,'simulation_engine.m'));
 lambda_all(end+1,1)=lam; J_all(end+1,1)=mean(all_metrics(:,2,1)); metrics_all{end+1}=all_metrics(:,:,1);
end
gaps=abs(J_all-J_target)./((J_all+J_target)/2); [relative_gap,best]=min(gaps); selected_lambda=lambda_all(best); ai_metrics=metrics_all{best};
grid_table=table(lambda_all,J_all,gaps,'VariableNames',{'lambda_c','Jrho_mean','relative_gap'}); writetable(grid_table,fullfile(data_dir,'performance_match_activation_grid.csv'));
save(fullfile(data_dir,'performance_match_workspace.mat'),'cp_metrics','ai_metrics','selected_lambda','relative_gap','J_target','grid_table');
if relative_gap>.005, error('STOP: activation matching remains above 0.5 percent.'); end
pair_metrics=cat(3,cp_metrics,ai_metrics); S=stats_table(pair_metrics); method=["CP-PETC";"Aperiodic AIPE-TC"]; lambda_c=[1;selected_lambda]; xi=[1;1]; Rc=[mean(cp_metrics(:,8));mean(ai_metrics(:,8))]; Ns=[mean(cp_metrics(:,9));mean(ai_metrics(:,9))];
writetable([table(method,lambda_c,xi) S table(Rc,Ns)],fullfile(data_dir,'performance_match_summary.csv'));
pd=zeros(3,3); for j=1:3, pd(j,:)=paired_ci(ai_metrics(:,j)-cp_metrics(:,j)); end
metric=["Delta_Nu";"Delta_Jrho";"Delta_Ju"]; writetable(table(metric,pd(:,1),pd(:,2),pd(:,3),'VariableNames',{'metric','mean_difference','ci_low','ci_high'}),fullfile(data_dir,'performance_match_paired.csv'));
case_par=base; case_par(6)=selected_lambda; case_names={'Matched AIPE'}; boundaries=[]; left=0; mc=0; n=round(Tsim/hEM);
while left<n
 width=max(4,4*round(.05*(1+.4*sin(sqrt(2)*mc))/(4*hEM))); width=min(width,n-left); right=left+width; boundaries(end+1,:)=[left left+ceil(selected_lambda*width) right]*hEM; left=right; mc=mc+1;
end
schedule_stats=[0 endpoint_deficit(boundaries,selected_lambda)]; run(fullfile(script_dir,'certificate_R4_remaining.m'));
periods=boundaries(:,3)-boundaries(:,1); a=boundaries(:,2)-boundaries(:,1);
cycle_fraction_ok=all(a+1e-12>=selected_lambda*periods);
all_window_bound=selected_lambda*(1-selected_lambda)*.07;
actual_max_rest_rise=max(selected_lambda*(periods-a));
average_control_rate_ok=cycle_fraction_ok && all_window_bound<=.05;
assert(average_control_rate_ok,'STOP: all-window rate verification failed.');
save(fullfile(data_dir,'performance_match_certificate.mat'),'cert_values','schedule_stats','schedule_check_pass','boundaries','all_window_bound','actual_max_rest_rise','cycle_fraction_ok','average_control_rate_ok');
fprintf('SELECTED_LAMBDA=%.12g GAP=%.9g CERTIFICATE=%d\n',selected_lambda,relative_gap,all(cert_values(10:12)) && schedule_check_pass);

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
