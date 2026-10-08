clear; clc;
script_dir=fileparts(mfilename('fullpath'));
addpath(script_dir);
dd=fullfile(fileparts(script_dir),'data');
tol=1e-12;
B=load(fullfile(dd,'baseline_reproduction.mat'));
S=load(fullfile(dd,'parameter_sensitivity_workspace.mat'));
baseline_pathwise_max_abs=0;
for ii=[2 5 8 11 14]
    baseline_pathwise_max_abs=max(baseline_pathwise_max_abs, ...
        max(abs(S.sensitivity_metrics(:,1:3,ii)-B.baseline_metrics(:,1:3,1)),[],'all'));
end
assert(baseline_pathwise_max_abs<tol,'Sensitivity baseline identity failed.');
G=readtable(fullfile(dd,'performance_match_activation_grid.csv'));
rates=unique([0.5;0.75;0.9;G.lambda_c]);
all_window_bound=rates.*(1-rates)*0.07;
actual_max_rest_rise=zeros(size(rates));
endpoint_max_deficit=zeros(size(rates));
cycle_fraction_ok=false(size(rates));
h=1e-5; n=150000;
for i=1:numel(rates)
    left=0; mc=0; bd=[];
    while left<n
        width=max(4,4*round(0.05*(1+0.4*sin(sqrt(2)*mc))/(4*h)));
        width=min(width,n-left);
        right=left+width;
        bd(end+1,:)=[left min(left+ceil(rates(i)*width),n) right]*h;
        left=right; mc=mc+1;
    end
    P=bd(:,3)-bd(:,1); a=bd(:,2)-bd(:,1);
    cycle_fraction_ok(i)=all(a+tol>=rates(i)*P);
    actual_max_rest_rise(i)=max(rates(i)*(P-a));
    t=unique([0;bd(:)]); A=zeros(size(t));
    for j=1:size(bd,1)
        A=A+max(0,min(t,bd(j,2))-bd(j,1));
    end
    d=rates(i)*t-A;
    endpoint_max_deficit(i)=max(d-cummin(d));
end
Nc=0.05*ones(size(rates));
pass=cycle_fraction_ok & all_window_bound<=Nc & endpoint_max_deficit<=all_window_bound+tol;
assert(all(pass),'All-window schedule verification failed.');
schedule=table(rates,all_window_bound,actual_max_rest_rise,endpoint_max_deficit,Nc,cycle_fraction_ok,pass, ...
    'VariableNames',{'lambda_c','all_window_bound','actual_max_rest_rise','endpoint_max_deficit','Nc','cycle_fraction_ok','pass'});
old_schedule=readtable(fullfile(dd,'schedule_all_window_checks.csv'));
assert(isequal(old_schedule.Properties.VariableNames,schedule.Properties.VariableNames));
assert(isequal(size(old_schedule),size(schedule)));
assert(max(abs(table2array(old_schedule)-table2array(schedule)),[],'all')<tol, ...
    'Existing schedule numerical values differ; no output written.');
C=load(fullfile(dd,'performance_match_certificate.mat'));
W=load(fullfile(dd,'performance_match_workspace.mat'));
F=load(fullfile(dd,'timestep_refinement_matched_workspace.mat'));
assert(C.average_control_rate_ok && C.cycle_fraction_ok && ...
    all(C.cert_values(10:12)) && C.schedule_check_pass,'Selected certificate failed.');
assert(abs(W.selected_lambda-0.99375)<tol && abs(F.selected_lambda-0.99375)<tol);
table9_table10_pathwise_max_abs=max([ ...
    max(abs(F.coarse_metrics(:,1:3,1)-W.cp_metrics(:,1:3)),[],'all'), ...
    max(abs(F.coarse_metrics(:,1:3,2)-W.ai_metrics(:,1:3)),[],'all')]);
assert(table9_table10_pathwise_max_abs<tol);
assert(F.bridge_max_abs_error<tol && F.coarse_identity_max_abs<tol);
Q=readtable(fullfile(dd,'performance_match_paired.csv'));
paired_ci_max_abs_error=0;
for j=1:3
    ci=paired_ci(W.ai_metrics(:,j)-W.cp_metrics(:,j));
    paired_ci_max_abs_error=max(paired_ci_max_abs_error, ...
        max(abs(ci-[Q.mean_difference(j),Q.ci_low(j),Q.ci_high(j)])));
end
assert(paired_ci_max_abs_error<tol);
M=readtable(fullfile(dd,'performance_match_summary.csv'),'TextType','string');
ic=find(M.method=="CP-PETC"); ia=find(M.method=="Aperiodic AIPE-TC");
assert(isscalar(ic) && isscalar(ia));
near_error_symmetric_gap=abs(M.Jrho_mean(ia)-M.Jrho_mean(ic))/ ...
    ((M.Jrho_mean(ia)+M.Jrho_mean(ic))/2);
assert(near_error_symmetric_gap<0.005 && abs(near_error_symmetric_gap-W.relative_gap)<tol);
active_time_reduction=1-M.Rc(ia)/M.Rc(ic);
sensing_reduction=1-M.Ns(ia)/M.Ns(ic);
assert(abs(active_time_reduction-0.006173333333333)<tol);
assert(abs(sensing_reduction-0.005766666666667)<tol);
selected=find(abs(rates-W.selected_lambda)<tol,1);
assert(~isempty(selected));
check=["baseline_pathwise_max_abs";"selected_all_window_bound"; ...
    "selected_endpoint_max_deficit";"table9_table10_pathwise_max_abs"; ...
    "bridge_pair_sum_max_abs";"coarse_identity_max_abs";"paired_ci_max_abs_error"; ...
    "near_error_symmetric_gap";"active_time_reduction";"sensing_reduction"];
value=[baseline_pathwise_max_abs;all_window_bound(selected); ...
    endpoint_max_deficit(selected);table9_table10_pathwise_max_abs; ...
    F.bridge_max_abs_error;F.coarse_identity_max_abs;paired_ci_max_abs_error; ...
    near_error_symmetric_gap;active_time_reduction;sensing_reduction];
tolerance=[tol;0.05;all_window_bound(selected)+tol;tol;tol;tol;tol;0.005;NaN;NaN];
pass=[value(1)<tol;value(2)<=tolerance(2);value(3)<=tolerance(3); ...
    value(4)<tol;value(5)<tol;value(6)<tol;value(7)<tol;value(8)<0.005;true;true];
assert(all(pass));
writetable(schedule,fullfile(dd,'schedule_all_window_checks.csv'));
writetable(table(check,value,tolerance,pass),fullfile(dd,'reproducibility_checks.csv'));
fprintf('REPRODUCIBILITY_VERIFICATION_PASS\n');
