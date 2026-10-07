function recompute_certificates(work_dir)
script_dir=fileparts(mfilename('fullpath'));
data_dir=fullfile(script_dir,'..','data');
if nargin<1, work_dir=tempname; end
if ~exist(work_dir,'dir'), mkdir(work_dir); end
pending_files={};
comparison=table;
for prefix={'final_main','final_bias','sampling_consistent'}
    name=prefix{1};
    S=load(fullfile(data_dir,[name '_results.mat']),'case_names','case_par','hEM');
    compute_sampling(S,name,work_dir,script_dir);
    filename=[name '_certificates.csv'];
    T=readtable(fullfile(work_dir,filename),'TextType','string','VariableNamingRule','preserve');
    old=readtable(fullfile(data_dir,filename),'TextType','string','VariableNamingRule','preserve');
    verify_csv(old,T);
    validate_bound(T);
    pending_files{end+1}=filename;
    comparison=[comparison; comparison_rows(T,S.case_par)];
end
S=load(fullfile(data_dir,'r4_remaining_results.mat'), ...
    'case_names','case_par','schedule_stats','cert_values');
[T,V]=compute_r4(S,script_dir);
verify_matrix(S.cert_values,V);
old=readtable(fullfile(data_dir,'r4_certificates.csv'),'TextType','string','VariableNamingRule','preserve');
T.schedule_all_window_verified=legacy_all_window(S.case_par);
T.certified=certificate_status(T,S.case_par);
verify_csv(old,T);
updated=old;
updated.R_ultimate=T.R_ultimate;
extras={'d','lambda_c','Nc','Lambda2','rbar','Tc','R_minrate','R_av','schedule_all_window_verified', ...
    'improvement_ratio','improvement_percent','certified'};
for j=1:numel(extras), updated.(extras{j})=T.(extras{j}); end
validate_bound(updated);
writetable(updated,fullfile(work_dir,'r4_certificates.csv'));
pending_files{end+1}='r4_certificates.csv';
comparison=[comparison; comparison_rows(T,S.case_par)];
S=load(fullfile(data_dir,'sparse_certified_addon_results.mat'), ...
    'case_names','case_par','schedule_stats','cert_values');
[T,V]=compute_r4(S,script_dir);
T.schedule_all_window_verified=legacy_all_window(S.case_par);
T.certified=certificate_status(T,S.case_par);
verify_matrix(S.cert_values,V);
validate_bound(T);
writetable(T,fullfile(work_dir,'sparse_certified_addon_certificates.csv'));
pending_files{end+1}='sparse_certified_addon_certificates.csv';
comparison=[comparison; comparison_rows(T,S.case_par)];
W=load(fullfile(data_dir,'performance_match_workspace.mat'),'selected_lambda');
C=load(fullfile(data_dir,'performance_match_certificate.mat'),'cert_values','schedule_stats','average_control_rate_ok');
B=load(fullfile(data_dir,'final_main_results.mat'),'case_names','case_par');
idx=find(strcmp(B.case_names,'AIPE_aperiodic_A0.4'));
assert(isscalar(idx),'Baseline metadata is not unique.');
S.case_names={'Matched AIPE'};
S.case_par=B.case_par(idx,:);
S.case_par(6)=W.selected_lambda;
S.schedule_stats=C.schedule_stats;
[T,V]=compute_r4(S,script_dir);
verify_matrix(C.cert_values,V);
T.schedule_all_window_verified=repmat(logical(C.average_control_rate_ok),height(T),1);
T.certified=certificate_status(T,S.case_par);
validate_bound(T);
writetable(T,fullfile(work_dir,'performance_match_certificates.csv'));
pending_files{end+1}='performance_match_certificates.csv';
comparison=[comparison; comparison_rows(T,S.case_par)];
S=load(fullfile(data_dir,'parameter_sensitivity_workspace.mat'), ...
    'case_par','schedule_stats','cert_values');
F=readtable(fullfile(data_dir,'parameter_sensitivity_full.csv'),'TextType','string');
S.case_names=cellstr(F.parameter+"_"+string(F.value));
[T,V]=compute_r4(S,script_dir);
verify_matrix(S.cert_values,V);
T.schedule_all_window_verified=logical(F.average_control_rate_ok);
T.certified=certificate_status(T,S.case_par);
validate_bound(T);
writetable(T,fullfile(work_dir,'parameter_sensitivity_certificates.csv'));
pending_files{end+1}='parameter_sensitivity_certificates.csv';
comparison=[comparison; comparison_rows(T,S.case_par)];
writetable(comparison,fullfile(work_dir,'certificate_bound_comparison.csv'));
pending_files{end+1}='certificate_bound_comparison.csv';
for j=1:numel(pending_files)
    copyfile(fullfile(work_dir,pending_files{j}),fullfile(data_dir,pending_files{j}),'f');
end
fid=fopen(fullfile(work_dir,'deterministic_verification_report.txt'),'w');
fprintf(fid,'All certificate invariance and average-rate consistency checks passed.\n');
fprintf(fid,'Matrix tests, muEhat, b, r, Lambda4, sampling limits and existing statuses unchanged.\n');
fprintf(fid,'All candidates validated before publication.\n');
fprintf(fid,'Monte Carlo simulations rerun: NO\nBrownian paths regenerated: NO\nsimulation_engine executed: NO\n');
fprintf(fid,'Embedded certificate arrays in precomputed stochastic MAT files remain untouched; current certificates are regenerated CSV files.\n');
fclose(fid);
disp(comparison(contains(comparison.case_name,'bias_q'),:));
fprintf('DETERMINISTIC_CERTIFICATES_VERIFIED: %d cases; %d CSV files.\n',height(comparison),numel(pending_files));
end

function compute_sampling(S,result_prefix,certificate_output_dir,script_dir)
case_names=S.case_names;
case_par=S.case_par;
hEM=S.hEM;
run(fullfile(script_dir,'sampling_certificate.m'));
end

function [T,cert_values]=compute_r4(S,script_dir)
case_names=S.case_names;
case_par=S.case_par;
schedule_stats=S.schedule_stats;
run(fullfile(script_dir,'certificate_R4_remaining.m'));
T=array2table(cert_values,'VariableNames',{'active_matrix_maxeig', ...
    'rest_matrix_maxeig','muEhat','b','r','rmin','Lambda4','R_ultimate', ...
    'h_certified_limit','active_matrix_ok','rest_matrix_ok','sampling_ok'});
T=addvars(T,string(case_names(:)),'Before',1,'NewVariableNames','case_name');
T.finite_horizon_schedule_max_deficit=schedule_stats(:,2);
T.finite_horizon_schedule_check=schedule_check_pass;
T.d=cert_values(:,6);
T.lambda_c=case_par(:,6);
T.Nc=repmat(Nc,height(T),1);
T.Lambda2=repmat(Lambda2,height(T),1);
T.rbar=rbar_values;
T.Tc=Tc_values;
T.R_minrate=R_minrate_values;
T.R_av=R_av_values;
T.improvement_ratio=improvement_ratio_values;
T.improvement_percent=improvement_percent_values;
end

function verify_csv(A,B)
assert(isequal(string(A{:,1}),string(B{:,1})),'Case identity/order changed.');
fields={'active_matrix_maxeig','rest_matrix_maxeig','muEhat','b','r','rmin', ...
    'Lambda4','h_certified_limit','certified','schedule_certified', ...
    'schedule_deficit_bound','aperiodic_amplitude','active_matrix_ok', ...
    'rest_matrix_ok','sampling_ok','finite_horizon_schedule_max_deficit', ...
    'finite_horizon_schedule_check'};
for j=1:numel(fields)
    f=fields{j};
    if ismember(f,A.Properties.VariableNames) && ismember(f,B.Properties.VariableNames)
        verify_numbers(A.(f),B.(f),f);
    end
end
end

function verify_matrix(A,B)
assert(isequal(size(A),size(B)),'Certificate matrix shape changed.');
verify_numbers(A(:,[1:7 9:12]),B(:,[1:7 9:12]),'fixed certificate matrix columns');
end

function verify_numbers(a,b,name)
a=double(a); b=double(b);
assert(isequal(isnan(a),isnan(b)) && isequal(isinf(a),isinf(b)), ...
    'Nonfinite status changed: %s',name);
finite=isfinite(a) & isfinite(b);
tol=5e-11*max(max(abs(a(finite)),abs(b(finite))),1e-12);
assert(all(abs(a(finite)-b(finite))<=tol,'all'),'Protected quantity changed: %s',name);
end

function validate_bound(T)
ok=isfinite(T.r) & isfinite(T.Lambda4);
r=T.r(ok); d=T.d(ok); rb=T.rbar(ok); L=T.Lambda4(ok);
tc=T.Tc(ok); av=T.R_av(ok); old=T.R_minrate(ok);
tol=1e-10*max(abs(old),1e-12);
assert(all(d>0 & rb>=d & tc>0));
assert(all(av<=old+tol));
same=r<=T.Lambda2(ok);
assert(all(abs(av(same)-old(same))<=tol(same)));
strict=r>T.Lambda2(ok) & L>0;
assert(all(av(strict)<old(strict)));
I=(1-exp(-d.*tc))./d+exp(-d.*tc)./rb;
positive=L>0;
assert(all(abs(av(positive)./(2*L(positive))-I(positive))<=1e-10*max(I(positive),1e-12)));
assert(all(abs(T.R_ultimate(ok)-av)<=tol));
end

function C=comparison_rows(T,case_par)
n=height(T);
case_name=string(T{:,1});
q_h=case_par(:,11);
lambda_c=case_par(:,6);
Nc=repmat(.05,n,1); Lambda2=repmat(.2,n,1);
r=T.r; d=T.d; rbar=T.rbar; Tc=T.Tc; Lambda4=T.Lambda4;
R_minrate=T.R_minrate; R_av=T.R_av;
ratio_Rav_to_Rmin=T.improvement_ratio;
reduction_percent=T.improvement_percent;
if ismember('certified',T.Properties.VariableNames)
    certified=logical(T.certified);
elseif ismember('schedule_all_window_verified',T.Properties.VariableNames)
    certified=logical(T.active_matrix_ok) & logical(T.rest_matrix_ok) & ...
        logical(T.sampling_ok) & T.schedule_all_window_verified & case_par(:,1)==3;
else
    error('Missing independent all-window schedule evidence.');
end
C=table(case_name,q_h,lambda_c,Nc,r,Lambda2,d,rbar,Tc,Lambda4, ...
    R_minrate,R_av,ratio_Rav_to_Rmin,reduction_percent,certified);
C.R_ultimate=R_av;
C.improvement_ratio=ratio_Rav_to_Rmin;
C.improvement_percent=reduction_percent;
end

function certified=certificate_status(T,case_par)
certified=logical(T.active_matrix_ok) & logical(T.rest_matrix_ok) & ...
    logical(T.sampling_ok) & T.schedule_all_window_verified & case_par(:,1)==3;
end

function verified=legacy_all_window(case_par)
lambda=case_par(:,6);
if size(case_par,2)>=12, amp=case_par(:,12); else, amp=repmat(.4,size(lambda)); end
Pmax=4e-5*round(.05*(1+amp)/(4e-5));
Pmax(case_par(:,2)==1)=.05;
exact_quarter_fraction=abs(4*lambda-round(4*lambda))<1e-12;
verified=exact_quarter_fraction & lambda.*(1-lambda).*Pmax<=.05+1e-14;
end
