clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
cd(script_dir);
required = {'simulation_engine_R4_remaining.m','certificate_R4_remaining.m'};
for k=1:numel(required)
    if exist(fullfile(script_dir,required{k}),'file')~=2
        error('Missing required file: %s',required{k});
    end
end
outdir = fullfile(script_dir,'..','data');
if exist(outdir,'dir')~=7, mkdir(outdir); end
Nmc = 100;
hEM = 1e-5;
Tsim = 1.5;
seed = 1;
Aaper = .4;
case_names = {'sparse_degree_p6_g30','sparse_leaf_p6_g30'};
case_par = [ ...
    3 2 5e-5 1e-4 5e-6 .75 .2 6 30 2 0 Aaper;
    3 2 5e-5 1e-4 5e-6 .75 .2 6 30 3 0 Aaper];
fprintf('--- Certified sparse add-on: %d cases, Nmc=%d ---\n',length(case_names),Nmc);
run(fullfile(script_dir,'simulation_engine_R4_remaining.m'));
run(fullfile(script_dir,'certificate_R4_remaining.m'));
if exist('tinv','file')==2
    tcrit=tinv(.975,Nmc-1);
else
    tcrit=1.98421695150868;
end
metric_names={'updates','J_error','E_input','tail_error'};
fid=fopen(fullfile(outdir,'sparse_certified_addon_summary.csv'),'w');
fprintf(fid,['case,updates_mean,updates_ci_low,updates_ci_high,' ...
    'J_mean,J_ci_low,J_ci_high,E_mean,E_ci_low,E_ci_high,' ...
    'tail_mean,tail_ci_low,tail_ci_high,active_matrix_maxeig,' ...
    'rest_matrix_maxeig,muEhat,h_certified_limit,' ...
    'active_matrix_ok,rest_matrix_ok,sampling_ok\n']);
for i=1:length(case_names)
    vals=zeros(1,12); pos=1;
    for im=1:4
        x=all_metrics(:,im,i);
        mu=mean(x); sd=std(x,0); half=tcrit*sd/sqrt(Nmc);
        vals(pos:pos+2)=[mu mu-half mu+half];
        pos=pos+3;
    end
    fprintf(fid,'%s',case_names{i});
    fprintf(fid,',%.12g',vals);
    fprintf(fid,',%.12g,%.12g,%.12g,%.12g,%d,%d,%d\n', ...
        cert_values(i,1),cert_values(i,2),cert_values(i,3),cert_values(i,9), ...
        round(cert_values(i,10)),round(cert_values(i,11)),round(cert_values(i,12)));
end
fclose(fid);
fid=fopen(fullfile(outdir,'sparse_certified_addon_paired.csv'),'w');
fprintf(fid,'metric,mean_difference,ci_low,ci_high\n');
for im=1:4
    delta=all_metrics(:,im,1)-all_metrics(:,im,2);
    mu=mean(delta); sd=std(delta,0); half=tcrit*sd/sqrt(Nmc);
    fprintf(fid,'%s,%.12g,%.12g,%.12g\n',metric_names{im},mu,mu-half,mu+half);
end
fclose(fid);
save(fullfile(outdir,'sparse_certified_addon_results.mat'), ...
    'all_metrics','case_names','case_par','schedule_stats','cert_values', ...
    'schedule_check_pass','Nmc','hEM','Tsim','seed','Aaper');
fprintf('\nExpected algebraic target for both cases:\n');
for i=1:length(case_names)
    fprintf('%s: active eig = %.9g, rest eig = %.9g, muEhat = %.9g, hcert = %.9g, flags = [%d %d %d]\n', ...
        case_names{i},cert_values(i,1),cert_values(i,2),cert_values(i,3), ...
        cert_values(i,9),round(cert_values(i,10)),round(cert_values(i,11)),round(cert_values(i,12)));
end
disp('SPARSE_CERTIFIED_ADDON_FINISHED');
