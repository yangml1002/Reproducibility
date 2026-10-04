clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
cd(script_dir);
required_files = { ...
    'simulation_engine_R4_remaining.m', ...
    'certificate_R4_remaining.m'};
for kk = 1:numel(required_files)
    f = fullfile(script_dir,required_files{kk});
    if exist(f,'file') ~= 2
        error('Required file not found: %s',f);
    end
end
outdir = fullfile(script_dir,'..','data');
if exist(outdir,'dir') ~= 7
    mkdir(outdir);
end
fprintf('Reviewer #4 remaining-experiment folder: %s\n',script_dir);
fprintf('Output folder: %s\n',outdir);
Nmc = 100;
hEM = 1e-5;
Tsim = 1.5;
seed = 1;
Aaper = 0.4;
base_hs = 5e-5;
base_sigma1 = 1e-4;
base_sigma0 = 5e-6;
base_duty = .75;
base_noise = .2;
base_pins = 4;
base_gain = 26;
base_qh = 0;
case_names = {};
case_group = {};
parameter_name = {};
parameter_value = [];
case_par = [];
case_names{end+1} = 'reference_complete';
case_group{end+1} = 'reference';
parameter_name{end+1} = 'baseline';
parameter_value(end+1) = NaN;
case_par(end+1,:) = [3 2 base_hs base_sigma1 base_sigma0 base_duty ...
    base_noise base_pins base_gain 1 base_qh Aaper];
case_names{end+1} = 'sparse_degree_p4';
case_group{end+1} = 'M13_sparse_pin_count';
parameter_name{end+1} = 'pins';
parameter_value(end+1) = 4;
case_par(end+1,:) = [3 2 base_hs base_sigma1 base_sigma0 base_duty ...
    base_noise 4 base_gain 2 base_qh Aaper];
case_names{end+1} = 'sparse_leaf_p4';
case_group{end+1} = 'M13_sparse_location';
parameter_name{end+1} = 'pinning_rule';
parameter_value(end+1) = 2;
case_par(end+1,:) = [3 2 base_hs base_sigma1 base_sigma0 base_duty ...
    base_noise 4 base_gain 3 base_qh Aaper];
p_sparse = [2 6];
for pp = p_sparse
    case_names{end+1} = sprintf('sparse_degree_p%d',pp);
    case_group{end+1} = 'M13_sparse_pin_count';
    parameter_name{end+1} = 'pins';
    parameter_value(end+1) = pp;
    case_par(end+1,:) = [3 2 base_hs base_sigma1 base_sigma0 base_duty ...
        base_noise pp base_gain 2 base_qh Aaper];
end
hs_list = [1e-5 5e-5 1e-4 5e-4 1e-3 5e-3 1e-2 2e-2 5e-2 1e-1];
for hs = hs_list
    case_names{end+1} = sprintf('hs_%g',hs);
    case_group{end+1} = 'M14_sampling';
    parameter_name{end+1} = 'hs';
    parameter_value(end+1) = hs;
    case_par(end+1,:) = [3 2 hs base_sigma1 base_sigma0 base_duty ...
        base_noise base_pins base_gain 1 base_qh Aaper];
end
sigma0_list = [1e-6 5e-6 2e-5];
for x = sigma0_list
    case_names{end+1} = sprintf('sigma0_%g',x);
    case_group{end+1} = 'M15_sigma0';
    parameter_name{end+1} = 'sigma0';
    parameter_value(end+1) = x;
    case_par(end+1,:) = [3 2 base_hs base_sigma1 x base_duty ...
        base_noise base_pins base_gain 1 base_qh Aaper];
end
sigma1_list = [2e-5 1e-4 3e-4];
for x = sigma1_list
    case_names{end+1} = sprintf('sigma1_%g',x);
    case_group{end+1} = 'M15_sigma1';
    parameter_name{end+1} = 'sigma1';
    parameter_value(end+1) = x;
    case_par(end+1,:) = [3 2 base_hs x base_sigma0 base_duty ...
        base_noise base_pins base_gain 1 base_qh Aaper];
end
duty_list = [.5 .75 .9];
for x = duty_list
    case_names{end+1} = sprintf('duty_%.2f',x);
    case_group{end+1} = 'M15_lambda_c';
    parameter_name{end+1} = 'lambda_c';
    parameter_value(end+1) = x;
    case_par(end+1,:) = [3 2 base_hs base_sigma1 base_sigma0 x ...
        base_noise base_pins base_gain 1 base_qh Aaper];
end
noise_list = [.1 .2 .4];
for x = noise_list
    case_names{end+1} = sprintf('noise_%.2f',x);
    case_group{end+1} = 'M15_noise';
    parameter_name{end+1} = 'noise';
    parameter_value(end+1) = x;
    case_par(end+1,:) = [3 2 base_hs base_sigma1 base_sigma0 base_duty ...
        x base_pins base_gain 1 base_qh Aaper];
end
gain_list = [18 26 34];
for x = gain_list
    case_names{end+1} = sprintf('gain_%d',x);
    case_group{end+1} = 'M15_gain';
    parameter_name{end+1} = 'gain';
    parameter_value(end+1) = x;
    case_par(end+1,:) = [3 2 base_hs base_sigma1 base_sigma0 base_duty ...
        base_noise base_pins x 1 base_qh Aaper];
end
fprintf('Total cases = %d, Nmc = %d\n',length(case_names),Nmc);
fprintf('Running memory-efficient paired simulations...\n');
run(fullfile(script_dir,'simulation_engine_R4_remaining.m'));
run(fullfile(script_dir,'certificate_R4_remaining.m'));
save(fullfile(outdir,'r4_remaining_results.mat'), ...
    'all_metrics','case_names','case_group','parameter_name','parameter_value', ...
    'case_par','schedule_stats','cert_values','schedule_check_pass', ...
    'Nmc','hEM','Tsim','seed','Aaper','case_elapsed');
if exist('tinv','file') == 2
    tcrit = tinv(.975,Nmc-1);
else
    tcrit = 1.98421695150868;
end
metric_names = { ...
    'updates','J_error','E_input','tail_error','final_error', ...
    'max_error','node_assignments','active_fraction','sensing_instants'};
fid = fopen(fullfile(outdir,'r4_cases.csv'),'w');
fprintf(fid,['case,group,parameter,parameter_value,strategy,schedule,hs,' ...
    'sigma1,sigma0,duty,noise,pins,gain,topology,q_h,aperiodic_amplitude\n']);
for i = 1:length(case_names)
    fprintf(fid,'%s,%s,%s,%.12g', ...
        case_names{i},case_group{i},parameter_name{i},parameter_value(i));
    fprintf(fid,',%.12g',case_par(i,:));
    fprintf(fid,'\n');
end
fclose(fid);
fid = fopen(fullfile(outdir,'r4_paths.csv'),'w');
fprintf(fid,['case,path,updates,J_error,E_input,tail_error,final_error,' ...
    'max_error,node_assignments,active_fraction,sensing_instants\n']);
for i = 1:length(case_names)
    for pth = 1:Nmc
        fprintf(fid,'%s,%d',case_names{i},pth);
        fprintf(fid,',%.12g',all_metrics(pth,:,i));
        fprintf(fid,'\n');
    end
end
fclose(fid);
fid = fopen(fullfile(outdir,'r4_summary.csv'),'w');
fprintf(fid,'case,group,parameter,parameter_value,metric,mean,sd,ci_low,ci_high,Nmc\n');
for i = 1:length(case_names)
    for im = 1:length(metric_names)
        x = all_metrics(:,im,i);
        mu = mean(x);
        sd = std(x,0);
        half = tcrit*sd/sqrt(Nmc);
        fprintf(fid,'%s,%s,%s,%.12g,%s,%.12g,%.12g,%.12g,%.12g,%d\n', ...
            case_names{i},case_group{i},parameter_name{i},parameter_value(i), ...
            metric_names{im},mu,sd,mu-half,mu+half,Nmc);
    end
end
fclose(fid);
fid = fopen(fullfile(outdir,'r4_certificates.csv'),'w');
fprintf(fid,['case,group,parameter,parameter_value,active_matrix_maxeig,' ...
    'rest_matrix_maxeig,muEhat,b,r,rmin,Lambda4,R_ultimate,' ...
    'h_certified_limit,active_matrix_ok,rest_matrix_ok,sampling_ok,' ...
    'finite_horizon_schedule_max_deficit,finite_horizon_schedule_check\n']);
for i = 1:length(case_names)
    fprintf(fid,'%s,%s,%s,%.12g', ...
        case_names{i},case_group{i},parameter_name{i},parameter_value(i));
    fprintf(fid,',%.12g',cert_values(i,1:9));
    fprintf(fid,',%d,%d,%d,%.12g,%d\n', ...
        round(cert_values(i,10)), ...
        round(cert_values(i,11)), ...
        round(cert_values(i,12)), ...
        schedule_stats(i,2), ...
        schedule_check_pass(i));
end
fclose(fid);
fid = fopen(fullfile(outdir,'r4_schedule_checks.csv'),'w');
fprintf(fid,['case,group,active_fraction,max_deficit_difference,ncycles,' ...
    'min_active_width,max_active_width,min_rest_width,max_rest_width\n']);
for i = 1:length(case_names)
    fprintf(fid,'%s,%s',case_names{i},case_group{i});
    fprintf(fid,',%.12g',schedule_stats(i,:));
    fprintf(fid,'\n');
end
fclose(fid);
sparse_names = {'reference_complete','sparse_degree_p4','sparse_leaf_p4', ...
    'sparse_degree_p2','sparse_degree_p6'};
fid = fopen(fullfile(outdir,'r4_sparse_summary.csv'),'w');
fprintf(fid,['case,updates_mean,updates_ci_low,updates_ci_high,' ...
    'J_mean,J_ci_low,J_ci_high,E_mean,E_ci_low,E_ci_high,' ...
    'tail_mean,active_matrix_maxeig,rest_matrix_maxeig,' ...
    'h_certified_limit,active_ok,rest_ok,sampling_ok\n']);
for ii = 1:length(sparse_names)
    idx = find(strcmp(case_names,sparse_names{ii}),1);
    vals = zeros(1,12);
    cols = [1 2 3 4];
    out = zeros(1,10);
    pos = 1;
    for jj = 1:length(cols)
        x = all_metrics(:,cols(jj),idx);
        mu = mean(x); sd = std(x,0); half=tcrit*sd/sqrt(Nmc);
        if jj <= 3
            out(pos:pos+2) = [mu mu-half mu+half];
            pos = pos+3;
        else
            tailmu = mu;
        end
    end
    fprintf(fid,'%s,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%d,%d,%d\n', ...
        sparse_names{ii}, ...
        out(1),out(2),out(3), ...
        out(4),out(5),out(6), ...
        out(7),out(8),out(9), ...
        tailmu, ...
        cert_values(idx,1),cert_values(idx,2),cert_values(idx,9), ...
        round(cert_values(idx,10)),round(cert_values(idx,11)),round(cert_values(idx,12)));
end
fclose(fid);
idxD = find(strcmp(case_names,'sparse_degree_p4'),1);
idxL = find(strcmp(case_names,'sparse_leaf_p4'),1);
fid = fopen(fullfile(outdir,'r4_sparse_paired_degree_minus_leaf.csv'),'w');
fprintf(fid,'metric,mean_difference,ci_low,ci_high\n');
for im = [1 2 3 4]
    delta = all_metrics(:,im,idxD)-all_metrics(:,im,idxL);
    mu = mean(delta);
    sd = std(delta,0);
    half = tcrit*sd/sqrt(Nmc);
    fprintf(fid,'%s,%.12g,%.12g,%.12g\n', ...
        metric_names{im},mu,mu-half,mu+half);
end
fclose(fid);
fid = fopen(fullfile(outdir,'r4_sensitivity_summary.csv'),'w');
fprintf(fid,['group,case,parameter,parameter_value,updates_mean,updates_ci_low,' ...
    'updates_ci_high,J_mean,J_ci_low,J_ci_high,E_mean,E_ci_low,E_ci_high,' ...
    'active_fraction,sensing_instants,active_matrix_ok,rest_matrix_ok,' ...
    'sampling_ok,h_certified_limit\n']);
for i = 1:length(case_names)
    is_sens = startsWith(case_group{i},'M15_') || ...
              strcmp(case_group{i},'M13_sparse_pin_count');
    if ~is_sens
        continue;
    end
    out = zeros(1,9);
    pos = 1;
    for im = [1 2 3]
        x = all_metrics(:,im,i);
        mu = mean(x); sd = std(x,0); half=tcrit*sd/sqrt(Nmc);
        out(pos:pos+2) = [mu mu-half mu+half];
        pos = pos+3;
    end
    fprintf(fid,'%s,%s,%s,%.12g', ...
        case_group{i},case_names{i},parameter_name{i},parameter_value(i));
    fprintf(fid,',%.12g',out);
    fprintf(fid,',%.12g,%.12g,%d,%d,%d,%.12g\n', ...
        mean(all_metrics(:,8,i)), ...
        mean(all_metrics(:,9,i)), ...
        round(cert_values(i,10)), ...
        round(cert_values(i,11)), ...
        round(cert_values(i,12)), ...
        cert_values(i,9));
end
fclose(fid);
disp('REVIEWER4_REMAINING_FINISHED');
disp('Use run_sampling_sweep_consistent.m for the final sampling-period sweep.');
