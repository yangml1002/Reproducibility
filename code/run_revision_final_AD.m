clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
cd(script_dir);
required_files = {'simulation_engine.m','write_results.m','sampling_certificate.m'};
for kk = 1:numel(required_files)
    f = fullfile(script_dir,required_files{kk});
    if exist(f,'file') ~= 2
        error('Required file not found: %s',f);
    end
end
outdir = fullfile(script_dir,'..','data');
if ~exist(outdir,'dir'), mkdir(outdir); end
fprintf('Final A-D folder: %s\n',script_dir);
fprintf('Using simulation engine: %s\n',which('simulation_engine.m'));
Nmc = 100;
hEM = 1e-5;
Tsim = 1.5;
seed = 1;
Aaper = 0.4;
case_names = {};
case_par = [];
case_names{end+1} = 'CP_original';
case_par(end+1,:) = [1 1 5e-5 3e-5 2e-6 .75 .2 4 26 1 0 0];
case_names{end+1} = 'AIP_periodic';
case_par(end+1,:) = [2 1 5e-5 1e-4 5e-6 .75 .2 4 26 1 0 0];
case_names{end+1} = 'AIPE_periodic';
case_par(end+1,:) = [3 1 5e-5 1e-4 5e-6 .75 .2 4 26 1 0 0];
case_names{end+1} = 'AIP_aperiodic_A0.4';
case_par(end+1,:) = [2 2 5e-5 1e-4 5e-6 .75 .2 4 26 1 0 Aaper];
case_names{end+1} = 'AIPE_aperiodic_A0.4';
case_par(end+1,:) = [3 2 5e-5 1e-4 5e-6 .75 .2 4 26 1 0 Aaper];
case_names{end+1} = 'CP_same_threshold';
case_par(end+1,:) = [1 1 5e-5 1e-4 5e-6 .75 .2 4 26 1 0 0];
threshold_scales = [.1 .3 3 10];
for s = threshold_scales
    case_names{end+1} = sprintf('CP_grid_x%.1f',s);
    case_par(end+1,:) = [1 1 5e-5 1e-4*s 5e-6*s .75 .2 4 26 1 0 0];
end
for s = threshold_scales
    case_names{end+1} = sprintf('AIPE_grid_x%.1f',s);
    case_par(end+1,:) = [3 2 5e-5 1e-4*s 5e-6*s .75 .2 4 26 1 0 Aaper];
end
fprintf('--- FINAL MAIN A-D BATCH: %d cases, Nmc=%d ---\n',length(case_names),Nmc);
run(fullfile(script_dir,'simulation_engine.m'));
if ~exist('all_metrics','var'), error('simulation_engine.m did not create all_metrics.'); end
save(fullfile(outdir,'final_main_results.mat'),'all_metrics','case_names','case_par', ...
    'curves','curve_time','Nmc','hEM','Tsim','seed','schedule_records', ...
    'state_example','event_example');
if ~exist(fullfile(script_dir,'results'),'dir'), mkdir(fullfile(script_dir,'results')); end
result_prefix = 'final_main';
run(fullfile(script_dir,'write_results.m'));
run(fullfile(script_dir,'sampling_certificate.m'));
movefile(fullfile(script_dir,'results','final_main_summary.csv'),fullfile(outdir,'final_main_summary.csv'),'f');
movefile(fullfile(script_dir,'results','final_main_paths.csv'),fullfile(outdir,'final_main_paths.csv'),'f');
movefile(fullfile(script_dir,'results','final_main_paired.csv'),fullfile(outdir,'final_main_paired.csv'),'f');
movefile(fullfile(script_dir,'results','final_main_curves.csv'),fullfile(outdir,'final_main_curves.csv'),'f');
movefile(fullfile(script_dir,'results','final_main_certificates.csv'),fullfile(outdir,'final_main_certificates.csv'),'f');
fid = fopen(fullfile(outdir,'final_main_cases.csv'),'w');
fprintf(fid,'case,strategy,schedule,hs,sigma1,sigma0,duty,noise,pins,gain,topology,q_h,aperiodic_amplitude\n');
for i=1:length(case_names)
    fprintf(fid,'%s',case_names{i}); fprintf(fid,',%.12g',case_par(i,:)); fprintf(fid,'\n');
end
fclose(fid);
lambda_c=.75; Nc=.05; muTheta=10; Theta_max=muTheta*Nc;
Nstep=round(Tsim/hEM); active=false(1,Nstep); left=0; mcycle=0; boundaries=[];
while left<Nstep
    width=4*round(.05*(1+Aaper*sin(sqrt(2)*mcycle))/(4*hEM));
    width=max(4,width); width=min(width,Nstep-left);
    right=left+width; on_width=round(lambda_c*width);
    active(left+1:left+on_width)=true;
    boundaries=[boundaries; left*hEM (left+on_width)*hEM right*hEM];
    left=right; mcycle=mcycle+1;
end
t=(0:Nstep)*hEM; active_time=[0 cumsum(active)*hEM]; D=lambda_c*t-active_time;
running_min=inf; max_deficit=0;
for k=1:length(D)
    running_min=min(running_min,D(k));
    max_deficit=max(max_deficit,D(k)-running_min);
end
J=-muTheta*D; runmax=zeros(size(J)); tmp=-inf;
for k=1:length(J), tmp=max(tmp,J(k)); runmax(k)=tmp; end
Theta=Theta_max+J-runmax;
fid=fopen(fullfile(outdir,'final_schedule_validation.csv'),'w');
fprintf(fid,'aperiodic_amplitude,active_fraction,max_deficit_difference,Nc,min_Theta,max_Theta,passes_Nc\n');
fprintf(fid,'%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%d\n',Aaper,mean(active),max_deficit,Nc,min(Theta),max(Theta),max_deficit<=Nc+1e-14);
fclose(fid);
dlmwrite(fullfile(outdir,'final_schedule_boundaries.csv'),boundaries,'delimiter',',','precision',12);
clear all_metrics curves curve_time schedule_records state_example event_example;
Tsim = 5;
q_list = [0 .02 .05 .1 .2];
case_names = cell(1,length(q_list));
case_par = zeros(length(q_list),12);
for j=1:length(q_list)
    q=q_list(j);
    case_names{j}=sprintf('bias_q%.2f',q);
    case_par(j,:)=[3 2 5e-5 1e-4 5e-6 .75 .2 4 26 1 q Aaper];
end
fprintf('--- FINAL BIAS BATCH: %d cases, Nmc=%d ---\n',length(case_names),Nmc);
run(fullfile(script_dir,'simulation_engine.m'));
if ~exist('all_metrics','var'), error('simulation_engine.m did not create all_metrics in bias batch.'); end
save(fullfile(outdir,'final_bias_results.mat'),'all_metrics','case_names','case_par', ...
    'curves','curve_time','Nmc','hEM','Tsim','seed','schedule_records', ...
    'state_example','event_example');
result_prefix='final_bias';
run(fullfile(script_dir,'write_results.m'));
run(fullfile(script_dir,'sampling_certificate.m'));
movefile(fullfile(script_dir,'results','final_bias_summary.csv'),fullfile(outdir,'final_bias_summary.csv'),'f');
movefile(fullfile(script_dir,'results','final_bias_paths.csv'),fullfile(outdir,'final_bias_paths.csv'),'f');
movefile(fullfile(script_dir,'results','final_bias_paired.csv'),fullfile(outdir,'final_bias_paired.csv'),'f');
movefile(fullfile(script_dir,'results','final_bias_curves.csv'),fullfile(outdir,'final_bias_curves.csv'),'f');
movefile(fullfile(script_dir,'results','final_bias_certificates.csv'),fullfile(outdir,'final_bias_certificates.csv'),'f');
fid=fopen(fullfile(outdir,'final_bias_cases.csv'),'w');
fprintf(fid,'case,strategy,schedule,hs,sigma1,sigma0,duty,noise,pins,gain,topology,q_h,aperiodic_amplitude\n');
for i=1:length(case_names)
    fprintf(fid,'%s',case_names{i}); fprintf(fid,',%.12g',case_par(i,:)); fprintf(fid,'\n');
end
fclose(fid);
disp('FINAL_AD_FINISHED');
disp('Final A-D results were written to ../data.');
