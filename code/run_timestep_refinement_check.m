clear; clc; close all;
script_dir=fileparts(mfilename('fullpath')); data_dir=fullfile(fileparts(script_dir),'data');
Nmc=100; hEM=1e-5; Tsim=1.5; seed=1; Aaper=.4;
base=[3 2 5e-5 1e-4 5e-6 .75 .2 4 26 1 0 .4];
S=load(fullfile(data_dir,'performance_match_workspace.mat'));
cp_metrics=S.cp_metrics; ai_metrics=S.ai_metrics; selected_lambda=S.selected_lambda;
pair_metrics=cat(3,cp_metrics,ai_metrics);
C=load(fullfile(data_dir,'performance_match_certificate.mat'),'boundaries');
boundaries=C.boundaries;
case_par=[1 1 base(3:end);base]; case_par(2,6)=selected_lambda; case_names={'CP-PETC','Aperiodic AIPE-TC'}; schedule_boundaries_override=boundaries;
h_fine=5e-6; h_coarse=1e-5; rng(seed,'twister');
dB_coarse=sqrt(h_coarse)*randn(round(Tsim/h_coarse),Nmc)';
bridge=sqrt(h_coarse/4)*randn(round(Tsim/h_coarse),Nmc)';
dB_fine=zeros(Nmc,round(Tsim/h_fine));
dB_fine(:,1:2:end)=.5*dB_coarse+bridge;
dB_fine(:,2:2:end)=.5*dB_coarse-bridge;
bridge_max_abs_error=max(abs(dB_fine(:,1:2:end)+dB_fine(:,2:2:end)-dB_coarse),[],'all');
assert(bridge_max_abs_error<1e-12,'STOP: Brownian bridge aggregation failed.');
hEM=h_fine; dB_all_override=dB_fine; run(fullfile(script_dir,'simulation_engine.m')); fine_metrics=all_metrics;
hEM=h_coarse; dB_all_override=dB_coarse; run(fullfile(script_dir,'simulation_engine.m')); coarse_metrics=all_metrics;
coarse_identity_max_abs=max(abs(coarse_metrics(:,1:3,:)-pair_metrics(:,1:3,:)),[],'all');
assert(coarse_identity_max_abs<1e-12,'STOP: Table9/Table10 pathwise identity failed.');
combined=cat(3,coarse_metrics(:,:,1),fine_metrics(:,:,1),coarse_metrics(:,:,2),fine_metrics(:,:,2)); method=repelem(["CP-PETC";"Aperiodic AIPE-TC"],2); step=repmat([h_coarse;h_fine],2,1);
writetable([table(method,step) stats_table(combined)],fullfile(data_dir,'timestep_refinement_matched_summary.csv'));
pd=zeros(6,3); relative_change=zeros(6,1); for m=1:2, for j=1:3, ii=(m-1)*3+j; pd(ii,:)=paired_ci(fine_metrics(:,j,m)-coarse_metrics(:,j,m)); relative_change(ii)=pd(ii,1)/mean(coarse_metrics(:,j,m)); end, end
method=repelem(["CP-PETC";"Aperiodic AIPE-TC"],3); metric=repmat(["Delta_Nu";"Delta_Jrho";"Delta_Ju"],2,1);
writetable(table(method,metric,pd(:,1),pd(:,2),pd(:,3),relative_change,'VariableNames',{'method','metric','mean_difference','ci_low','ci_high','relative_change'}),fullfile(data_dir,'timestep_refinement_matched_paired.csv'));
fine_J=squeeze(mean(fine_metrics(:,2,:),1)); coarse_J=squeeze(mean(coarse_metrics(:,2,:),1)); gap_fine=abs(diff(fine_J))/mean(fine_J); gap_coarse=abs(diff(coarse_J))/mean(coarse_J);
save(fullfile(data_dir,'timestep_refinement_matched_workspace.mat'),'coarse_metrics','fine_metrics','boundaries','selected_lambda','gap_fine','gap_coarse','coarse_identity_max_abs','bridge_max_abs_error');
assert(gap_fine<=.01,'STOP: fine-grid matched gap exceeds 1 percent.');
assert(sign(diff(fine_J))==sign(diff(coarse_J)),'STOP: error ranking reversal.');
assert(sign(mean(fine_metrics(:,3,2))-mean(fine_metrics(:,3,1)))==sign(mean(coarse_metrics(:,3,2))-mean(coarse_metrics(:,3,1))),'STOP: effort ranking reversal.');
fprintf('MATCHED_TIMESTEP_DONE coarse_gap=%.9g fine_gap=%.9g\n',gap_coarse,gap_fine);

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
