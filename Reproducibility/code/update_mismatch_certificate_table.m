function update_mismatch_certificate_table(work_dir)
script_dir=fileparts(mfilename('fullpath'));
data_dir=fullfile(script_dir,'..','data');
table_file=fullfile(script_dir,'..','tables','table_persistent_mismatch.tex');
if nargin<1, work_dir=tempname; end
if ~isfolder(work_dir), mkdir(work_dir); end
C=readtable(fullfile(data_dir,'certificate_bound_comparison.csv'));
S=readtable(fullfile(data_dir,'final_bias_summary.csv'));
q=[0;.02;.05;.10];
Rold=zeros(4,1); Rav=Rold; measured=Rold; ratio=Rold; display=cell(4,1);
txt=fileread(table_file);
for j=1:4
    name=sprintf('bias_q%.2f',q(j));
    ci=find(strcmp(string(C.case_name),name));
    si=find(strcmp(string(S{:,1}),name)&strcmp(string(S.metric),'tail_error'));
    assert(isscalar(ci)&&isscalar(si));
    Rold(j)=C.R_minrate(ci); Rav(j)=C.R_av(ci);
    measured(j)=S.mean(si); ratio(j)=Rav(j)/measured(j);
    display{j}=sprintf('%.6g',Rav(j));
    pattern=['(?m)^(' sprintf('%.2f',q(j)) ' & [^\r\n]* & )[^&\r\n]+( \\\\[^\r\n]*)$'];
    lines=regexp(txt,pattern,'tokens');
    assert(numel(lines)==1,'Expected one mismatch table row.');
    replacement=[lines{1}{1} display{j} lines{1}{2}];
    matched=regexp(txt,pattern,'match');
    txt=strrep(txt,matched{1},replacement);
end
old='Persistent-mismatch results and certified ultimate bounds.';
new=[old ' $R_{\rm cert}=R_{\rm av}$ denotes the average-rate ultimate certificate in Theorem~1.'];
if ~contains(txt,new), txt=strrep(txt,old,new); end
fid=fopen(table_file,'w'); assert(fid>=0); fprintf(fid,'%s',txt); fclose(fid);
rows=table(q,Rold,Rav,measured,ratio,display);
writetable(rows,fullfile(work_dir,'mismatch_certificate_regression.csv'));
baseline=C(strcmp(string(C.case_name),'AIPE_aperiodic_A0.4'),:);
report=struct('q',q,'old',Rold,'Rav',Rav,'measured',measured,'ratio',ratio, ...
    'display',{display},'baseline',table2struct(baseline));
fid=fopen(fullfile(work_dir,'mismatch_certificate_regression.json'),'w');
fprintf(fid,'%s',jsonencode(report)); fclose(fid);
disp(rows);
end
