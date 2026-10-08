if Nmc==100
    tcrit=1.98421695150868;
elseif Nmc==20
    tcrit=2.09302405440826;
elseif Nmc==10
    tcrit=2.26215716274099;
else
    tcrit=1.96;
end
metric_names={'updates','J_error','E_input','tail_error','final_error','max_error','pin_assignments','active_fraction','sensing_instants'};
fid=fopen(fullfile('results',[result_prefix '_summary.csv']),'w');
fprintf(fid,'case,metric,mean,sd,ci_low,ci_high,Nmc\n');
fidraw=fopen(fullfile('results',[result_prefix '_paths.csv']),'w');
fprintf(fidraw,'case,path,updates,J_error,E_input,tail_error,final_error,max_error,pin_assignments,active_fraction,sensing_instants\n');
for icase=1:length(case_names)
    for j=1:9
        values=all_metrics(:,j,icase);
        mu=mean(values);
        sd=std(values);
        hw=tcrit*sd/sqrt(Nmc);
        fprintf(fid,'%s,%s,%.12g,%.12g,%.12g,%.12g,%d\n', ...
            case_names{icase},metric_names{j},mu,sd,mu-hw,mu+hw,Nmc);
    end
    for mc=1:Nmc
        fprintf(fidraw,'%s,%d',case_names{icase},mc);
        fprintf(fidraw,',%.12g',all_metrics(mc,:,icase));
        fprintf(fidraw,'\n');
    end
end
fclose(fid);
fclose(fidraw);
fid=fopen(fullfile('results',[result_prefix '_paired.csv']),'w');
fprintf(fid,'case_a,case_b,metric,mean_a_minus_b,sd_difference,ci_low,ci_high,Nmc\n');
for ia=1:length(case_names)
    for ib=ia+1:length(case_names)
        for j=[1 2 3 4]
            delta=all_metrics(:,j,ia)-all_metrics(:,j,ib);
            mu=mean(delta);
            sd=std(delta);
            hw=tcrit*sd/sqrt(Nmc);
            fprintf(fid,'%s,%s,%s,%.12g,%.12g,%.12g,%.12g,%d\n', ...
                case_names{ia},case_names{ib},metric_names{j},mu,sd,mu-hw,mu+hw,Nmc);
        end
    end
end
fclose(fid);
dlmwrite(fullfile('results',[result_prefix '_curves.csv']),[curve_time' curves'],'precision','%.12g');
