data = readtable('RSPS_1000_RIS_100(29_04_2020)_senstivity_analysis_parameter_set_446.csv');
par_sets = table2cell(data);

diff = zeros(size(data,1),1); % stores the diffence of steady states
per = zeros(size(data,1),1); % stores the percentage change between the sets
f_states=string.empty;
for i = 1:size(data,1)
    a=split(table2cell(data(i,26)),',');
    a=str2num(cell2mat(a));
    if (length(a)<2)
        diff(i,1)=a;
        f_states(i) = string(a(1));
    else
        diff(i,1) = a(end) - a(1); % calculate the diffence of steady states
        f_states(i) = (strjoin([string(a(1)),string(a(end))],','));
    end
    
    per(i,1) = ((diff(1,1) - diff(i,1))/diff(1,1))*100; % calculate the percentage change between the sets
end
per_d=reshape(per(2:end),[],2); % reshalpe to fit 2 columns with + 10% and - 10% percentage changes
% names of the y-axis
par_names = categorical({'total AMPK';'total AKT';'total PHLPP2';'total PP2Ca';...
    'k_{ac} AMPK';'k_{ac} AKT';'k_{ac} PHLPP2';'k_{ac} PP2Ca';...
    'k_{dac} AMPK';'k_{dac} AKT';'k_{dac} PHLPP2';'k_{dac} PP2Ca';...
    'l _{PP2Ca}';'l _{PHLPP2}';'l _{AKT}';'l _{AMPK}';...
    'n _{PP2Ca}';'n _{PHLPP2}';'n _{AKT}';'n _{AMPK}';...
    'PP2Ca^0 _{AMPK}';'PHLPP2^0 _{AKT}';'AMPK^0 _{PHLPP2}';'AKT^0 _{PP2Ca}'});
%generating bar plot
fig=figure();
set(fig,'PaperUnits','inches');
set(fig,'PaperPosition',[0 0 8 10]);
barh(par_names,per_d);
ax = gca;
ax.FontSize = 12; 
xlabel('% change between AMPK Steady states','FontSize', 12);
xlim([min(per)-10 max(per)+10]);
legend(' 10% Up',' 10% Down','Location','southoutside','Orientation','horizontal','FontSize', 12)
title({'Parameter Sensitivity Analysis',...
    ['parameter set : ',char(string(par_sets(1,1)))]},'FontSize', 14);
fig_name=['RSPS_1000_RIS_100(29_04_2020)_senstivity_analysis_plot_parameter_set_',char(string(par_sets(1,1)))];
saveas(fig,fig_name,'jpeg');
%saving the files
% dift=table(diff,per,'VariableNames',{'delta_AMPK_Final_states','percentage_change_orginal'});
difs=[repmat(diff(1,1),1,2);reshape(diff(2:end),[],2)]; % reshalpe to fit 2 columns with + 10% and - 10% percentage changes
pers=[repmat(per(1,1),1,2);per_d];
f_state=[repmat(f_states(1,1),1,2);reshape(f_states(2:end),[],2)];
par_name=[[char(string(par_sets(1,1))),'_core'],data.Properties.VariableNames(2:25)]';
col_names={'Parameters','Final_AMPK_par_10_percent_up','Final_AMPK_par_10_percent_Down','Diff_States_AMPK_par_10_percent_up','Diff_States_AMPK_par_10_percent_Down','Diff_States_percent_AMPK_par_10_percent_up','Diff_States_percent_AMPK_par_10_percent_Down'};
%data1=[par_name,f_state,difs,pers];
var_par=table(par_name,f_state(:,1),f_state(:,2),difs(:,1),difs(:,2),pers(:,1),pers(:,2),'VariableNames',col_names);
var_par_f=['RSPS_1000_RIS_100(29_04_2020)_senstivity_analysis_plot_data_parameter_set_',char(string(par_sets(1,1))),'.txt'];
var_par_fc=['RSPS_1000_RIS_100(29_04_2020)_senstivity_analysis_plot_data_parameter_set_',char(string(par_sets(1,1))),'.csv'];
writetable(var_par,var_par_f,'WriteVariableNames',true,'Delimiter','\t');
writetable(var_par,var_par_fc,'WriteVariableNames',true,'Delimiter',',');