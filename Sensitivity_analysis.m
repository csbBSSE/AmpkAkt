%% Nullclines for AMPK and AKT with parameters sets
file_name='/path/to/bistable/parameter/set/file';

data_path=fileparts(file_name);
cur_dir=pwd;
cd(data_path);
if (~isfolder('Sensitivity_Analysis'))
    mkdir Sensitivity_Analysis;
end
cd Sensitivity_Analysis;
data_dir=pwd;
cd(cur_dir);
prefix_name=string(split(data_path,filesep));
prefix_name=char(prefix_name(end));

data = readtable(file_name);
par_sets = table2cell(data);
par_set = cell2mat(par_sets(:,2:25));
col_names=[data.Properties.VariableNames(2:25)]';
disp_name = {'total AMPK';'total AKT';'total PHLPP2';'total PP2Ca';...
    'k_{ac} AMPK';'k_{ac} AKT';'k_{ac} PHLPP2';'k_{ac} PP2Ca';...
    'k_{dac} AMPK';'k_{dac} AKT';'k_{dac} PHLPP2';'k_{dac} PP2Ca';...
    'l _{PP2Ca}';'l _{PHLPP2}';'l _{AKT}';'l _{AMPK}';...
    'n _{PP2Ca}';'n _{PHLPP2}';'n _{AKT}';'n _{AMPK}';...
    'PP2Ca^0 _{AMPK}';'PHLPP2^0 _{AKT}';'AMPK^0 _{PHLPP2}';'AKT^0 _{PP2Ca}'};

for m = 1:size(par_set,1)
    cd(data_dir);
    dir_n=['parameter_set_',char(string(par_sets(m,1)))];
    if (~isfolder(eval('dir_n')))
        feval('mkdir',eval('dir_n'));
    end
    feval('cd',eval('dir_n'));
    data_dir_name=pwd;
    cd(cur_dir);
    
    var_pars=repmat(par_set(m,:),size(par_set,2)*2,1);
    par_names=cell(size(par_set,2)*2,1);
    disp_names=cell(size(par_set,2)*2,1);
    for i = 1:size(var_pars,2)
        var_pars(i,i) = var_pars(i,i)*1.1;
        par_names(i)={[char(string(par_sets(m,1))),'_',char(string(col_names(i))),'-10_percent_up']};
        disp_names(i)={[char(string(par_sets(m,1))),' , ',char(string(disp_name(i))),' - 10 % up']};
        var_pars(i+size(var_pars,2),i) = var_pars(i+size(var_pars,2),i)*0.9;
        par_names(i+size(var_pars,2))={[char(string(par_sets(m,1))),'_',char(string(col_names(i))),'-10_percent_down']};
        disp_names(i+size(var_pars,2))={[char(string(par_sets(m,1))),' , ',char(string(disp_name(i))),' - 10 % down']};
    end
    var_pars= [par_set(m,:);var_pars];
    par_name=[[char(string(par_sets(m,1))),'_(core)'];[par_names]];
    disp_names=[[char(string(par_sets(m,1))),' (core)'];[disp_names]];
    var_par =[[str2num(string(par_sets(m,1))),1:size(var_pars,1)-1]',var_pars];
    
    parameter_set=zeros(size(var_par)); %stores the parameter set
    state_vector=zeros(1,9); % stores all the initial and final states of the simulations
    
    parameter_names=strsplit('set_no,t_AMPK,t_AKT,t_PHLPP2,t_PP2Ca,kac_AMPK,kac_AKT,kac_PHLPP2,kac_PP2Ca,kdac_AMPK,kdac_AKT,kdac_PHLPP2,kdac_PP2Ca,l_PP2Ca,l_PHLPP2,l_AKT,l_AMPK,n_PP2Ca,n_PHLPP2,n_AKT,n_AMPK,PP2Ca_0_AMPK,PHLPP2_0_AKT,AMPK_0_PHLPP2,AKT_0_PP2Ca',',');
    var_name=strsplit('set_no,i_AMPK,i_AKT,i_PHLPP2,i_PP2Ca,f_AMPK,f_AKT,f_PHLPP2,f_PP2Ca',',')';
    final_states=cell(size(var_par,1),4);
    
    for i = 1:size(var_par,1)
        %this code generates plot of 100 initial conditions for a given set of parameters
        stress = 0.5;
        
        t_AMPK = var_par(i,2); % total AMPK
        t_AKT = var_par(i,3); % total AKT
        t_PHLPP2 = var_par(i,4); % total PHLPP2
        t_PP2Ca = var_par(i,5); % total PP2Ca
        
        % Activation rate: (per min)
        kac_AMPK = var_par(i,6); % background activation constant of AMPK
        kac_AKT = var_par(i,7); % background activation constant of AKT
        %kac_AKT=kac_AMPK;
        kac_PHLPP2 = var_par(i,8); % background activation constant of PHLPP2
        kac_PP2Ca = var_par(i,9); % background activation constant of PP2Ca
        %kac_PP2Ca = kac_PHLPP2;
        %Jcp = 1; % Michaelis-constant of PHLPP2, PP2Ca
        
        % Inactivation constants: (per min)
        kdac_AMPK = var_par(i,10); % background inactivation constant of AMPK
        kdac_AKT = var_par(i,11); % background inactivation constant of AKT
        %kdac_AKT=kdac_AMPK;
        kdac_PHLPP2 = var_par(i,12); % background inactivation constant of PHLPP2
        kdac_PP2Ca = var_par(i,13); % background inactivation constant of PP2Ca
        %kdac_PP2Ca = kdac_PHLPP2;
        
        % Lambdas
        l_PP2Ca = var_par(i,14);
        l_PHLPP2 = var_par(i,15);
        l_AKT = var_par(i,16);
        l_AMPK = var_par(i,17);
        
        % co-operativity
        n_PP2Ca = round(var_par(i,18));
        n_PHLPP2 = round(var_par(i,19));
        n_AKT = round(var_par(i,20));
        n_AMPK = round(var_par(i,21));
        
        %Thresholds
        PP2Ca_0_AMPK = var_par(i,22);
        PHLPP2_0_AKT = var_par(i,23);
        AMPK_0_PHLPP2 = var_par(i,24);
        AKT_0_PP2Ca = var_par(i,25);
        
        %parameters = [t_AMPK,t_AKT,t_PHLPP2,t_PP2Ca,kac_AMPK,kac_AKT,kac_PHLPP2,kac_PP2Ca,kdac_AMPK,kdac_AKT,kdac_PHLPP2,kdac_PP2Ca,l_PP2Ca,l_PHLPP2,l_AKT,l_AMPK,n_PP2Ca,n_PHLPP2,n_AKT,n_AMPK,PP2Ca_0_AMPK,PHLPP2_0_AKT,AMPK_0_PHLPP2,AKT_0_PP2Ca]';
        tspan = 0:0.1:100;
        % [AMPK, AKT, PHLPP2, PP2Ca]
        %SV = [0 10 1 1]; % initial Concentrations
        
        Initial=zeros(100,4);
        Final=zeros(100,4);
        
        %% solving the equations
        for n = 1:1:100
            SV = [rand(1)*t_AMPK,rand(1)*t_AKT,rand(1)*t_PHLPP2,rand(1)*t_PP2Ca];
            Initial(n,:) = SV;
            [t,x] = ode23(@(t,SV) ampk_akt_dynamics_equations(t,SV,stress,t_AMPK,t_AKT,t_PHLPP2,t_PP2Ca,kac_AMPK,kac_AKT,kac_PHLPP2,kac_PP2Ca,kdac_AMPK,kdac_AKT,kdac_PHLPP2,kdac_PP2Ca,l_PP2Ca,l_PHLPP2,l_AKT,l_AMPK,n_PP2Ca,n_PHLPP2,n_AKT,n_AMPK,PP2Ca_0_AMPK,PHLPP2_0_AKT,AMPK_0_PHLPP2,AKT_0_PP2Ca),tspan,SV);
            Final(n,:) = x(end,:);
        end
        
        
        parameters = [var_par(i,1),t_AMPK,t_AKT,t_PHLPP2,t_PP2Ca,kac_AMPK,kac_AKT,kac_PHLPP2,kac_PP2Ca,kdac_AMPK,kdac_AKT,kdac_PHLPP2,kdac_PP2Ca,l_PP2Ca,l_PHLPP2,l_AKT,l_AMPK,n_PP2Ca,n_PHLPP2,n_AKT,n_AMPK,PP2Ca_0_AMPK,PHLPP2_0_AKT,AMPK_0_PHLPP2,AKT_0_PP2Ca];
        %     par_file=array2table(parameters,'RowNames',parameter_names);
        %     par_f=[data_dir,'/parameter_set_',num2str(j),'.txt'];
        %writetable(par_file,par_f,'WriteRowNames',true,'Delimiter','\t');
        
        parameter_set(i,:) = parameters; % updating Parameter set list
        
        set_no=ones(size(Initial,1),1)*j;
        test_initial=[set_no,Initial,Final];
        
        state_vector=[state_vector;test_initial]; % updating state vector list
        
        %     test_initial=array2table(test_initial,'VariableNames',var_name);
        %     test_f=[data_dir,'/SV_RIS_(100)_parameter_set_',num2str(j),'.txt'];
        %writetable(test_initial,test_f,'WriteVariableNames',true,'Delimiter','\t');
        final_states{i,1} = strjoin(string(unique(round(Final(:,1)))),',');
        final_states{i,2} = strjoin(string(unique(round(Final(:,2)))),',');
        final_states{i,3} = strjoin(string(unique(round(Final(:,3)))),',');
        final_states{i,4} = strjoin(string(unique(round(Final(:,4)))),',');
 
%% sensitivity analysis
diff = zeros(size(var_par,1),1); % stores the diffence of steady states
per = zeros(size(var_par,1),1); % stores the percentage change between the sets
f_states=string.empty;
for i = 1:size(diff,1)
    a=split(string(final_states(i,1)),',');
    a=cellfun(@str2num,a);
    %a=str2num(cell2mat(a));
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
difs=[repmat(diff(1,1),1,2);reshape(diff(2:end),[],2)]; % reshalpe to fit 2 columns with + 10% and - 10% percentage changes
pers=[repmat(per(1,1),1,2);per_d];
f_state=[repmat(f_states(1,1),1,2);reshape(f_states(2:end),[],2)];
col_names_fc_data={'Parameters',['Final_AMPK_par_10_percent_up'],['Final_AMPK_par_10_percent_Down'],...
    ['Diff_States_AMPK_par_10_percent_up'],['Diff_States_AMPK_par_10_percent_Down'],...
    ['Diff_States_percent_AMPK_par_10_percent_up'],['Diff_States_percent_AMPK_par_10_percent_Down']};
parameters=[[char(string(par_sets(m,1))),'_core'];col_names];
dif_par=table(parameters,f_state(:,1),f_state(:,2),difs(:,1),difs(:,2),pers(:,1),pers(:,2),'VariableNames',col_names_fc_data);

dif_par_f=[data_dir_name,filesep,'parameter_set_',char(string(par_sets(m,1))),'_Sensitivity_analysis_plot_data_',prefix_name,'.txt'];
dif_par_fc=[data_dir_name,filesep,'parameter_set_',char(string(par_sets(m,1))),'_Sensitivity_analysis_plot_data_',prefix_name,'.csv'];
writetable(dif_par,dif_par_f,'WriteVariableNames',true,'Delimiter','\t');
writetable(dif_par,dif_par_fc,'WriteVariableNames',true,'Delimiter',',');


pars=array2table(var_par(:,2:end),'VariableNames',parameter_names(2:end));
par_nm=table(char(string(par_name)),'VariableNames',parameter_names(1));
fc_data=table([f_states]',diff,per,'VariableNames',{'Final_State_AMPK','Delta_Final_State_AMPK','Percentage_difference_from_core'});
var_par_set=[par_nm,pars,finals,fc_data];
var_par_f=[data_dir_name,filesep,'parameter_set_',char(string(par_sets(m,1))),'_Sensitivity_analysis_',prefix_name,'.txt'];
var_par_fc=[data_dir_name,filesep,'parameter_set_',char(string(par_sets(m,1))),'_Sensitivity_analysis_',prefix_name,'.csv'];
writetable(var_par_set,var_par_f,'WriteVariableNames',true,'Delimiter','\t');
writetable(var_par_set,var_par_fc,'WriteVariableNames',true,'Delimiter',',');

%% plotting sensetivity analysis
fig1=figure();
set(fig1,'PaperUnits','inches');
set(fig1,'PaperPosition',[0 0 8 10]);
barh(categorical(disp_name),per_d);
ax = gca;
ax.FontSize = 12;
xlabel('% change between AMPK Steady states','FontSize', 12);
xlim([min(per)-10 max(per)+10]);
legend(' 10% Up',' 10% Down','Location','southoutside','Orientation','horizontal','FontSize', 12)
title({'Parameter Fold change Sensitivity Analysis',...
    ['parameter set : ',char(string(par_sets(m,1)))]},'FontSize', 14);
fig_name1=[data_dir_name,filesep,'parameter_set_',char(string(par_sets(m,1))),'_Sensitivity_analysis_plot_',prefix_name];
saveas(fig1,fig_name1,'jpeg');
close;
end