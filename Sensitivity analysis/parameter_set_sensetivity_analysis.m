file_name = '/path/to/bistable/parameter/sets';
par_sets = readtable(file_name);
cur_dir=pwd;
par_sets = table2cell(par_sets);
par_set = cell2mat(par_sets(1:2,2:25));
for m = 1:size(par_set,1)
    var_pars=repmat(par_set(m,:),size(par_set,2)*2,1);
    for i = 1:size(var_pars,2)
        var_pars(i,i) = var_pars(i,i)*1.1;
        var_pars(i+size(var_pars,2),i) = var_pars(i+size(var_pars,2),i)*0.9;
    end
    var_pars= [par_set(m,:);var_pars];
    var_par =[[str2num(string(par_sets(m,1))),1:size(var_pars,1)-1]',var_pars];
    parameter_set=zeros(size(var_par,1),size(var_par,2)); %stores the parameter set
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
        
        
        parameters = [i,t_AMPK,t_AKT,t_PHLPP2,t_PP2Ca,kac_AMPK,kac_AKT,kac_PHLPP2,kac_PP2Ca,kdac_AMPK,kdac_AKT,kdac_PHLPP2,kdac_PP2Ca,l_PP2Ca,l_PHLPP2,l_AKT,l_AMPK,n_PP2Ca,n_PHLPP2,n_AKT,n_AMPK,PP2Ca_0_AMPK,PHLPP2_0_AKT,AMPK_0_PHLPP2,AKT_0_PP2Ca]';
        %     par_file=array2table(parameters,'RowNames',parameter_names);
        %     par_f=[data_dir,'/parameter_set_',num2str(j),'.txt'];
        %writetable(par_file,par_f,'WriteRowNames',true,'Delimiter','\t');
        
        parameter_set = [parameter_set; parameters']; % updating Parameter set list
        
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
    end
    finals=cell2table(final_states,'VariableNames',{'final_AMPK','final_AKT','final_PHLPP2','final_PP2Ca'});
    pars=array2table(var_par,'VariableNames',parameter_names);
    var_par_set=[pars,finals];
    var_par_f=['E:/Adithya_C/RSPS_1000_RIS_100(29_04_2020)_senstivity_analysis_RIS_100_parameter_set_',char(string(par_sets(m,1))),'.txt'];
    var_par_fc=['E:/Adithya_C/RSPS_1000_RIS_100(29_04_2020)_senstivity_analysis_RIS_100_parameter_set_',char(string(par_sets(m,1))),'.csv'];
    writetable(var_par_set,var_par_f,'WriteVariableNames',true,'Delimiter','\t');
    writetable(var_par_set,var_par_fc,'WriteVariableNames',true,'Delimiter',',');
    
end
    
        