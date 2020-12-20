%{
This script 10000 generates random parameter sets, simulates them for 1000 random
initial conditions and gets saves the parameters and steady states. 
The ODE system is defined in the file "AMPK_AKT_dynamics.m"
%}
%% this code is used to generate 10000 random parameter sets and test the final steady states of the system for 1000 random initial conditions for each set
nps=10000; %number of perameter sets
nis=1000; % initial state vectors
%total concentration of molecules in cells
total_molecules = 100; 
%range of kac and kdac
k_ac_range = [0.02 0.2];
k_deac_range = [0.02 0.2];

%Thresholds ranges
threshold_range = [0.25 0.75]; 

% Lambdas ranges
l_range = [5 10];

% co-operativity 
n_range = [4 6];

%RSPS - Random Symmetric Parameter Sets
%RIS - Random Initial States
formatOut = 'dd_mm_yyyy';
sim_date=datestr(now,formatOut);
cur_dir=pwd;
dir_name=['RPS_',num2str(nps),'_RIS_',num2str(nis),'_(',sim_date,')'];
if (~isfolder(eval('dir_name')))
    feval('mkdir',eval('dir_name'));
end
feval('cd',eval('dir_name'));
data_dir=pwd;
cd(cur_dir);


parameter_names=strsplit('set_no,t_AMPK,t_AKT,t_PHLPP2,t_PP2Ca,kac_AMPK,kac_AKT,kac_PHLPP2,kac_PP2Ca,kdac_AMPK,kdac_AKT,kdac_PHLPP2,kdac_PP2Ca,l_PP2Ca,l_PHLPP2,l_AKT,l_AMPK,n_PP2Ca,n_PHLPP2,n_AKT,n_AMPK,PP2Ca_0_AMPK,PHLPP2_0_AKT,AMPK_0_PHLPP2,AKT_0_PP2Ca',',');
var_name=strsplit('set_no,i_AMPK,i_AKT,i_PHLPP2,i_PP2Ca,f_AMPK,f_AKT,f_PHLPP2,f_PP2Ca',',')';    

final_states=string(zeros(nps,4));

parameter_set=zeros(nps,25); %stores the parameter set
state_vector=zeros(1,9); % stores all the initial and final states of the simulations

sel_par=[];

for j = 1:1:nps
%% randamization of parameters
    %Total molecules
    t_AMPK = total_molecules; % total AMPK
    t_AKT = total_molecules; % total AKT
    t_PHLPP2 = total_molecules; % total PHLPP2
    t_PP2Ca = total_molecules; % total PP2Ca

    % Activation rate: (per min)
    kac_AMPK = (rand(1) * (k_ac_range(2) - k_ac_range(1))) + k_ac_range(1); % background activation constant of AMPK
    kac_AKT = (rand(1) * (k_ac_range(2) - k_ac_range(1))) + k_ac_range(1); % background activation constant of AKT
    %kac_AKT=kac_AMPK;
    kac_PHLPP2 = (rand(1) * (k_ac_range(2) - k_ac_range(1))) + k_ac_range(1); % background activation constant of PHLPP2
    kac_PP2Ca = (rand(1) * (k_ac_range(2) - k_ac_range(1))) + k_ac_range(1); % background activation constant of PP2Ca
    %kac_PP2Ca = kac_PHLPP2;

    % Inactivation constants: (per min)
    kdac_AMPK = (rand(1) * (k_deac_range(2) - k_deac_range(1))) + k_deac_range(1); % background inactivation constant of AMPK
    kdac_AKT = (rand(1) * (k_deac_range(2) - k_deac_range(1))) + k_deac_range(1); % background inactivation constant of AKT
    %kdac_AKT=kdac_AMPK;
    kdac_PHLPP2 = (rand(1) * (k_deac_range(2) - k_deac_range(1))) + k_deac_range(1); % background inactivation constant of PHLPP2
    kdac_PP2Ca = (rand(1) * (k_deac_range(2) - k_deac_range(1))) + k_deac_range(1); % background inactivation constant of PP2Ca
    %kdac_PP2Ca = kdac_PHLPP2;
	
    %Thresholds
    PP2Ca_0_AMPK = ((rand(1)* (threshold_range(2) - threshold_range(1))) + threshold_range(1))*t_PP2Ca; 
    PHLPP2_0_AKT = ((rand(1)* (threshold_range(2) - threshold_range(1))) + threshold_range(1))*t_PHLPP2;
    %PHLPP2_0_AKT = PP2Ca_0_AMPK;
    AMPK_0_PHLPP2 = ((rand(1) * (threshold_range(2) - threshold_range(1))) + threshold_range(1))*t_AMPK;
    AKT_0_PP2Ca = ((rand(1) * (threshold_range(2) - threshold_range(1))) + threshold_range(1))*t_AKT;
    %AKT_0_PP2Ca = AMPK_0_PHLPP2;

    % Lambdas
    l_PP2Ca = (rand(1) * (l_range(2) - l_range(1))) + l_range(1);
    l_PHLPP2 = (rand(1) * (l_range(2) - l_range(1))) + l_range(1);
    %l_PHLPP2 = l_PP2Ca;
    l_AMPK = (rand(1) * (l_range(2) - l_range(1))) + l_range(1);
    l_AKT = (rand(1) * (l_range(2) - l_range(1))) + l_range(1);
    %l_AKT = l_AMPK;

    % co-operativity
    n_PP2Ca = round((rand(1) * (n_range(2) - n_range(1)))) + n_range(1);
    n_PHLPP2 = round((rand(1) * (n_range(2) - n_range(1)))) + n_range(1);
    %n_PHLPP2 = n_PP2Ca;
    n_AMPK = round((rand(1) * (n_range(2) - n_range(1)))) + n_range(1);
    n_AKT = round((rand(1) * (n_range(2) - n_range(1)))) + n_range(1);
    %n_AKT = n_AMPK;

    tspan = 0:0.1:1000;
    % [AMPK, AKT, PHLPP2, PP2Ca]
    %SV = [0 10 1 1]; % initial Concentrations

    Initial=zeros(nis,4);
    Final=zeros(nis,4);
    %% solving the equations
    parfor n = 1:1:nis
        SV = [rand(1,4)*total_molecules];
        Initial(n,:) = SV;
        [t,x] = ode23(@(t,SV) AMPK_AKT_dynamics(t,SV,t_AMPK,t_AKT,t_PHLPP2,t_PP2Ca,kac_AMPK,kac_AKT,kac_PHLPP2,kac_PP2Ca,kdac_AMPK,kdac_AKT,kdac_PHLPP2,kdac_PP2Ca,l_PP2Ca,l_PHLPP2,l_AKT,l_AMPK,n_PP2Ca,n_PHLPP2,n_AKT,n_AMPK,PP2Ca_0_AMPK,PHLPP2_0_AKT,AMPK_0_PHLPP2,AKT_0_PP2Ca),tspan,SV);
        Final(n,:) = x(end,:);
    end
    
    %% writing current parameter set simulation data into files
    %cd(data_dir);
    
    parameters = real([j,t_AMPK,t_AKT,t_PHLPP2,t_PP2Ca,kac_AMPK,kac_AKT,kac_PHLPP2,kac_PP2Ca,kdac_AMPK,kdac_AKT,kdac_PHLPP2,kdac_PP2Ca,l_PP2Ca,l_PHLPP2,l_AKT,l_AMPK,n_PP2Ca,n_PHLPP2,n_AKT,n_AMPK,PP2Ca_0_AMPK,PHLPP2_0_AKT,AMPK_0_PHLPP2,AKT_0_PP2Ca]');
    parameter_set(j,:) = parameters;
    par_file=array2table(parameters,'RowNames',parameter_names);
    par_f=[data_dir,'/parameter_set_',num2str(j),'.txt'];
    writetable(par_file,par_f,'WriteRowNames',true,'Delimiter','\t');   
    
    %parameter_set = [parameter_set; parameters']; % updating Parameter set list
    
    set_no=real(ones(size(Initial,1),1)*j);
    test_initial=[set_no,Initial,Final];
    
    state_vector=[state_vector;test_initial]; % updating state vector list
    
    test_initial=array2table(test_initial,'VariableNames',var_name);
    test_f=[data_dir,'/SV_RIS_(',num2str(nis),')_parameter_set_',num2str(j),'.txt']; 
    writetable(test_initial,test_f,'WriteVariableNames',true,'Delimiter','\t');
        
    %cd(cur_dir);
    final_states(j,1) = strjoin(string(unique(round(Final(:,1)))),',');
    final_states(j,2) = strjoin(string(unique(round(Final(:,2)))),',');
    final_states(j,3) = strjoin(string(unique(round(Final(:,3)))),',');
    final_states(j,4) = strjoin(string(unique(round(Final(:,4)))),',');
    if (length(unique(round(Final(:,1))))>1)
        sel_par = [sel_par,j];
    end
    
end
%% saving all parameter sets and state vectors
    par_set_file=array2table(parameter_set,'VariableNames',parameter_names');
    par_set_filename=[data_dir,'/parameter_sets_',dir_name,'.txt'];
    writetable(par_set_file, par_set_filename,'WriteVariableNames',true,'Delimiter','\t');
    
    
    par_sets=[parameter_set,final_states];
    a=[string([parameter_names]')]';
    col_names=[[string([parameter_names]')]','final_AMPK','final_AKT','final_PHLPP2','final_PP2Ca'];
    par_sets_file=array2table(par_sets);
    par_sets_file.Properties.VariableNames = col_names;
    par_sets_filename=[data_dir,'/parameter_sets_and_final_states_',dir_name,'.txt'];
    writetable(par_sets_file, par_sets_filename,'WriteVariableNames',true,'Delimiter','\t');
    
    sel_sets=par_sets(sel_par,:);
    sel_2=[];
    for i = 1:(size(sel_sets,1))
       a=split(sel_sets(i,26),',');
       b=zeros(size(a));
       for j =1:length(a)
           b(j)=str2num(a(j));
       end
       % a=str2num(cell2mat(a));
       if((b(end) - b(1))>=20)
           sel_2=[sel_2,i];
       end
    end
    sel_2_sets=sel_sets(sel_2,:);
    sel_2_file=array2table(sel_2_sets);
    sel_2_file.Properties.VariableNames = col_names;
    sel_2_sets_filename=[data_dir,'/parameter_sets_with_2_final_states_',dir_name,'.txt'];
    writetable(sel_2_file, sel_2_sets_filename,'WriteVariableNames',true,'Delimiter','\t');
    
    state_vector(1,:)=[];
    state_vector_file=array2table(state_vector,'VariableNames',var_name);
    state_vector_filename=[data_dir,'/state_vectors_',dir_name,'.txt'];
    writetable(state_vector_file, state_vector_filename,'WriteVariableNames',true,'Delimiter','\t');
     
