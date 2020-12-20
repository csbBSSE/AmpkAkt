%file_name='/media/csb/D/Adithya/AMPK_AKT/RSPS_1000_RIS_100(05_05_2020)/parameter_sets_with_2_final_states_RSPS_1000_RIS_100(05_05_2020).txt';
%file_name='/media/csb/D/Adithya/AMPK_AKT/RPS_10000_RIS_100(11_05_2020)/parameter_sets_with_2_final_states_RPS_10000_RIS_100(11_05_2020).txt';
%'

file_name = '/path/to/parameter/set/file';


data_path=fileparts(file_name);
cur_dir=pwd;
cd(data_path);
if (~isfolder('Noise_perturbation'))
    mkdir Noise_perturbation;
end
cd Noise_perturbation;
data_dir=pwd;
cd(cur_dir);
prefix_name=string(split(data_path,filesep));
prefix_name=char(prefix_name(end));

data = readtable(file_name);
par_sets = table2cell(data);
par_set = cell2mat(par_sets(:,1:25));


for i = 1:size(par_set,1)
    
%this code generates plot of 100 initial conditions for a given set of parameters 

    t_AMPK = par_set(i,2); % total AMPK
    t_AKT = par_set(i,3); % total AKT
    t_PHLPP2 = par_set(i,4); % total PHLPP2
    t_PP2Ca = par_set(i,5); % total PP2Ca

    % Activation rate: (per min)
    kac_AMPK = par_set(i,6); % background activation constant of AMPK
    kac_AKT = par_set(i,7); % background activation constant of AKT
    %kac_AKT=kac_AMPK;
    kac_PHLPP2 = par_set(i,8); % background activation constant of PHLPP2
    kac_PP2Ca = par_set(i,9); % background activation constant of PP2Ca
    %kac_PP2Ca = kac_PHLPP2;
    %Jcp = 1; % Michaelis-constant of PHLPP2, PP2Ca

    % Inactivation constants: (per min)
    kdac_AMPK = par_set(i,10); % background inactivation constant of AMPK
    kdac_AKT = par_set(i,11); % background inactivation constant of AKT
    %kdac_AKT=kdac_AMPK;
    kdac_PHLPP2 = par_set(i,12); % background inactivation constant of PHLPP2
    kdac_PP2Ca = par_set(i,13); % background inactivation constant of PP2Ca
    %kdac_PP2Ca = kdac_PHLPP2;

    % Lambdas
    l_PP2Ca = par_set(i,14);
    l_PHLPP2 = par_set(i,15);
    l_AKT = par_set(i,16);
    l_AMPK = par_set(i,17);
    
    % co-operativity
    n_PP2Ca = par_set(i,18);
    n_PHLPP2 = par_set(i,19);
    n_AKT = par_set(i,20);
    n_AMPK = par_set(i,21);
    
    %Thresholds
    PP2Ca_0_AMPK = par_set(i,22); 
    PHLPP2_0_AKT = par_set(i,23);
    AMPK_0_PHLPP2 = par_set(i,24);
    AKT_0_PP2Ca = par_set(i,25);
    
parameters = [t_AMPK,t_AKT,t_PHLPP2,t_PP2Ca,kac_AMPK,kac_AKT,kac_PHLPP2,kac_PP2Ca,kdac_AMPK,kdac_AKT,kdac_PHLPP2,kdac_PP2Ca,l_PP2Ca,l_PHLPP2,l_AKT,l_AMPK,n_PP2Ca,n_PHLPP2,n_AKT,n_AMPK,PP2Ca_0_AMPK,PHLPP2_0_AKT,AMPK_0_PHLPP2,AKT_0_PP2Ca]';
%   tspan = 0:0.1:200;
    % [AMPK, AKT, PHLPP2, PP2Ca]
    %SV = [0 10 1 1]; % initial Concentrations

    %Initial=zeros(100,4);
    %Final=zeros(100,4);
%     tspan1 = 0:0.1:100;
%     tspan = 100:0.1:1000;
    Tmax = 5000;
    t1 = 100;
    t_init = 0;
    t_step = 0.1;
    %% solving the equations
    for amp = 20:10:40
        %SV = [45,45,70,50];
        SV=[0,0,0,0];
        %SV = [50,50,50,50];
        init = SV;
        t0 = t_init;
        tspan = 0;
        dat = zeros(1,4);
        for counter = 1:Tmax/t1 % counter upto Tmax/t1
            tsp = t0:t_step:(t0 + t1 - t_step); % simulation time for current iteration
            [t,x] = ode23(@(t,SV) AMPK_AKT_dynamics(t,SV,...
                t_AMPK,t_AKT,t_PHLPP2,t_PP2Ca,kac_AMPK,...
                kac_AKT,kac_PHLPP2,kac_PP2Ca,kdac_AMPK,kdac_AKT,...
                kdac_PHLPP2,kdac_PP2Ca,l_PP2Ca,l_PHLPP2,l_AKT,l_AMPK,...
                n_PP2Ca,n_PHLPP2,n_AKT,n_AMPK,PP2Ca_0_AMPK,PHLPP2_0_AKT,...
                AMPK_0_PHLPP2,AKT_0_PP2Ca),tsp,init); % simulation
            tspan = [tspan;t]; % save the simulation time
            dat = [dat;x]; % save the simulated data
            t0 = t0 + t1; % update the starting time for next iteration
            init = x(end,:) + randn(1,4)*amp; % update the initial condition by adding noise
        end
        noise_data=[tspan,dat];
        noise_data(1,:)=[];
        np_file=array2table(noise_data,'VariableNames',{'time','AMPK','AKT','PHLPP2','PP2Ca'});
        np_filename=[data_dir,filesep,'parameter_set_',char(string(par_sets(i,1))),'_noise_(amplitude_',num2str(amp),')_induced_perturbation_data_',prefix_name,'.txt'];
        writetable(np_file,np_filename,'WriteVariableNames',true,'Delimiter','\t');
  
   %%
   fig = figure();
   subplot(2,1,1);
   plot(tspan(1:100:size(tspan,1)), dat(1:100:size(tspan,1),1),'-r','LineWidth',2);
   xlabel('time');
   ylabel('AMPK');
   ylim([0 t_AMPK]);
   yticks(0:20:t_AMPK);
   xlim([0 Tmax]);
    set(gca,'FontSize',16)
   %legend('AMPK');
   hold on;
   subplot(2,1,2);
   plot(tspan(1:100:size(tspan,1)), dat(1:100:size(tspan,1),2),'-g','LineWidth',2);
   %legend('AKT');
   xlabel('time');
   ylabel('AKT');
   ylim([0 t_AKT]);
   xlim([0 Tmax]);
   yticks(0:20:t_AKT);
   sgtitle({'{\bf\fontsize{12}Noise induced Dynamics of AMPK and AKT}',['parameter set : ',char(string(par_sets(i,1))),', Noise amplitude : ',num2str(amp)]},'FontWeight','Normal');
    set(gca,'FontSize',16)
   fig_name=[data_dir,filesep,'parameter_set_',char(string(par_sets(i,1))),'_noise_(amplitude_',num2str(amp),')_induced_Dynamics_of_AMPK_AKT_plot_',prefix_name];
   saveas(fig,fig_name,'jpeg');
   saveas(fig,fig_name,'png');
   savefig(fig_name);
   close;
  
   fig1 = figure();
   plot(tspan(1:100:size(tspan,1)), dat(1:100:size(tspan,1),1),'-r','LineWidth',2);
   hold on;
   plot(tspan(1:100:size(tspan,1)), dat(1:100:size(tspan,1),2),'-g','LineWidth',2);
   xlabel('time');
   ylabel('AMPK, AKT levels');
   ylim([0 t_AKT]);
   xlim([0 Tmax]);
   yticks(0:20:t_AKT);
   set(gca,'FontSize',16)
   legend('AMPK','AKT');
   title('Noise induced Dynamics of AMPK and AKT');
   fig_name1=[data_dir,filesep,'parameter_set_',char(string(par_sets(i,1))),'_Noise_(amplitude_',num2str(amp),')_induced_Dynamics_of_AMPK_and_AKT_plot_',prefix_name];
   saveas(fig1,fig_name1,'jpeg');
   saveas(fig1,fig_name1,'png');
   savefig(fig_name1)
   close;
     end
end