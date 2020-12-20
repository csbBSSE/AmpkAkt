%% Nullclines for AMPK and AKT with parameters sets
file_name = '/path/to/bistable/parameter/sets';
data_path=fileparts(file_name);
cur_dir=pwd;
cd(data_path);
if (~isfolder('Nullclines'))
    mkdir Nullclines;
end
cd Nullclines;
data_dir=pwd;
cd(cur_dir);
prefix_name=string(split(data_path,filesep));
prefix_name=char(prefix_name(end));

par_sets = readtable(file_name(m,:));
par_sets = table2cell(par_sets);
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

    % parameters = [t_AMPK,t_AKT,t_PHLPP2,t_PP2Ca,kac_AMPK,kac_AKT,kac_PHLPP2,kac_PP2Ca,kdac_AMPK,kdac_AKT,kdac_PHLPP2,kdac_PP2Ca,l_PP2Ca,l_PHLPP2,l_AKT,l_AMPK,n_PP2Ca,n_PHLPP2,n_AKT,n_AMPK,PP2Ca_0_AMPK,PHLPP2_0_AKT,AMPK_0_PHLPP2,AKT_0_PP2Ca]';
    tspan = 0:0.1:1000;
    % [AMPK, AKT, PHLPP2, PP2Ca]
    %SV = [0 10 1 1]; % initial Concentrations
    %% solving the equations for AMPK with constant AKT
    AKTc = 0:0.1:t_AKT;
    AMPK_data = zeros(length(AKTc),3);
    for n = 1:1:length(AKTc)
        SV = [0,0,0];
        %Initial(n,:) = SV;
        [t,x] = ode23(@(t,SV) akt_nullcline(t,SV,AKTc(n),t_AMPK,t_AKT,t_PHLPP2,t_PP2Ca,kac_AMPK,kac_AKT,kac_PHLPP2,kac_PP2Ca,kdac_AMPK,kdac_AKT,kdac_PHLPP2,kdac_PP2Ca,l_PP2Ca,l_PHLPP2,l_AKT,l_AMPK,n_PP2Ca,n_PHLPP2,n_AKT,n_AMPK,PP2Ca_0_AMPK,PHLPP2_0_AKT,AMPK_0_PHLPP2,AKT_0_PP2Ca),tspan,SV);
        AMPK_data(n,:) = x(end,:);
        %plot(t,x(:,1),'-r',t,x(:,2),'-g')
        %hold on;
    end
    %test_initial=[Initial,Final];
    %% solving the equations for AKT with conctant AMPK
    AMPKc = 0:0.1:t_AMPK;
    AKT_data = zeros(length(AMPKc),3);
    for n = 1:1:length(AMPKc)
        SV = [0,0,0];
        %Initial(n,:) = SV;
        [t,x] = ode23(@(t,SV) ampk_nullcline(t,SV,AMPKc(n),t_AMPK,t_AKT,t_PHLPP2,t_PP2Ca,kac_AMPK,kac_AKT,kac_PHLPP2,kac_PP2Ca,kdac_AMPK,kdac_AKT,kdac_PHLPP2,kdac_PP2Ca,l_PP2Ca,l_PHLPP2,l_AKT,l_AMPK,n_PP2Ca,n_PHLPP2,n_AKT,n_AMPK,PP2Ca_0_AMPK,PHLPP2_0_AKT,AMPK_0_PHLPP2,AKT_0_PP2Ca),tspan,SV);
        AKT_data(n,:) = x(end,:);
        %plot(t,x(:,1),'-r',t,x(:,2),'-g')
        %hold on;
    end
    %test_initial=[Initial,Final];
    %%
    fig = figure();
    plot(AMPK_data(:,1),AKTc,'-r',AMPKc,AKT_data(:,1),'-g');
    xlabel('AMPK');
    ylabel('AKT');
    xlim([0 100]);
    xticks(0:10:par_set(i,2));
    ylim([0 100]);
    yticks(0:10:par_set(i,3));
    %hold on;
    legend('AKT nullclines','AMPK nullclines')%,'PHLPP2','PP2Ca');
    title({'Nullclines of AMPK and AKT',['parameter set : ',char(string(par_sets(i,1)))]});
    fig_name=[data_dir,filesep,'parameter_set_',char(string(par_sets(i,1))),'_Nullclines_of_AMPK_and_AKT_',prefix_name];
    saveas(fig,fig_name,'jpeg')
    saveas(fig,fig_name,'png')
    close;
    %% saving nullclines data
    aktc_data=[AMPK_data(:,1),AKTc',AMPK_data(:,2:3)];
    ampkc_data=[AMPKc',AKT_data];
    aktc_file=array2table(aktc_data,'VariableNames',{'AMPK_data','AKTc','PHLPP2_data','PP2Ca_data'});
    ampkc_file=array2table(ampkc_data,'VariableNames',{'AMPKc','AKT_data','PHLPP2_data','PP2Ca_data'});
    aktc_f=[data_dir,filesep,'parameter_set_',char(string(par_sets(i,1))),'_AKT_nullcline_data_',prefix_name,'.txt'];
    ampkc_f=[data_dir,filesep,'parameter_set_',char(string(par_sets(i,1))),'_AMPK_nullcline_data_',prefix_name,'.txt'];

    writetable(aktc_file,aktc_f,'WriteVariableNames',true,'Delimiter','\t');
    writetable(ampkc_file,ampkc_f,'WriteVariableNames',true,'Delimiter','\t');
end
