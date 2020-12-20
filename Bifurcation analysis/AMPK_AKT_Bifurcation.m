function [x,v,s,h,f] = AMPK_AKT_Bifurcation(parameters,ap)

addpath 'C:\Users\Asus\Desktop\Github\Projects\MatCont'; 
curdir = pwd;
init;
cd(curdir);

opt = contset;
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',20000);
opt=contset(opt,'MinStepsize',0.001);
opt=contset(opt,'MaxStepsize',0.01);
opt=contset(opt,'Eigenvalues',1);
%opt=contset(opt,'Backward',1);



    t_AMPK = parameters(1); % total AMPK
    t_AKT = parameters(2); % total AKT
    t_PHLPP2 = parameters(3); % total PHLPP2
    t_PP2Ca = parameters(4); % total PP2Ca

    % Activation rate: (per min)
    kac_AMPK = parameters(5);%par_set(i,6); % background activation constant of AMPK
    kac_AKT = parameters(6); % background activation constant of AKT
    %kac_AKT=kac_AMPK;
    kac_PHLPP2 = parameters(7); % background activation constant of PHLPP2
    kac_PP2Ca = parameters(8); % background activation constant of PP2Ca
    %kac_PP2Ca = kac_PHLPP2;
    %Jcp = 1; % Michaelis-constant of PHLPP2, PP2Ca

    % Inactivation constants: (per min)
    kdac_AMPK = parameters(9); % background inactivation constant of AMPK
    kdac_AKT = parameters(10); % background inactivation constant of AKT
    %kdac_AKT=kdac_AMPK;
    kdac_PHLPP2 = parameters(11); % background inactivation constant of PHLPP2
    kdac_PP2Ca = parameters(12); % background inactivation constant of PP2Ca
    %kdac_PP2Ca = kdac_PHLPP2;

    % Lambdas
    l_PP2Ca = parameters(13);
    l_PHLPP2 = parameters(14);
    l_AKT = parameters(15);
    l_AMPK = parameters(16);
    
    % co-operativity
    n_PP2Ca = parameters(17);
    n_PHLPP2 = parameters(18);
    n_AKT = parameters(19);
    n_AMPK = parameters(20);
    
    %Thresholds
    PP2Ca_0_AMPK = parameters(21); 
    PHLPP2_0_AKT = parameters(22);
    AMPK_0_PHLPP2 = parameters(23);
    AKT_0_PP2Ca = parameters(24);

% all the notation used here is from first to second, that the value
% corresponds to the effect of the first gene on the second


%ap = 6; %describes the index of parameter for which the bifurcation is drawn using the init_EP_EP function. 

handles = feval(@AMPK_AKT);
tspan = 0:100:50000;

% initial condition
x_start = [0 0 0 0 ];

%calculating steady state for given initial condition 
[t,x_time] = ode23(@(t,SV)handles{2}(t,SV,t_AMPK,t_AKT,t_PHLPP2,t_PP2Ca,kac_AMPK,kac_AKT,kac_PHLPP2,kac_PP2Ca,kdac_AMPK,kdac_AKT,kdac_PHLPP2,kdac_PP2Ca,l_PP2Ca,l_PHLPP2,l_AKT,l_AMPK,n_PP2Ca,n_PHLPP2,n_AKT,n_AMPK,PP2Ca_0_AMPK,PHLPP2_0_AKT,AMPK_0_PHLPP2,AKT_0_PP2Ca),tspan,x_start);
x_init = x_time(end,:)';

%drawing bifurcation using a continuation method
[x0,v0] = init_EP_EP(@AMPK_AKT,x_init,[t_AMPK,t_AKT,t_PHLPP2,t_PP2Ca,kac_AMPK,kac_AKT,kac_PHLPP2,kac_PP2Ca,kdac_AMPK,kdac_AKT,kdac_PHLPP2,kdac_PP2Ca,l_PP2Ca,l_PHLPP2,l_AKT,l_AMPK,n_PP2Ca,n_PHLPP2,n_AKT,n_AMPK,PP2Ca_0_AMPK,PHLPP2_0_AKT,AMPK_0_PHLPP2,AKT_0_PP2Ca],ap);
[x,v,s,h,f] = cont(@equilibrium, x0, v0,opt);

    end