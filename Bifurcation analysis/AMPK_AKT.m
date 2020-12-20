function out = AMPK_AKT
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];
%-----------------------------------------------------------------------------------------------------
function dydt = fun_eval(t,SV,t_AMPK,t_AKT,t_PHLPP2,t_PP2Ca,kac_AMPK,kac_AKT,kac_PHLPP2,kac_PP2Ca,kdac_AMPK,kdac_AKT,kdac_PHLPP2,kdac_PP2Ca,l_PP2Ca,l_PHLPP2,l_AKT,l_AMPK,n_PP2Ca,n_PHLPP2,n_AKT,n_AMPK,PP2Ca_0_AMPK,PHLPP2_0_AKT,AMPK_0_PHLPP2,AKT_0_PP2Ca)
%hill functions
H_PHLPP2_on_AKT = (1+l_PHLPP2*(SV(3)/PHLPP2_0_AKT)^n_PHLPP2)/(1+(SV(3)/PHLPP2_0_AKT)^n_PHLPP2);
H_PP2Ca_on_AMPK = (1+l_PP2Ca*(SV(4)/PP2Ca_0_AMPK)^n_PP2Ca)/(1+(SV(4)/PP2Ca_0_AMPK)^n_PP2Ca);
H_AKT_on_PP2Ca = (1+l_AKT*(SV(2)/AKT_0_PP2Ca)^n_AKT)/(1+(SV(2)/AKT_0_PP2Ca)^n_AKT);
H_AMPK_on_PHLPP2 = (1+l_AMPK*(SV(1)/AMPK_0_PHLPP2)^n_AMPK)/(1+(SV(1)/AMPK_0_PHLPP2)^n_AMPK);
%H_Stress_on_AMPK = (1+l_PP2Ca*(SV(4)/PP2Ca_0_AMPK)^n_PP2Ca)/(1+(SV(4)/PP2Ca_0_AMPK)^n_PP2Ca);

%equations

dydt = [
        kac_AMPK*(t_AMPK - SV(1)) - kdac_AMPK*SV(1)*H_PP2Ca_on_AMPK; % AMPK 
        kac_AKT*(t_AKT - SV(2)) - kdac_AKT*SV(2)*H_PHLPP2_on_AKT; % AKT
        kac_PHLPP2*(t_PHLPP2 - SV(3))*H_AMPK_on_PHLPP2 - kdac_PHLPP2*SV(3); % PHLPP2
        kac_PP2Ca*(t_PP2Ca - SV(4))*H_AKT_on_PP2Ca - kdac_PP2Ca*SV(4); % PP2Ca
       ];
    
    % --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(AMPK_AKT);
y0=[0,0,0,0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% % --------------------------------------------------------------------------
% function jac = jacobian(t,kmrgd,S,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,S10u200,S20mz,u2000,nzu200,nS1u200,nzmz,nS2mz,nu200,lamdazu200,lamdaS1u200,lamdazmz,lamdaS2mz)
% % --------------------------------------------------------------------------
% function jacp = jacobianp(t,kmrgd,S,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,S10u200,S20mz,u2000,nzu200,nS1u200,nzmz,nS2mz,nu200,lamdazu200,lamdaS1u200,lamdazmz,lamdaS2mz)
% %--------------------------------------------------------------------------
% function hess = hessians(t,kmrgd,S,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,S10u200,S20mz,u2000,nzu200,nS1u200,nzmz,nS2mz,nu200,lamdazu200,lamdaS1u200,lamdazmz,lamdaS2mz)
% %--------------------------------------------------------------------------
% function hessp = hessiansp(t,kmrgd,S,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,S10u200,S20mz,u2000,nzu200,nS1u200,nzmz,nS2mz,nu200,lamdazu200,lamdaS1u200,lamdazmz,lamdaS2mz)
% %---------------------------------------------------------------------------
% function tens3  = der3(t,kmrgd,S,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,S10u200,S20mz,u2000,nzu200,nS1u200,nzmz,nS2mz,nu200,lamdazu200,lamdaS1u200,lamdazmz,lamdaS2mz)
% %---------------------------------------------------------------------------
% function tens4  = der4(t,kmrgd,S,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,S10u200,S20mz,u2000,nzu200,nS1u200,nzmz,nS2mz,nu200,lamdazu200,lamdaS1u200,lamdazmz,lamdaS2mz)
% %---------------------------------------------------------------------------
% function tens5  = der5(t,kmrgd,S,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,S10u200,S20mz,u2000,nzu200,nS1u200,nzmz,nS2mz,nu200,lamdazu200,lamdaS1u200,lamdazmz,lamdaS2mz)
