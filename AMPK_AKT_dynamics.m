function dydt = AMPK_AKT_dynamics(t,SV,t_AMPK,t_AKT,t_PHLPP2,t_PP2Ca,kac_AMPK,kac_AKT,kac_PHLPP2,kac_PP2Ca,kdac_AMPK,kdac_AKT,kdac_PHLPP2,kdac_PP2Ca,l_PP2Ca,l_PHLPP2,l_AKT,l_AMPK,n_PP2Ca,n_PHLPP2,n_AKT,n_AMPK,PP2Ca_0_AMPK,PHLPP2_0_AKT,AMPK_0_PHLPP2,AKT_0_PP2Ca)
%hill functions
H_PHLPP2_on_AKT = (1+l_PHLPP2*(SV(3)/PHLPP2_0_AKT)^n_PHLPP2)/(1+(SV(3)/PHLPP2_0_AKT)^n_PHLPP2);
H_PP2Ca_on_AMPK = (1+l_PP2Ca*(SV(4)/PP2Ca_0_AMPK)^n_PP2Ca)/(1+(SV(4)/PP2Ca_0_AMPK)^n_PP2Ca);
H_AKT_on_PP2Ca = (1+l_AKT*(SV(2)/AKT_0_PP2Ca)^n_AKT)/(1+(SV(2)/AKT_0_PP2Ca)^n_AKT);
H_AMPK_on_PHLPP2 = (1+l_AMPK*(SV(1)/AMPK_0_PHLPP2)^n_AMPK)/(1+(SV(1)/AMPK_0_PHLPP2)^n_AMPK);

%equations
dydt = [
        kac_AMPK*(t_AMPK - SV(1)) - kdac_AMPK*SV(1)*H_PP2Ca_on_AMPK; % AMPK 
        kac_AKT*(t_AKT - SV(2)) - kdac_AKT*SV(2)*H_PHLPP2_on_AKT; % AKT
        kac_PHLPP2*(t_PHLPP2 - SV(3))*H_AMPK_on_PHLPP2 - kdac_PHLPP2*SV(3); % PHLPP2
        kac_PP2Ca*(t_PP2Ca - SV(4))*H_AKT_on_PP2Ca - kdac_PP2Ca*SV(4); % PP2Ca
       ];
end