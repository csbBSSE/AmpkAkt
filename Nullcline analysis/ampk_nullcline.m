function dydt = ampk_nullcline(t,SV,AMPKc,t_AMPK,t_AKT,t_PHLPP2,t_PP2Ca,kac_AMPK,kac_AKT,kac_PHLPP2,kac_PP2Ca,kdac_AMPK,kdac_AKT,kdac_PHLPP2,kdac_PP2Ca,l_PP2Ca,l_PHLPP2,l_AKT,l_AMPK,n_PP2Ca,n_PHLPP2,n_AKT,n_AMPK,PP2Ca_0_AMPK,PHLPP2_0_AKT,AMPK_0_PHLPP2,AKT_0_PP2Ca)
%hill functions
H_PHLPP2_on_AKT = (1+l_PHLPP2*(SV(2)/PHLPP2_0_AKT)^n_PHLPP2)/(1+(SV(2)/PHLPP2_0_AKT)^n_PHLPP2);
H_PP2Ca_on_AMPK = (1+l_PP2Ca*(SV(3)/PP2Ca_0_AMPK)^n_PP2Ca)/(1+(SV(3)/PP2Ca_0_AMPK)^n_PP2Ca);
H_AKT_on_PP2Ca = (1+l_AKT*(SV(1)/AKT_0_PP2Ca)^n_AKT)/(1+(SV(1)/AKT_0_PP2Ca)^n_AKT);
H_AMPK_on_PHLPP2 = (1+l_AMPK*(AMPKc/AMPK_0_PHLPP2)^n_AMPK)/(1+(AMPKc/AMPK_0_PHLPP2)^n_AMPK);
%H_Stress_on_AMPK = (1+l_PP2Ca*(SV(3)/PP2Ca_0_AMPK)^n_PP2Ca)/(1+(SV(3)/PP2Ca_0_AMPK)^n_PP2Ca);

%equations

dydt = [
        %kac_AMPK*(t_AMPK - AMPKc) - kdac_AMPK*AMPKc*H_PP2Ca_on_AMPK; % AMPK 
        kac_AKT*(t_AKT - SV(1)) - kdac_AKT*SV(1)*H_PHLPP2_on_AKT; % AKT
        kac_PHLPP2*(t_PHLPP2 - SV(2))*H_AMPK_on_PHLPP2 - kdac_PHLPP2*SV(2); % PHLPP2
        kac_PP2Ca*(t_PP2Ca - SV(3))*H_AKT_on_PP2Ca - kdac_PP2Ca*SV(3); % PP2Ca
       ];
end