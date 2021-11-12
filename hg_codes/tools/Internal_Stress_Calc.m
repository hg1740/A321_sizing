function Stresses=Internal_Stress_Calc(Param, Box_CrossSec, Box_dimensions, Loads)

%% Material properties

    E=70e9;
    sigma_Y=5.2E8;
    
%% Wing-box cross sectional parameters 

if isfield(Param,'FWT')
    
    Spar_Cap_Thickness=Param.Wing.SparCap_Thickness;
    Spar_Web_Thickness=Param.Wing.SparWeb_Thickness;
    Skin_Thickness=Param.Wing.Skin_Thickness;
    Stringer_Area=Param.Wing.Stringer_Area;
    
    FWT_Spar_Cap=Param.FWT.SparCap_Thickness;
    FWT_Spar_Web=Param.FWT.SparWeb_Thickness;
    FWT_Skin=Param.FWT.Skin_Thickness;
    FWT_Stringer=Param.FWT.Stringer_Area;
    
    % Combine
    t_sc=[Spar_Cap_Thickness, FWT_Spar_Cap(2:end)];
    t_sw=[Spar_Web_Thickness, FWT_Spar_Web(2:end)];
    t_skin=[Skin_Thickness, FWT_Skin(2:end)];
    A_strg=[Stringer_Area,FWT_Stringer(2:end)];
    
    Stringer_length=sqrt(A_strg/0.36);
    Stringer_thickness=0.12*Stringer_length;
    
    % short hand
    ds=Stringer_length;
    ts=Stringer_thickness;
    
    NumSec=35;
    
else
    
    Spar_Cap_Thickness=Param.Wing.SparCap_Thickness;
    Spar_Web_Thickness=Param.Wing.SparWeb_Thickness;
    Skin_Thickness=Param.Wing.Skin_Thickness;
    Stringer_Area=Param.Wing.Stringer_Area;
    
    Stringer_length=sqrt(Stringer_Area/0.36);
    Stringer_thickness=0.12*Stringer_length;
    
    % short hand
    t_sc=Spar_Cap_Thickness;
    t_sw=Spar_Web_Thickness;
    t_skin=Skin_Thickness;
    A_strg=Stringer_Area;
    ds=Stringer_length;
    ts=Stringer_thickness;
    
    NumSec=25;
    
end


    
    %% Internal loads 
    
    % initialise stress values
    
    % skin
    sigma_skin=ones(1,NumSec);
    tau_skin=ones(1,NumSec);
    sigma_pp=ones(1,NumSec);
    tau_pp=ones(1,NumSec);
    Von_skin=ones(1,NumSec);
    
    skin_buck_sigma=ones(1,NumSec);
    skin_buck_tau=ones(1,NumSec);

    % spar cap
    sigma_spr1=ones(1,NumSec); 
    Von_spr1=ones(1,NumSec);
    sparcap_buck_sigma=ones(1,NumSec);
    sparcap_sheetbuck_sigma=ones(1,NumSec);
    
    % spar web
    tau_spr2=ones(1,NumSec);
    Von_spr2=ones(1,NumSec);
    spr2_buck_tau=ones(1,NumSec);
    
    N_crush=ones(1,NumSec);

    % stringer
    sigma_strg=ones(1,NumSec);
    sigma_crip=ones(1,NumSec);
    sigma_col=ones(1,NumSec);
    
    
    for jj=1:NumSec
        
        % External loads
        M_P2=Loads.Moment_P2(jj);
        T=Loads.Torque(jj);
        S_P2=Loads.Shear_P2(jj);
        
        % box dimensions        
        hs=Box_dimensions.Inboard.Height(jj);
        ws=Box_dimensions.Inboard.Width(jj);
        
        lcw=Param.Wing.CapEta_width*ws; % spar cap top length
        lch=Param.Wing.CapEta_height*hs; % spar cap web length
        
        % bending stiffness
        Iyy=Box_CrossSec.Iyy(jj);
        
        % stringer pitch
        strg_n=Param.Wing.StringerPitch;
        rib_n=0.60;

        % skin
        
        % normal stress
        sigma_skin(jj)=0.5*M_P2*hs/Iyy;
        
        % shear stress
        tau_skin(jj)= T/(2*hs*ws*t_skin(jj));
        
        % skin VonMise
        
        Von_skin(jj)=sqrt(sigma_skin(jj)^2 + 3*tau_skin(jj)^2);
        
        % princial stresses
        sigma_pp(jj) = 0.5*sigma_skin(jj) + sqrt((0.5*sigma_skin(jj))^2 + tau_skin(jj)^2);
        tau_pp(jj) = sqrt((0.5*sigma_skin(jj))^2+tau_skin(jj)^2);
        
        % critical buckling stress - for skin--------------------------
        skin_buck_sigma(jj)=4.0*pi^2*E*(t_skin(jj)/strg_n)^2/(12*(1-0.33^2));
        skin_buck_tau(jj)=5.6*pi^2*E*(t_skin(jj)/strg_n)^2/(12*(1-0.33^2));
        
        
        % spar cap: VonMise
        sigma_spr1(jj)=0.5*M_P2*hs/Iyy;
        Von_spr1(jj)=sigma_spr1(jj);
         
        % spar cap: column buckling
        sparcap_Iyy = t_sc(jj)*lch^3/12 + lch^3*t_sc(jj)/4 + (lch + lcw)*t_sc(jj)*(lch^2/(2*lch + 2*lcw))^2;
        
        sparcap_area=(lch+lcw)*t_sc(jj);
        
        R=sqrt(sparcap_Iyy/sparcap_area);
        
        sparcap_buck_sigma(jj)=pi^2*E/(rib_n/R);
        
        % spar cap: sheet buckling
        sparcap_sheetbuck_sigma(jj)=1.2*E*(t_sc(jj)/strg_n)^2;
        

        % spar web: take vertical shear force
        Q=ws*t_skin(jj)*hs/2 + hs^2*t_sw(jj)/4 + (lcw*t_sc(jj)*(0.5*hs) + lch*t_sc(jj)*(0.5*hs - 0.5*lch))*2;
        
        tau_spr2(jj)=S_P2*Q/(2*Iyy*t_sw(jj)) + T/(2*hs*ws*t_sw(jj));
        
        Von_spr2(jj)=sqrt(3*tau_spr2(jj)^2);
        
        
        % spar web: crushing force 
        N_crush(jj)=0.5*M_P2^2*hs*t_sc(jj)/(Iyy^2*E);
        

        % spar web: critical shear buckling stress
        ar=rib_n/hs;
        
        Ks=(0.121*ar^2 + 3.71*ar +2)/(ar-0.29);
        
        spr2_buck_tau(jj)=Ks*E*(t_sw(jj)/hs);
        
        
        % stringer
        sigma_strg(jj)=M_P2*(0.5*hs-0.5*ds(jj))/Iyy;
        
        
        % striner crippling stress--------------------------------------
        A1=-0.7885; B1=0.6194; A2=-0.8046; B2=1.2117;
        
        Rc1_=B1*(sqrt(sigma_Y/E)*(ds(jj)/ts(jj)))^A1;
        Rc2_=B2*(sqrt(sigma_Y/E)*(ds(jj)/ts(jj)))^A2;
        
        sigmac1_=min(Rc1_,1.45)*sigma_Y;
        sigmac2_=min(Rc2_,1.45)*sigma_Y;
        
        sigma_crip(jj)=(ds(jj)*ts(jj)*sigmac1_*2 + ds(jj)*ts(jj)*sigmac2_)/(ts(jj)*ds(jj)*3);
        %-------------------------------------------------------------
        
        % column buckling

        A=3*ds(jj)*ts(jj);
        
        I=(ts(jj)*ds(jj)^3/12)+(ds(jj)*ts(jj)^3/12 + ds(jj)*ts(jj)*(ds(jj)/2)^2)*2;
        
        rg=sqrt(I/A);

        K=rib_n/rg;
        
        Kcr=pi*sqrt(2*E/sigma_crip(jj));
        
        if K < Kcr
            
            sigma_col(jj)=sigma_crip(jj)*(1-sigma_crip(jj)*K^2/(4*pi^2*E));
            
        elseif K > Kcr
            
            sigma_col(jj)=pi^2*E/K^2;
            
        end
          
        
    end
    
    %% output results
    
    Stresses.skin_normal=sigma_skin;
    Stresses.skin_shear=tau_skin;
    Stresses.skin_principal_normal=sigma_pp;
    Stresses.skin_principal_shear=tau_pp;
    Stresses.skin_VonMise=Von_skin;
    
    Stresses.skin_critical_buck_normal=skin_buck_sigma;
    Stresses.skin_critical_buck_shear=skin_buck_tau;
    
    Stresses.sparcap_VonMise=Von_spr1;
    Stresses.sparcap_ColumnBuckling=sparcap_buck_sigma;
    Stresses.sparcap_SheetBuckling=sparcap_sheetbuck_sigma;
    
    Stresses.sparweb_VonMise=Von_spr2;
    Stresses.sparweb_Shear=tau_spr2;
    Stresses.sparweb_ShearBuckling=spr2_buck_tau;
    Stresses.sparweb_Crush=N_crush;
    
    Stresses.Strg_compression=sigma_strg;
    Stresses.Strg_crippling=sigma_crip;
    Stresses.Strg_column_buckle=sigma_col;
    
    Stresses.Y=Loads.Y;
    


end