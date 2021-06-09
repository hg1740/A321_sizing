
% 1-25 thickness 1 spar
% 26-50 thickness 2 skin
% 51-75 Astrg
% initial values
x=[ones(1,25)*0.003,ones(1,25)*0.003,ones(1,25)*0.00001];

%initial condition 
cond_set1=[1.1,1.1,1.1,1.1];
cond_set2=[1.1,1.1,1.1,1.1];

record=zeros(100,75);
counter=1;

while max(cond_set1)>1 || min(cond_set2)<0.95
    
    
    record(counter,:)=x(1,:);
    
    [RVon_skn, RVon_spr, Rsigmab_skn, Rtaub_skn, Rsigma_strg, Rsigma_crip, Rsigma_col]=BeamStress_calc(x);
    
    % skin check 
    Check_von_skin1=RVon_skn/5e8;
    Check_von_skin2=RVon_skn./Rsigmab_skn;
    
    [~,index1]=find(Check_von_skin1>1);
    [~,index2]=find(Check_von_skin1<0.95);
    [~,index3]=find(Check_von_skin2>1);
    [~,index4]=find(Check_von_skin2<0.95);
    
    skin_inc=unique([index1,index3]); % one of the condition suggest increase then need to increase t
    
    skin_dec=intersect(index2,index4); % only both condition suggest to decrease then decrease
    
    % find the max. coeff. to increase 
    increase_coeff1_=Check_von_skin1(skin_inc);
    increase_coeff2_=Check_von_skin2(skin_inc);
    increase_coeff=max([increase_coeff1_; increase_coeff2_]);
     
    % find max. to decrease, decrease slowly 
    
    decrease_coeff1_=Check_von_skin1(skin_dec);
    decrease_coeff2_=Check_von_skin2(skin_dec);
    decrease_coeff=max([decrease_coeff1_; decrease_coeff2_]);
    
    x(1,25+skin_inc)=x(1,25+skin_inc)*1.1;
    x(1,25+skin_dec)=x(1,25+skin_dec)*0.9;
    
    % spar check 
    Check_von_spar=RVon_spr/5e8;
   
    [~,sp_index1]=find(Check_von_spar>1);
    [~,sp_index2]=find(Check_von_spar<0.95);
       
    spar_inc=sp_index1;
    spar_dec=sp_index2;
    
    sp_inc_co=Check_von_spar(spar_inc);
    sp_dec_co=Check_von_spar(spar_dec);
    
    x(1,spar_inc)=x(1,spar_inc)*1.1;
    x(1,spar_dec)=x(1,spar_dec)*0.9;
    
    % stringers check 
    Check_strg=RVon_spr./Rsigma_crip;
    
    [~,sg_index1]=find(Check_strg>1);
    [~,sg_index2]=find(Check_strg<0.95);
    
    strg_inc=sg_index1;
    strg_dec=sg_index2;
    
    strg_inc_co=Check_strg(strg_inc);
    strg_dec_co=Check_strg(strg_dec);
    
    x(1,50+strg_inc)=x(1,50+strg_inc)*1.1;
    x(1,50+strg_dec)=x(1,50+strg_dec)*0.9;
 
    % end condition check 
    cond1=max(Check_von_skin1);
    cond2=max(Check_von_skin2);
    cond3=max(Check_von_spar);
    cond4=max(Check_strg);
    
    cond5=min(Check_von_skin1);
    cond6=min(Check_von_skin2);
    cond7=min(Check_von_spar);
    cond8=min(Check_strg);
    
    cond_set1=[cond1,cond2,cond3,cond4];
    cond_set2=[cond6,cond7,cond8];
    
    disp(cond_set1);
    disp(cond_set2);
    
    counter=counter+1;
    
    
end



% %tensile
% cond1=max(RVon_skn/5e8); 
% cond2=max(RVon_spr/5e8);
% 
% %compressive
% cond3=max(RVon_skn./Rsigmab_skn);
% cond4=max(Rsigma_strg./Rsigma_crip);