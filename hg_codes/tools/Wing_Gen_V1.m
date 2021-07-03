

function [Wing, varargout]= Wing_Gen_V1(Param)

    
%% Wing configurations for starboard wing
  
    Aspect_ratio=Param.Wing.AR; % Aspect ratio = 10.172 for A321 model
    
    Total_area=Param.Wing.TotalArea;         % include two wing surface areas + floor size on the fuselage
    Fuselage_width=Param.Layout.Fuselage_Width;       % dimeter of the fuselage
  
    Wing_span = sqrt(Aspect_ratio*Total_area);
    BeamLoc = 0.4;          % choose a spar location: 0 --> 1
    Semi_span=(Wing_span-Fuselage_width)/2; % length of one wing: 16m for A321 model
    
    Root_chord =  Total_area/(1.064*Semi_span + 4);
    LE_sweep=27;            % deg
    
    Wing_area = (Total_area - Fuselage_width*Root_chord)/2;
    
    Mid_chord=0.63685*Root_chord;
    Tip_chord=0.2248*Root_chord;
    
    X0=Root_chord; 
    X1=0.27*Semi_span*tan(27*pi/180) + Mid_chord;
    X2=Semi_span*tan(27*pi/180) + Tip_chord;
    
    tan_TE_sweep1=(X1-X0)/(0.27*Semi_span);
    tan_TE_sweep2=(X2-X1)/(0.73*Semi_span);
    
    TE_sweep1=atan(tan_TE_sweep1)*180/pi; % deg
    TE_sweep2=atan(tan_TE_sweep2)*180/pi; % deg
    
    %% output
    
    Wing.Root_Chord=Root_chord;
    Wing.LE_Sweep=27;
    Wing.TE_Sweep1=TE_sweep1;
    Wing.TE_Sweep2=TE_sweep2;
    Wing.Span=Wing_span;
    Wing.Semi_Span=Semi_span;
    
    %% obtain wingbox geometric properties
    
    if isfield(Param,'FWT') 
        
        fold_eta=Param.FWT.Fold_eta;

        FWT__root_tc=interp1([0,Param.Wing.Kink,1],[Param.Wing.ThicknessToChord_Root,Param.Wing.ThicknessToChord_kink,Param.Wing.ThicknessToChord_tip],fold_eta);
        
        FWT__root_chord=interp1([0,Param.Wing.Kink,1 ],[Root_chord,Mid_chord,Tip_chord],fold_eta);
        
        FWT__root_height=FWT__root_chord*FWT__root_tc;
        
        FWT__tip_tc=Param.Wing.ThicknessToChord_tip;
        
        FWT__tip_chord=Tip_chord;
        
        FWT__tip_height=Tip_chord*FWT__tip_tc;

        %% outputs   
        
        FWT_param.Root_Chord=FWT__root_chord;
        FWT_param.Root_Height=FWT__root_height;
        
        FWT_param.Tip_Chord=FWT__tip_chord;
        FWT_param.Tip_Height=FWT__tip_height;
        
        varargout{1}=FWT_param;
        
    end
    

  
    
    
end