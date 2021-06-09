function [Moment, Torque, Shear]=Gust_Thistory(WingNodes,LoadCase,run_folder,file_name,NumSteps)

    
    % data extraction 
    Gust_force=h5read(strcat(run_folder,file_name),'/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BEAM');
    
    % Find number of total nodes in the model 
    [~,ind]=find(Gust_force.GRID(1,:)==Gust_force.GRID(1));
    NumGrid=ind(2)-ind(1)+1;
    
    % Find the number of gust in the load case
    NumGust=numel(LoadCase.GustLength);
    
    % index for the node at wing root
    All_index=ismember([Gust_force.GRID(1,:)],[WingNodes(1).GID]);
    
    % loads at all wing nodes
    All_Moment2=Gust_force.BM2(1,All_index);
    All_Torque=Gust_force.TTRQ(1,All_index);
    All_Shear2=Gust_force.TS2(1,All_index);
    
    % Time history of the forces

    Moment=reshape(All_Moment2,NumSteps, NumGust);
    Torque=reshape(All_Torque, NumSteps, NumGust);
    Shear=reshape(All_Shear2,NumSteps, NumGust);
    

   

end