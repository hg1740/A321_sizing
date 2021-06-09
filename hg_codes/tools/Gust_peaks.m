function [Max_Moment,Min_Moment, Max_Torque, Min_Torque, Max_Shear, Min_Shear]=Gust_peaks(WingNodes,LoadCase,run_folder,file_name,NumSteps)

    
    % data extraction 
    Gust_force=h5read(strcat(run_folder,file_name),'/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BEAM');
    
    % Find number of total nodes in the model 
    [~,ind]=find(Gust_force.GRID(1,:)==Gust_force.GRID(1));
    NumGrid=ind(2)-ind(1)+1;
    
    % Find the number of gust in the load case
    NumGust=numel(LoadCase.GustLength);
    
    % index for all the wing nodes
    All_index=ismember([Gust_force.GRID(1,:)],[WingNodes.GID]);
    
    % loads at all wing nodes
    All_Moment2=Gust_force.BM2(1,All_index);
    All_Torque=Gust_force.TTRQ(1,All_index);
    All_Shear2=Gust_force.TS2(1,All_index);
    
    Moment2_matrix=reshape(All_Moment2,numel(WingNodes)-1,NumGust*NumSteps);
    Torque_matrix=reshape(All_Torque,numel(WingNodes)-1,NumGust*NumSteps);
    Shear_matrix=reshape(All_Shear2,numel(WingNodes)-1,NumGust*NumSteps);

    % delta values for bending moments, torque and shear forces
    Max_Moment = max( Moment2_matrix, [], 2);
    Min_Moment = min( Moment2_matrix, [], 2);

    Max_Torque = max(Torque_matrix, [], 2);
    Min_Torque = min(Torque_matrix, [], 2);

    Max_Shear = max(Shear_matrix, [], 2);
    Min_Shear = min(Shear_matrix, [], 2);
    

end