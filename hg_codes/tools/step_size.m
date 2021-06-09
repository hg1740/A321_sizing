function delta=step_size(x)

coefficient = ones(1,length(x));

for i=1:numel(x)
    
    if x(i)>2
        
        coefficient(i)=1.3;
        
    elseif x(i)<=2 && x(i)>1
        
        coefficient(i)=0.3*(x(i)-1)+1;
        
    elseif x(i)<1 && x(i)>=0.5
        
        coefficient(i)=1-0.6*(1-x(i));
        
    elseif x(i)<0.5
        
        coefficient(i)=0.7;
        
    end
    
    delta=coefficient;
      
end


