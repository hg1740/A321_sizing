function output=test(varargin)

    
    p=inputParser;

    
    addParameter(p,'Load_Factor',1,@(x)validateattributes(x,{'numeric'},{'scalar', 'nonnan', 'finite', 'real','nonempty'}))
    addParameter(p,'Name','file_name',@(x)validateattributes(x,{'char'},{'nonempty'}))
    addParameter(p,'Hinge_Lock','on',@(x)any(validatestring(x,{'on','off'})))
    
    parse(p,varargin{:})
    
    disp(p.Results.Load_Factor)
    disp(p.Results.Name)
    disp(p.Results.Hinge_Lock)
    
    output=p.Results;
    
    
end