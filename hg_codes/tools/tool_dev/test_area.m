function Area = test_area(constant,varargin)



p=inputParser;

addParameter(p,'width',1,@(x)validateattributes(x,{'numeric'},{'scalar', 'nonnan', 'finite', 'real','nonempty'}))
addParameter(p,'length',1,@(x)validateattributes(x,{'numeric'},{'scalar', 'nonnan', 'finite', 'real','nonempty'}))
addParameter(p,'disp','s',@(x)validateattributes(x,{'char'},{'nonempty'}))

parse(p,varargin{:})

Area=constant + p.Results.width * p.Results.length;

disp(p.Results.disp)


end