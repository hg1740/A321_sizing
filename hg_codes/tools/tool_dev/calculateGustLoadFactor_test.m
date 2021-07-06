function varargout = calculateGustLoadFactor_test(obj, Aircraft, CL_alfa) 
            %calculateGustLoadFactor Calculates the vertical load factor
            %using the equation for a Pratt Gust.                        
            
            validateattributes(Aircraft, {'awi.model.Aircraft'}, ...
                {'scalar'}, 'calculateGustLoadFactor', 'Aircraft');
            if nargin < 3
                %Assume the aircraft lift-curve slope is equal to 2pi
                CL_alfa = 2 * pi;
            else
                validateattributes(CL_alfa, {'numeric'}, {'row', 'positive'}, 'calculateGustLoadFactor', 'CL_alfa');
            end
            
            %Only proceed using LoadCases of type 'Pratt Gust'
            obj = obj(ismember({obj.LoadCaseType}, 'Pratt Gust'));
            if isempty(obj) %Escape route
                return
            end
            
            %Check we have the correct inputs
            props = get(Aircraft, {'RefArea', 'RefSpan'});
            if any(cellfun(@isempty, props)) %Escape route 
                return
            end
            area = props{1};
            span = props{2};
            chrd = area / span;
            props = get(obj, {'Altitude', 'AcMass', 'Mach'});
            if any(any(cellfun(@isempty, props))) || ...
                    any(any(cellfun(@isnan, props))) %Escape route 
                return
            end
            FP   = getFlightPointData(obj);
            rho  = [FP.Density];
            FPv  = [FP.AcVelocity];
            alt  = horzcat(props{:, 1});
            mass = horzcat(props{:, 2});
            
            %Sometimes the aircraft velocity is explicitly defined. If that
            %is the case then use this data instead of the flight point
            %data.
            vel      = [obj.AcVelocity];
            idx      = isnan(vel);
            vel(idx) = FPv(idx);
            
            %Get mass ratio 'muG' and gust load alleviation factor 'kG'
            muG = mass ./ (0.5 .* rho .* area .* chrd .* CL_alfa);
            kG  = (0.88 .* muG) ./ (5.3 + muG);
            
            %Calculate gust velocity
            Uref = [ ...
                0    , 15000, 60000 ; ...
                17.07, 13.41, 6.36 ];
            Ug = interp1(Uref(1, :), Uref(2, :), alt);  %[m/s], EAS
            Ug = Ug ./ sqrt(rho ./ obj(1).RefDensity);  %[m/s], TAS
            
            %Calculate incremental load factor
            nZ = 1 + (kG .* CL_alfa .* rho .* area .* vel .* Ug ./ (2 .* mass * 9.81));           
            %VT = kG .*(1 + Ug ./ vel) .* (1 + Ug ./ vel) .^2;
            
            %User asked for it back?
            if nargout > 0
                %origin
                varargout{1} = nZ;
                %charles:add
                varargout{2} = Ug;

            else
                %Assign to the object
                set(obj, {'LoadFactor'}, num2cell(nZ'));
            end            
            
        end