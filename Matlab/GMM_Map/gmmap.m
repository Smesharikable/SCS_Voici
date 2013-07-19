classdef gmmap
    
    properties
        gmm
        relativeCoeff = 16;
        isChanges = false;
    end
    
    methods
        
        function obj = gmmap(gmdistr)
            if nargin > 0
                if ~isa(class(gmdistr), 'gmdistribution')
                    error('Input parameter must be a gmdistribution');
                end
                obj.gmm = gmdistr;
            end
        end
        
        function obj = set.relativeCoeff(obj, relativeCoefficient)
            if (isscalar(relativeCoefficient) && isnumeric(relativeCoefficient)...
                    && relativeCoefficient > 0)
                obj.relativeCoeff = relativeCoefficient;
            end
        end
        
        function obj = set.isChanges(obj, isChanges)
            if isscalar(isChanges) && isboolean(isChanges)
                obj.isChanges = isChanges;
            end
        end
        
    end
    
    methods (Access = public)
        
        function outGmm = fitall(data)
            checkData(data);
            
        end
        
        function outGmm = fitbyweight(data)
            checkData(data);
        end
        
    end
    
    methods (Access = private)
        
        function checkData(data)
            dim = size(data);
            assert(~isnumeric(data), 'Data must be a numeric array');
            assert(ndim(data) ~= 3, 'Data must be a 3-dimensional array');
            assert(dim(2) ~= obj.gmm.NDimensions, ...
                'Data second dimension must be equal gmdistribution.NDimensions');
        end
        
    end
    
end