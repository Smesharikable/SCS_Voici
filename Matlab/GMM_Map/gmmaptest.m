classdef gmmaptest
    
    properties (Access = public)
        relativeCoeff = 16;
        visualizing = false;
    end
    
    properties (GetAccess = public, SetAccess = protected)
        gmmin
        gmmout
        pstat = [];
        Exstat = [];
        newp = [];
        newmu = [];
        coef =[];
    end
    
    methods
        
        function obj = gmmaptest(gmdistr)
            if nargin > 0
                if ~isa(gmdistr, 'gmdistribution')
                    error('Input parameter must be a gmdistribution');
                end
                obj.gmmin = gmdistr;
                obj.pstat = zeros(1, gmdistr.NComponents);
                obj.Exstat = zeros(gmdistr.NComponents, gmdistr.NDimensions);
            end
        end
        
        function obj = set.gmmout(obj, gmm)
            if ~isa(gmm, 'gmdistribution')
                error('Input parameter must be a gmdistribution');
            end
            obj.gmmout = gmm;
        end
        
        function obj = set.relativeCoeff(obj, relativeCoefficient)
            if (isscalar(relativeCoefficient) && isnumeric(relativeCoefficient)...
                    && relativeCoefficient > 0)
                obj.relativeCoeff = relativeCoefficient;
            end
        end
        
        function obj = set.visualizing(obj, visualize)
            if isscalar(visualize) && islogical(visualize)
                obj.visualizing = visualize;
            end
        end
        
    end
    
    methods (Access = public)
        
        % data - matrix n-by-d, where n is the number of observations,
        %   and d is the dimensions of the data
        % equal coefficients
        function test = fitbymeans(obj, data, iterations)
            obj.checkData(data);
            outGmm = obj.gmmin;
            T = length(data(:, 1));
            oldp = obj.gmmin.PComponents;
            oldmu = obj.gmmin.mu;
            for k = 1:iterations
                Pr = posterior(outGmm, data);
                obj.newp = zeros(1, obj.gmmin.NComponents);
                obj.newmu = zeros(obj.gmmin.NComponents, obj.gmmin.NDimensions);
                obj.coef = obj.pstat / (obj.pstat + obj.relativeCoeff);
                for i = 1:obj.gmmin.NComponents
                    obj.pstat(i) = sum(Pr(:, i));
                    obj.Exstat(i, :) = (Pr(:, i)' * data) / obj.pstat(i);
                    obj.newp(i) = obj.coef * obj.pstat(i) / T + (1 - obj.coef) * oldp(i);
                    obj.newmu(i, :) = obj.coef * obj.Exstat(i, :) + (1 - obj.coef) * oldmu(i, :);
                end
                outGmm = gmdistribution(obj.newmu, obj.gmmin.Sigma, oldp);
                oldmu = obj.newmu;
            end
            if obj.visualizing
                obj.visualize(data, outGmm);
            end
            obj.gmmout = outGmm;
            obj.gmmin = outGmm;
            test = obj;
        end
        
        % different coefficients
        function test = fitbymeans2(obj, data, iterations)
            obj.checkData(data);
            outGmm = obj.gmmin;
            T = length(data(:, 1));
            oldp = obj.gmmin.PComponents;
            oldmu = obj.gmmin.mu;
            for k = 1:iterations
                Pr = posterior(outGmm, data);
                obj.newp = zeros(1, obj.gmmin.NComponents);
                obj.newmu = zeros(obj.gmmin.NComponents, obj.gmmin.NDimensions);
                obj.coef = zeros(1, obj.gmmin.NComponents);
                for i = 1:obj.gmmin.NComponents
                    obj.pstat(i) = sum(Pr(:, i));
                    obj.coef(i) = obj.pstat(i) / (obj.pstat(i) + obj.relativeCoeff);
                    obj.Exstat(i, :) = (Pr(:, i)' * data) / obj.pstat(i);
                    obj.newp(i) = obj.coef(i) * obj.pstat(i) / T + (1 - obj.coef(i)) * oldp(i);
                    obj.newmu(i, :) = obj.coef(i) * obj.Exstat(i, :) + (1 - obj.coef(i)) * oldmu(i, :);
                end
                outGmm = gmdistribution(obj.newmu, obj.gmmin.Sigma, oldp);
                oldmu = obj.newmu;
            end
            if obj.visualizing
                obj.visualize(data, outGmm);
            end
            obj.gmmout = outGmm;
            obj.gmmin = outGmm;
            test = obj;
        end
        
        % fitting given gmm by weight using equal coefficients
        function test = fitbyweight(obj, data, iterations)
            obj.checkData(data);
            outGmm = obj.gmmin;
            T = length(data(:, 1));
            oldp = obj.gmmin.PComponents;
            for k = 1:iterations
                Pr = posterior(outGmm, data);
                obj.newp = zeros(1, obj.gmmin.NComponents);
                obj.coef = zeros(1, obj.gmmin.NComponents);
                obj.coef(1, :) = obj.pstat / (obj.pstat + obj.relativeCoeff);
                for i = 1:obj.gmmin.NComponents
                    obj.pstat(i) = sum(Pr(:, i));
                    obj.newp(i) = obj.coef(i) * obj.pstat(i) / T + (1 - obj.coef(i)) * oldp(i);
                end
                outGmm = gmdistribution(obj.gmmin.mu, obj.gmmin.Sigma, obj.newp);
            end
            if obj.visualizing
                obj.visualize(data, outGmm);
            end
            obj.gmmout = outGmm;
            test = obj;
        end
        
        % fitting given gmm by weight using different coefficients
        function test = fitbyweight2(obj, data, iterations)
            obj.checkData(data);
            outGmm = obj.gmmin;
            T = length(data(:, 1));
            oldp = obj.gmmin.PComponents;
            for k = 1:iterations
                Pr = posterior(outGmm, data);
                obj.newp = zeros(1, obj.gmmin.NComponents);
                obj.coef = zeros(1, obj.gmmin.NComponents);
                for i = 1:obj.gmmin.NComponents
                    obj.pstat(i) = sum(Pr(:, i));
                    obj.coef(i) = obj.pstat(i) / (obj.pstat(i) + obj.relativeCoeff);
                    obj.newp(i) = obj.coef(i) * obj.pstat(i) / T + (1 - obj.coef(i)) * oldp(i);
                end
                outGmm = gmdistribution(obj.gmmin.mu, obj.gmmin.Sigma, obj.newp);
            end
            if obj.visualizing
                obj.visualize(data, outGmm);
            end
            obj.gmmout = outGmm;
            test = obj;
        end
        
    end
    
    methods (Access = private)
        
        function checkData(obj, data)
            dim = size(data);
            assert(isnumeric(data), 'Data must be a numeric array');
            assert(ndims(data) == 2, 'Data must be a 2-dimensional array');
            assert(dim(2) == obj.gmmin.NDimensions, ...
                'Data second dimension must be equal gmdistribution.NDimensions');
        end
        
        function visualize(obj, data, outGmm)
            figure;
            hold on;
            scatter(data(:, 1), data(:, 2))
            oldmu = obj.gmmin.mu;
            for i = 1:obj.gmmin.NComponents
                plot([oldmu(i, 1), obj.Exstat(i, 1)], [oldmu(i, 2), obj.Exstat(i, 2)]);
            end
            scatter(oldmu(:, 1), oldmu(:, 2))
            scatter(obj.Exstat(:, 1), obj.Exstat(:, 2))
            ezcontour(@(x, y)pdf(outGmm, [x y]), [-15 15], [-15 15], 300)
        end
        
    end
    
end