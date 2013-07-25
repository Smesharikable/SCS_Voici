classdef gmmaptest
    
    properties (Access = public)
        relativeCoeff = 16;
        isChanges = false;
    end
    
    properties (GetAccess = public, SetAccess = protected)
        gmmin
        gmmout
        pstat = [];
        Exstat = [];
        newp = [];
        newmu = [];
        newsigma = [];
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
        
        function obj = set.isChanges(obj, isChanges)
            if isscalar(isChanges) && islogical(isChanges)
                obj.isChanges = isChanges;
            end
        end
        
    end
    
    methods (Access = public)
        
        % data - matrix n-by-d, where n is the number of observations,
        %   and d is the dimensions of the data
        function test = fitbymeans(obj, data, iterations)
            outGmm = obj.gmmin;
            T = length(data(:, 1));
            oldp = obj.gmmin.PComponents;
            oldmu = obj.gmmin.mu;
            %oldSigma = obj.gmmin.Sigma;
            Pr = posterior(outGmm, data);
            for k = 1:iterations
%                 obj.newp = zeros(1, obj.gmmin.NComponents);
%                 obj.newmu = zeros(obj.gmmin.NComponents, obj.gmmin.NDimensions);
%                 obj.coef = zeros(1, obj.gmmin.NComponents);
%                 p = zeros(1, obj.gmmin.NComponents);
%                 Pr = zeros(T, obj.gmmin.NComponents);
%                 compDensity = zeros(1, obj.gmmin.NComponents);
%                 compDensityMult = zeros(1, obj.gmmin.NComponents);
%                 for t = 1:T
%                     for i = 1:obj.gmmin.NComponents
%                         for j = 1:obj.gmmin.NComponents
%                             compDensityMult(j) = 1 / ((2 * pi)^(obj.gmmin.NDimensions/2) * (det(oldSigma(j)^(0.5))));
%                             compDensity(j)= compDensityMult(j) * exp(-0.5 * (data(t,:) - oldmu(j,:)) * inv(oldSigma(j)) * (data(t,:) - oldmu(j,:))');
%                         end;
%                         summ = sum(oldp(:) .* compDensity(:));
%                         p(i) = (oldp(i) * compDensity(i)) / summ;
%                         Pr(t,i) = p(i);
%                     end
%                 end
                
                for i = 1:obj.gmmin.NComponents
                    obj.pstat(i) = sum(Pr(:,i));
                    obj.Exstat(i, :) = (Pr(:, i)' * data) / obj.pstat(i);
                    obj.coef(i) = obj.pstat(i)/ (obj.pstat(i) + obj.relativeCoeff);
                    obj.newp(i) = obj.coef(i) * obj.pstat(i) / T + (1 - obj.coef(i)) * oldp(i);
                    obj.newmu(i, :) = obj.coef(i) * obj.Exstat(i, :) + (1 - obj.coef(i)) * oldmu(i, :);
                end
                outGmm = gmdistribution(obj.newmu, obj.gmmin.Sigma, oldp);
                oldmu = obj.newmu;
            end
            obj.gmmout = outGmm;
            obj.gmmin = outGmm;
            test = obj;
        end
        
        function test = fitbyweight(obj, data, iterations)
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
                oldp = obj.newp;
            end
            test = outGmm;
        end
        
        function [outG, outD] = makeOut(obj, data, intervalLength, iteration)
            amount = fix(length(data(:,1)) / intervalLength);
            outGmms = {};
            for i = 1:amount
                start = (i - 1)*intervalLength + 1;
                stop = i * intervalLength;
                indata = data(start:stop,:);
                outGmms{i} = obj.fitbyweight(indata, iteration);
                outData{i} = indata;
            end
            outG = outGmms;
            outD = outData;
        end
        
        function outGmm = EM(obj, data, iterations)
            outGmm = obj.gmmin;
            T = length(data(:, 1));
            N = obj.gmmin.NComponents;
            oldp = obj.gmmin.PComponents;
            oldmu = obj.gmmin.mu;
            oldSigma = obj.gmmin.Sigma;
            for k = 1:iterations
                obj.newp = zeros(1, obj.gmmin.NComponents);
                obj.newmu = zeros(obj.gmmin.NComponents, obj.gmmin.NDimensions);
                obj.coef = zeros(1, obj.gmmin.NComponents);
                obj.newsigma = zeros(1, obj.gmmin.NDimensions, obj.gmmin.NComponents);
                p = zeros(1, obj.gmmin.NComponents);
                Pr = zeros(T, obj.gmmin.NComponents);
                compDensity = zeros(1, obj.gmmin.NComponents);
                compDensityMult = zeros(1, obj.gmmin.NComponents);
                for i =1:N
                    for t = 1:T
                        for j = 1:N
                            compDensityMult(j) = 1 / ((2 * pi)^(obj.gmmin.NDimensions/2) * (det(oldSigma(j)^(0.5))));
                            compDensity(j)= compDensityMult(j) * exp(-0.5 * (data(t,:) - oldmu(j,:)) * inv(oldSigma(j)) * (data(t,:) - oldmu(j,:))');
                        end;
                        summ = sum(oldp(:) .* compDensity(:));
                        p(i) = (oldp(i) * compDensity(i)) / summ;
                        obj.newp(i) = obj.newp(i) + p(i);
                        obj.newmu(i,:) = obj.newmu(i,:) + p(i) * (data(t,:));
                        obj.newsigma(:,:,i) = obj.newsigma(:,:,i) + p(i) * (data(t,:) * (data(t,:)'));
                        Pr(t,i) = p(i);
                    end;
                    obj.newp(i) = obj.newp(i) / T;
                    obj.newmu(i,:) = obj.newmu(i,:) / sum(Pr(:,i));
                    obj.newsigma(:,:,i) = obj.newsigma(:,:,i) / sum(Pr(:,i)) - obj.newmu(i,:) * obj.newmu(i,:)'; 
                end 
                oldp = obj.newp
                oldmu = obj.newmu;
                oldSigma = obj.newsigma;
            end
            outGmm = gmdistribution(obj.newmu, obj.newsigma, obj.newp);
        end
    end
        
    methods
        function value = get.newp(obj)
            value = obj.newp;
        end
             
        function value = get.newmu(obj)
            value = obj.newmu;
        end 
        
        function value = get.gmmout(obj)
            value = obj.gmmout;
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
            ezcontour(@(x, y)pdf(outGmm, [x y]), [-15 15], [-15 15], 200)
        end
        
    end
    
end