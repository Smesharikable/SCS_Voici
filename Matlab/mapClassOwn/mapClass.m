classdef mapClass
    
    properties (Constant)
        RelativeCoeff = 16;
        Components = 16;
        Dimension = 2;
        Obs = 300;
    end
    
    properties (GetAccess = public, SetAccess = protected)
        mData;
        mUBM;
        mOutGMM;
        mWeights;
        mMeans;
        mCoef;
        mpstat;
        mExstat;
    end
    
    methods (Static)
        function ubmModel = createUBM
            mu = rand(mapClass.Components, mapClass.Dimension) * 20 - 10;
            sigma = [2 0;0 .5];
            p = rand(1, mapClass.Components);
            ubmModel = gmdistribution(mu,sigma,p);
        end
        
        function input = createInput
            mu = rand(mapClass.Components, mapClass.Dimension) * 20 - 10;
            sigma = [2 0;0 .5];
            p = rand(1, mapClass.Components);
            gmm = gmdistribution(mu,sigma,p);
            input = random(gmm, mapClass.Obs);
        end   
    end
    
    methods
        function obj = mapClass
            obj.mWeights = zeros(1, mapClass.Components);
            obj.mCoef = zeros(1, mapClass.Components);
            obj.mExstat = zeros(mapClass.Components, mapClass.Dimension);
            obj.mpstat = zeros(1, mapClass.Components);
            obj.mMeans = zeros(mapClass.Components, mapClass.Dimension);
            obj.mUBM = mapClass.createUBM;
            obj.mData = mapClass.createInput;
        end
        
        function value = get.mWeights(obj)
            value = obj.mWeights;
        end
             
        function value = get.mMeans(obj)
            value = obj.mMeans;
        end
        
        function obj = set.mWeights(obj,value)
             obj.mWeights = value;
        end
         
        function obj = set.mMeans(obj,value)
             obj.mMeans = value;
        end 
        
        function value = get.mUBM(obj)
            value = obj.mUBM;
        end
             
        function value = get.mData(obj)
            value = obj.mData;
        end
        
        function obj = set.mUBM(obj,value)
             obj.mUBM = value;
        end
         
        function obj = set.mData(obj,value)
             obj.mData = value;
        end 
        
        function test = fitByMW(obj, data, iterations)
            outGmm = obj.mUBM;
            data = obj.mData;
            T = length(data(:, 1));
            oldp = obj.mUBM.PComponents;
            oldmu = obj.mUBM.mu;
            for k = 1:iterations
                Pr = posterior(outGmm, data);
                obj.mWeights = zeros(1, mapClass.Components);
                obj.mMeans = zeros(mapClass.Components, mapClass.Dimension);
                obj.mCoef = zeros(1, mapClass.Components);
                for i = 1:mapClass.Components
                    obj.mpstat(i) = sum(Pr(:, i));
                    obj.mExstat(i, :) = (Pr(:, i)' * data) / obj.mpstat(i);
                    obj.mCoef(i) = obj.mpstat(i) / (obj.mpstat(i) + mapClass.RelativeCoeff);
                    obj.mWeights(i) = obj.mCoef(i) * obj.mpstat(i) / T + (1 - obj.mCoef(i)) * oldp(i);
                    obj.mMeans(i, :) = obj.mCoef(i) * obj.mExstat(i, :) + (1 - obj.mCoef(i)) * oldmu(i, :);
                end
                outGmm = gmdistribution(obj.mMeans, obj.mUBM.Sigma, obj.mWeights);
                oldmu = obj.mMeans;
                oldp = obj.mWeights;
            end
            obj.visualize(data, outGmm);
            obj.mOutGMM = outGmm;
            obj.mUBM = outGmm;
            test = obj;
        end
        
        function visualize(obj, data, outGmm)
            figure;
            hold on;
            scatter(data(:, 1), data(:, 2))
            oldmu = obj.mUBM.mu;
            for i = 1:mapClass.Components
                plot([oldmu(i, 1), obj.mExstat(i, 1)], [oldmu(i, 2), obj.mExstat(i, 2)]);
            end
            scatter(oldmu(:, 1), oldmu(:, 2))
            scatter(obj.mExstat(:, 1), obj.mExstat(:, 2))
            ezcontour(@(x, y)pdf(outGmm, [x y]), [-15 15], [-15 15], 200)
        end
        
    end
end