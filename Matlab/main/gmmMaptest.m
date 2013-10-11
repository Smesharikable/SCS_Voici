% TODO: change name to gmmMap
classdef gmmMaptest
    
    properties (Access = public)
        % Proofed coefficient
        relativeCoeff = 16;
    end
    
    properties (GetAccess = public)
        gmmin
        gmmout
        pstat = [];
        Exstat = [];
        newp = [];
        newmu = [];
        newsigma = [];
        coef =[];
        NComponents
    end
    
    
    methods
        % gmdist - GMM model(UBM), which keeps start values (mu, sigma, weight, Number of Components) 
        function obj = gmmMaptest(gmdistr)
            if nargin > 0
                if ~isa(gmdistr, 'gmdistribution')
                    error('Input parameter must be a gmdistribution');
                end
                obj.gmmin = gmdistr;
                obj.pstat = zeros(1, gmdistr.NComponents);
                obj.NComponents = gmdistr.NComponents;
                obj.Exstat = zeros(gmdistr.NComponents, gmdistr.NDimensions);
            end
        end
         
        function obj = set.gmmout(obj, gmm)
            if ~isa(gmm, 'gmdistribution')
                error('Input parameter must be a gmdistribution');
            end
            obj.gmmout = gmm;
        end
    end
    
    methods (Access = public)
        % Data - matrix n-by-d, where "n" is the number of observations,
        % and "d" is the dimension of the data vectors
        function test = fitByMeans(obj, data, iterations)
            % set's gmm from class
            outGmm = obj.gmmin;
            T = length(data(:, 1));
            oldp = obj.gmmin.PComponents;
            oldmu = obj.gmmin.mu;
            Pr = posterior(outGmm, data);
            for k = 1:iterations
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
        
        % Data - matrix n-by-d, where "n" is the number of observations,
        % and "d" is the dimension of the data vectors
        function test = fitByWeights(obj, data, iterations)
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
            obj.gmmout = outGmm;
            test = obj;
        end
        
        % Uses weight fit for every feature vector block
        % outGmm - modified gmm list, where gmm is different for every block
        % iteration - iteration amount for fitByWeights
        function [outGmm] = fitBlockWeights(obj, data, blockAmount, iteration)
            outData = obj.divideBlock(data, blockAmount); 
            for i = 1:blockAmount
                obj = obj.fitByWeights(outData{i}, iteration); 
                outGmms{i} = obj.gmmout;
            end
            outGmm = outGmms;
        end
        
        % Divides data into blocks = "intervalAmount" of feature vectors
        % outData - divided data blocks' list
        function [outData] = divideBlock(~, data, blocklAmount)
            amountOfVectors = fix(length(data(:,1)) / blocklAmount);
            outData = {};
            for i = 1:blocklAmount
                start = (i - 1) * amountOfVectors + 1;
                stop = i * amountOfVectors;
                indata = data(start : stop,:);
                outData{i} = indata;
            end
        end
        
        % EM algorithm testing
        % data - feature vectors
        % iterations - iterations of algorithm
        % outGmm - new trained GMM
        function outGmm = EM(obj, data, iterations)
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
                oldp = obj.newp;
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
        % Visualization for 2 dimensional case
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