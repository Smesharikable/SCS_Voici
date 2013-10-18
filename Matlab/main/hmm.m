classdef hmm
    
    properties (GetAccess = public)
        % Transition probability matrix of
        % a(i,j) from i-state to j-state
        transitA
        % State observation matrix of
        % b(i,j) of likelihood of i-vectorBlock in j state
        emisB
        % Viterbi matrix
        matVit
        % Viterbi path
        pathVit
        % GMM - array of gmdistribution for model training
        GMM
        % DATA - block divided data for training
        DATA
        % Amount of possible states 
        nStates
    end
    
    methods (Access = public)
        % GMM - block fitted GMMS cell-array
        % DATA - feature vectores divided by block
        function obj = hmm(GMMS, DATA)
            obj.GMM = GMMS;
            obj.DATA = DATA;
            obj.nStates = length(DATA);
            % Initialization - look at the initA() description
            obj.transitA = obj.initA();
            % Initialization, using normal probability density function
            obj.emisB = obj.initB();
            obj.matVit = zeros(obj.nStates, obj.nStates);
            [obj.pathVit, obj.matVit] = obj.findPath();
        end
        
        % Create transition probability matrix
        % with probability 0.5 to taking a self-loop
        % and probability 0.5 to going to the next state
        function matA = initA(obj)
            matA = .5 * eye(obj.nStates);
            for i = 1:obj.nStates - 1
                matA(i, i + 1) = 0.5;
            end
        end
        
        % Counting likelihood matrix matB consisting
        % b(i,j) i - GMM(i) on vectorBlock(j), j - vectorBlock number
        function matB = initB(obj)
            matB = zeros(obj.nStates);
            for i = 1:obj.nStates
                for j = 1:obj.nStates
                    matB(i,j) = obj.likelihoodCalc(obj.GMM{i}, obj.DATA{j});
                end
            end   
        end
        
        % State observation likelihood calculation
        function likelihood = likelihoodCalc(~, gmmDistr, Data)
            vectAmount = length(Data(:,1));
            vectProb = 0;
            counter = 0;
            for i = 1:vectAmount
                for j = 1:gmmDistr.NComponents
                    vectProb = vectProb + gmmDistr.PComponents(j) * ...
                    mvnpdf(Data(i, :), gmmDistr.mu(j,:), gmmDistr.Sigma(:,:,1));
                end
                counter = counter + log(vectProb);
            end
            likelihood = counter;
        end
        
        % Finding Viterbi path for current data
        function [path, Viterbi] = findPath(obj)
            path = zeros(obj.nStates);
            Viterbi = zeros(obj.nStates);
            % Initialization step
            for i = 1:obj.nStates
                 Viterbi(1, i) = obj.transitA(1, i) * obj.emisB(1, i); 
            end
            % Recursion step
            for time = 1:obj.nStates
                for state = 1:obj.nStates
                    fprintf('Cycle time/state = %d / %d \n', time, state);
                    [Viterbi(state, time), path(state, time)]  = obj.maxPrevVit(Viterbi, state, time);
                    Viterbi
                end
            end
        end
        
        % Recousrion step in Viterbi
        function [maxValue, maxArg] = maxPrevVit(obj, Viterbi, state, time)
            fprintf('Recursion time/state = %d / %d \n', time, state);
            maxValue = -inf;
            maxArg = 0;
            if time <= 1
                maxValue = obj.transitA(1, state) * obj.emisB(state, 1);
                maxArg = 0;
                fprintf('Out recursion\n');
            else
                b = obj.emisB(state, time);
                for s = 1:obj.nStates
                    a = obj.transitA(s, state);
                    [vPrev, maxArg] = obj.maxPrevVit(Viterbi, s, time - 1);
                    vPrev = vPrev * a * b;
                    fprintf('vPrev = %d a = %d b = %d \n', vPrev, a, b);
                    if maxValue < vPrev 
                        maxValue = vPrev;
                    end
                end
            end
        end
        
    end
    
end

