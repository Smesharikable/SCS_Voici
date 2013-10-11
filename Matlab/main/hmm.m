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
        end
        
        % Create transition probability matrix
        % with probability 0.5 to taking a self-loop
        % and probability 0.5 to going to the next state
        function matA = initA(obj)
            matA = eye(obj.nStates);
            for i = 1:obj.nStates - 1
                matA(i, i + 1) = 0.5;
            end
        end
        
        % Counting likelihood matrix matB consisting
        % b(i,j) i - vectorBlock, j state
        function matB = initB(obj)
            matB = zeros(obj.nStates);
            for i = 1:obj.nStates
                for j = 1:obj.nStates
                    % Check this later
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
    end
    
end

