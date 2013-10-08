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
        % GMM - arrray of gmdistribution for model training
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
            % Initialization - up-triangle matrix filled by 0.5
            obj.transitA = triu(ones(obj.nStates));
            % Initialization, using normal probability density function
            obj.emisB = obj.initB(); 
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
        function likelihood = likelihoodCalc(obj, gmmDistr, Data)
            vectAmount = length(Data(:,1));
            vectProb = zeros(vectAmount, 1);
            gaussProb = zeros(gmmDistr.NComponents, 1);
            counter = 0;
            for i = 1:vectAmount
                for j = 1:gmmDistr.NComponents
                    gaussProb(j) =  gmmDistr.PComponents(j) * ...
                    mvnpdf(Data(i, :), gmmDistr.mu(j,:), gmmDistr.Sigma(:,:,1));
                    vectProb(i,1) = vectProb(i,1) + gaussProb(j,1);
                end
                counter = counter + log(vectProb(i,1));
            end
            likelihood = counter;
        end
    end
    
end

