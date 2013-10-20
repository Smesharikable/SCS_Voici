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
        % Initial & last step probaility
        pIL
    end
    
    methods (Access = public)
        % GMM - block fitted GMMS cell-array
        % DATA - feature vectores divided by block
        function obj = hmm(GMMS, DATA)
            obj.GMM = GMMS;
            obj.DATA = DATA;
            obj.nStates = length(DATA);
            obj.pIL = 1 / obj.nStates;
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
            obj.emisB
            obj.transitA
            % Initialization step
            for i = 1:obj.nStates
                Viterbi(i, 1) = -Inf;
                path(i, 1) = 0;
            end
            Viterbi(1, 1) = obj.emisB(1, 1);
            % Fulfill table step
            for time = 2:obj.nStates
                for state = 1:obj.nStates
                    curProb = Viterbi(state, time - 1) + ...
                            log(obj.transitA(state, state)) + obj.emisB(state, time);
                    if (state > 1)
                        prevProb = Viterbi(state - 1, time - 1) +...
                            log(obj.transitA(state - 1, state)) + obj.emisB(state, time);
                        if (prevProb > curProb)
                            Viterbi(state, time) = prevProb;
                            path(state, time) = state - 1;
                        else
                            Viterbi(state, time) = curProb;
                            path(state, time) = state;
                        end
                    else
                        Viterbi(state, time) = curProb;
                        path(state, time) = state;
                    end
                        
                    %maxProb = Viterbi(state, time - 1) + log(obj.transitA(1, state)) + obj.emisB(state, time);
                    %sMax = 1;
                    %for curState = 2:obj.nStates
                    %    curProb = Viterbi(state, time - 1) + log(obj.transitA(curState, state)) + obj.emisB(state, time);
                    %    if (curProb > maxProb)
                    %        maxProb = curProb;
                    %        sMax = curState;
                    %    end
                    %end
                    %Viterbi(state, time) = maxProb;
                    %path(state, time) = sMax;
                end
            end
            
            % Termination step
            maxProb = Viterbi(1, obj.nStates);
            lastS = 1;
            for s = 2:obj.nStates
                curProb = Viterbi(s, obj.nStates);
                if (curProb > maxProb)
                    maxProb = curProb;
                    lastS = s;
                end
            end
            
            Viterbi
            
            maxProb
            path
            path = obj.backtrace(path, lastS);
            path
        end
        
        % Replace to private
        function [backpath] = backtrace(obj, path, state)
            time = obj.nStates;
            backpath = zeros(1, time);
            backpath(1, time)= state;
            curPath = state;
            for i = time : -1 : 2
                curPath = path(curPath, i);
                backpath(1, i - 1) = curPath;
            end
        end
        
    end
    
end

