classdef recognition
    
    properties (GetAccess = public)
        % GMM of voice, which will be compared
        gmmMain
        % Block - weight fitted gmms
        GMMS = {}
        % Block divided data
        DATA = {}
        % Amount of data blocks(states in HMM) 
        dataStates
    end
    
    methods (Access = public)
        % Shows Viterbi way for original person
        % dataMatrix - feature vectores, which describes original person
        % mapIteration - algorithm intensity
        % stateAmount -  Amount of data blocks(states in HMM) 
        function obj = recognition(dataMatrix, mapIteration, stateAmount)
             ubm = importdata('c:\Users\AirVan\Documents\MATLAB\SCS_Voici\UBM\mainUbmObj.mat');
             obj.gmmMain = gmmMaptest(ubm);
             obj.dataStates = stateAmount;
             obj.gmmMain = obj.gmmMain.fitByMeans(dataMatrix, mapIteration);
             obj.DATA = obj.gmmMain.divideBlock(dataMatrix,  obj.dataStates);
             obj.GMMS = obj.gmmMain.fitBlockWeights(dataMatrix,  obj.dataStates, mapIteration);
             %ViterbiD(obj.GMMS, obj.DATA);
             for i = 1:obj.dataStates
                fprintf('Identification [%i]: %i\n', i, simpleIdent(obj.GMMS{i}, dataMatrix));
             end
        end
        
        % Verifies data, using fitted GMMS from original person
        % dataMatrix - testing feature vectors
        function test = compare(obj, dataMatrix)
            otherDATA = obj.gmmMain.divideBlock(dataMatrix,  obj.dataStates);
            %ViterbiD(obj.GMMS, otherDATA);
             for i = 1:obj.dataStates
                fprintf('Identification [%i]: %i\n', i, simpleIdent(obj.GMMS{i}, dataMatrix));
             end
        end
    end
end