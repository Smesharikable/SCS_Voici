classdef recognition
    
    properties (GetAccess = public)
        % GMM of voice, which will be compared
        gmmMain
        % Block - weight&mean fitted gmms
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
            % TODO: change constant string to parameter
            ubm = importdata('UBM\mainUbmObj.mat');
            obj.gmmMain = gmmMap(ubm);
            obj.dataStates = stateAmount;
            obj.gmmMain = obj.gmmMain.fitByMeans(dataMatrix, mapIteration);
            obj.DATA = obj.gmmMain.divideBlock(dataMatrix,  obj.dataStates);
            obj.GMMS = obj.gmmMain.fitBlockWeights(dataMatrix,  obj.dataStates, mapIteration);
            HMM = hmm(obj.GMMS, obj.DATA);
            fprintf('Emis matrix : GMMs x Data \n');
            HMM.emisB
        end
        
        % Verifies data, using fitted GMMS from original person
        % dataMatrix - testing feature vectors
        function dataToTest = compare(obj, dataMatrix)
            dataToTest = obj.gmmMain.divideBlock(dataMatrix,  obj.dataStates);
        end
    end
end