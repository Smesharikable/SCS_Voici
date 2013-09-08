function ViterbiD = ViterbiD(GMMS, DATA)    
    % Modeling for D-dimensional case

    Start = cputime; % Measure of time

    % ------ What we have: ------- %
    D = GMMS{1}.NDimensions;
    component = max(size(DATA)); % count parts
    feature = length(DATA{1}); % count of vector-features in single MAP-model 
    countFeature = component*feature; % count of all vector-features
    gaussian = GMMS{1}.NComponents; % count of gaussians in single MAP-model
    % a = 0; b = 15; % limit of the range (expectation value)

    fprintf('Calculating...\n');

    % ------ initialize of MAP-models ------- %
    % Map(component,1) = struct('set', '', 'features', '');
    % set(gaussian,1) = struct('mu', '', 'sigma', '', 'w', '');
    % features(feature,1) = struct('value', '');
    % 
    % fprintf('Waiting...\n');
    % 
    % for i = 1:1:component
    %     for j = 1:1:gaussian
    %         for k = 1:1:gaussian
    %             Map(i).set(j).mu = a + b*rand(D,1); % MU is D-dimensional vector
    %         end
    %         for k = 1:1:gaussian % matrix of covariance
    %            Map(i).set(j).sigma = diag(cov(-1 + 2 * randn(D)))';
    %         end
    %         for k = 1:1:gaussian
    %             Map(i).set(j).w = 1/gaussian; % equal weights
    %         end
    %     end
    %     for j = 1:feature
    %         Map(i).features(j).value = a + b*rand(D,1); % vector-features
    %     end
    % end

    % probability vector of initial states
    pi = zeros(component,1);
    pi(1,1) = 1;

    % matrix of transitions
    A = .5 * eye(component);
    for i = 1:(component-1)
        A(i,i+1)=.5;
    end

    % calculation of emission probability
    B = zeros(component); 
    for i = 1:component
       B(i,i) = B_emis(i, GMMS{i}, DATA);
    end


    % ------ algorithm Viterbi -------%
    Pr = zeros(component); % target probability
    T = component; % time moment
    TIndex = zeros(component);
    aaa = zeros(component,1);

    for i = 1:component
        Pr(i,1) = pi(i)* B(i,1);
        TIndex(1,1) = 1;
    end

    for i = 2:T
        for j = 1:component
            for k = 1:component
                    aaa(k,1) = Pr(k,i-1) + log(A(k,j));
                    % aaa(k,1) = Pr(k,i-1) * A(k,j)
            end
            if (i == 2) && (j == 1)
    %             fprintf(' %i \n', aaa);
            end
                Pr(j,i) = myMax(aaa) + B(j,i); % on i-step in j-state
                TIndex(j,i) = argmax(aaa);
        end
    end

    for k = 1:component
       aaa(k,1) = Pr(k,T);
    end

    X(T) = argmax(aaa);

    for i = (T-1):-1:2
        X(i) = TIndex(X(i+1), i+1);
    end

    fprintf('States: 1 ');
    fprintf('%i ', X(1:component-1));
    fprintf('\n');
    fprintf('Probability: %i ', myMax(aaa));
    fprintf('\n');

    % End measure of time
    fprintf('Time: %i sec \n', cputime - Start);
end