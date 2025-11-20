% Import util functions
thisDir     = fileparts(mfilename('fullpath'));    % ...\examples\green_product
projectRoot = fileparts(fileparts(thisDir));       % go up two levels
addpath(fullfile(projectRoot, 'utils'));

% In this module we implement the main loop for simulating the performance
% of both approaches. This module consists of the main simulation loop and
% the necessary aiding functions. 

m = 100;                 % Number of runs

% Initialize jointpd storages
Joint_pd_Odds     = zeros(177147, m);
Joint_pd_Original = zeros(177147, m);

% initialize ECDFs
CDF_Original_all = zeros(177147, m);
CDF_Odds_all     = zeros(177147, m);
ET_run           = zeros(1, m);

% The main simulation loop that runs m-times
%
% Each simulation run we generate a random Bayes-network with 11 
% uncertainty factors each having 3 outcomes. The maximum number of 
% child-parent connections is 6, and the connections are randomly
% generated. The cross-impact statements are also randomly generated

for H = 1:m % run m-simulations
    
    % initialize a random Bayes-network
    [cross_impacts, probs, nodes, states] = random_bayes_network();
    
    % initialize join probability density storages for both methods
    jointpd_Original = [];
    jointpd_Odds = [];
    
    ETO = zeros(1, length(nodes)); % Store run times for both methods
    ETOD = zeros(1, length(nodes));

    % Perform the optimization and store results
    for i = 1:length(nodes)

        % Retrieve Cross-impacts with parents
        CI = parentCIs(i, nodes, cross_impacts, states);
        
        % Retrieve joint probability distribution for parents
        pr = parentPD(i, nodes, jointpd_Original, states); % A)
        pr_Odds = parentPD(i, nodes, jointpd_Odds, states); % B)
        
        % Original method A)
        tic; % start performance clock

        % Calculate the conditional probability distributions for 
        % the i-th uncertainty factor given the states of its parents
        nodes(i).cdist = ls_bayes_org(CI, probs(3*i-2:3*i), pr, states([i, nodes(i).parents]));    
        
        % Append time storage vector
        ETO(i) = toc;
        
        % New method B) (Using ls_bayes30 for some reason????)
        tic; % start performance clock

        % Calculate the conditional probability distributions for 
        % the i-th uncertainty factor given the states of its parents
        nodes(i).cdist_Odds = ls_bayes_odds(CI, probs(3*i-2:3*i), pr_Odds, states([i, nodes(i).parents]));
        
        % Append time storage vector
        ETOD(i) = toc;
        
        % Store resulting joint probability distributions
        jointpd_Original = updatePD(i, nodes, jointpd_Original, states, nodes(i).cdist);
        jointpd_Odds = updatePD(i, nodes, jointpd_Odds, states, nodes(i).cdist_Odds);
    end
    % Formulate the CDFs
    [CDF_Odds, CDF_Original] = empirical_cumulative_distribution_functions(jointpd_Original, jointpd_Odds);
    
    % Store CDFs
    CDF_Original_all(:, H) = CDF_Original(:);
    CDF_Odds_all(:, H)     = CDF_Odds(:);
    ET_run(H)              = mean(ETOD ./ ETO);

    % Store jointpds
    Joint_pd_Original(:, H) = jointpd_Original(:)';
    Joint_pd_Odds(:, H) = jointpd_Odds(:)';
end

% Funtion for generating a random_bayes_network with 11 
% uncertainty factors each having 3 outcomes. The maximum number of 
% child-parent connections is 6, and the connections are randomly
% generated. The cross-impact statements are also randomly generated
function [cross_impacts, probs, nodes, states] = random_bayes_network()
    
    % Output parameters
    % crossimpacts: Initial cross-impacts converted using formula (7) on
    % page (14)
    % probs: generated marginal probabilites
    % nodes: All node (struct) elements in the Bayesian network
    % states: A vector storing the number of states for each node
    
    % initalize states
    states = ones(1, 11) * 3; % Each uncertainty factor has 3 states

    % Simulating 11 uncertainty factors with 3 outcomes for each
    cross_impacts = zeros(33, 33);
    
    % Generating random cross-impact statements between [-3,3]
    for i = 1:11
        for j = i+1:11            
            for k = 1:3 
                for l = 1:3
                    cross_impacts((i-1)*3 + k, (j-1)*3 + l) = randi([-3, 3]);
                end
            end
        end
    end
    cross_impacts = cross_impacts + cross_impacts';
    
    % Generate random probabilities for 11 uncertainty factors with 3 outcomes each
    probs = zeros(1, 33);
    for i = 1:11
        a = rand(1, 3);
        a = a / sum(a);  % Normalize probabilities
        probs(3*i-2:3*i) = a;  % Assign probabilities to probRNG
    end
    
    % Initialize node structures
    nodes = struct("name", "", "parents", [], "children", [], "CIindex", 1);
    
    for i = 1:11
        nodes(i).CIindex = i;
        nodes(i).name = i;
        
        % Determine parents and children, ensuring max 6 connections
        max_connections = 6;
        if i == 1
            num_parents = 0;
        else
            num_parents = randi([1, min(max_connections, i-1)]);
        end
        num_children = randi([0, max_connections - num_parents]);
        
        if num_parents > 0
            parents = sort(randperm(i-1, num_parents)); % Choose random parents
            nodes(i).parents = parents;
            
            % Add this node as a child to its parents
            for p = parents
                nodes(p).children = [nodes(p).children, i];
            end
        end
    end
end


% Function for generating the cumulative empirical density functions for
% both approaches
function [CDF_Odds, CDF_Original] = ...
    empirical_cumulative_distribution_functions(jointpd_Original, ...
    jointpd_Odds)

    % Input parameters
    % jointpd_Original: The obtained joint probability distribution for the
    % original method
    % jointpd_Odds: The obtained joint probability distribution for the new
    % method

    % Output parameters
    % CDF_Odds: The empirical cumulative distribution function for the new
    % method
    % CDF_Original: The empirical cumulative distribution function for the
    % original method
    
    Odds = jointpd_Odds(:)';
    Original = jointpd_Original(:)';
    
    % Sort Org from largest probability to smallest probability
    [sorted_Org, sortOrder] = sort(Original, 'descend');
    Odds = Odds(sortOrder); % org dist is the reference one for comparison
    Original=sorted_Org;
    x = (1:length(Odds));

    % Cumulative Distribution Functions (CDF)
    CDF_Odds = cumsum(Odds);
    CDF_Original = cumsum(Original);
end