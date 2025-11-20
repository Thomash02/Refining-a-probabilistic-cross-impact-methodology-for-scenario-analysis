% Import util functions
thisDir     = fileparts(mfilename('fullpath'));    % ...\examples\green_product
projectRoot = fileparts(fileparts(thisDir));       % go up two levels
addpath(fullfile(projectRoot, 'utils'));

% Script for the illustrative example of exploring the possibilities of
% green products (see Table 3 on page 17).

% Define nodes in the Bayesian network
nodes = struct("name", "1. Geopolitics", "parents", [], "children", [2 3], "CIindex", 1);
nodes(2) = struct("name", "2. Economic development", "parents", [1], "children", [3], "CIindex", 2);
nodes(3) = struct("name", "3. Regulation", "parents", [1 2], "children", [], "CIindex", 3);

% Read the cross-impact matrix from a CSV file
crossimpacts = readmatrix('green-product_example.csv', "FileType", "text");
% Process the cross-impact matrix for further use
crossimpacts = 2 * log2(crossimpacts + crossimpacts');

% Probabilities for each state in the same order as the nodes
% Taken from table 2
probs = {[0.4 0.5 0.1], [0.3 0.4 0.3], [0.3 0.25 0.45]};

% States is a vector indicating the number of possible states for each uncertainty factor
states = [3 3 3];

% Initialize the joint probability distribution
jointpd = [];

% Iterate through each uncertainty factor one by one
for i = 1:length(nodes)
    % Retrieve cross impacts with the parents
    CI = parentCIs(i, nodes, crossimpacts, states);
    % Calculate the joint probability distribution for the parents only
    pr = parentPD(i, nodes, jointpd, states);
    % Calculate the conditional probability distributions for the i-th uncertainty factor
    % given the states of its parents
    nodes(i).cdist = ls_bayes_odds(CI, probs{i}, pr, states([i, nodes(i).parents]));
    % Update the joint probability distribution to include the i-th node
    % based on the calculated conditional distribution
    jointpd = updatePD(i, nodes, jointpd, states, nodes(i).cdist);
end

% The line below can be used to produce a Bayesian network in the GeNIe software
% genie_parser(nodes, states, "filename.xdsl");