% Import util functions
thisDir     = fileparts(mfilename('fullpath'));    % ...\examples\green_product
projectRoot = fileparts(fileparts(thisDir));       % go up two levels
addpath(fullfile(projectRoot, 'utils'));

% Script for printing network example on page 28.

% Define nodes in the Bayesian network
nodes=struct("name","1. Global industry growth","parents",[],"children",[2 3 4 5 7 8],"CIindex",1);
nodes(2)=struct("name","2. Printing speed compared to present","parents",[1],"children",[3 4 11],"CIindex",4);
nodes(3)=struct("name","3. Printing cost compared to present","parents",[1 2],"children",[4 8 9 10],"CIindex",3);
nodes(4)=struct("name","4. Finnish industry growth","parents",[1 2 3],"children",[5 6 7 8 9 10 11],"CIindex",2);
nodes(5)=struct("name","5. Number of graduates with 3D-printing expertise","parents",[1 4],"children",[9],"CIindex",5);
nodes(6)=struct("name","6. Legal regulation of 3D-printing in Finland","parents",[1 4],"children",[8 9 10],"CIindex",6);
nodes(7)=struct("name","7. Standardization of processes and models","parents",[1 4],"children",[9 11],"CIindex",7);
nodes(8)=struct("name","8. Use of 3D-printed objects in FDF","parents",[3 4 6],"children",[9 10 11],"CIindex",8);
nodes(9)=struct("name","9. FDF access to 3D-printing model files","parents",[3 4 5 6 7 8],"children",[10 11],"CIindex",9);
nodes(10)=struct("name","10. FDF 3D-printing spare parts in peacetime","parents",[3 4 6 8 9],"children",[11],"CIindex",10);
nodes(11)=struct("name","11. FDF 3D-printing spare parts in crisis times","parents",[2 4 7 8 9 10],"children",[],"CIindex",11);

% Read the cross-impact matrix from a CSV file
crossimpacts = readmatrix('cross-impacts-3d.csv',"FileType","text");
crossimpacts = crossimpacts+crossimpacts';

% Probabilities for each state in the same order as the nodes
probs={[0.3 0.5 0.2], [0.4 0.5 0.1], [0.5 0.4 0.1], [0.3 0.5 0.2], [0.2 0.6 0.2], [0.05 0.9 0.05], [0.35 0.45 0.2], [0.1 0.5 0.4], [0.2 0.7 0.1], [0.7 0.29 0.01], [0.45 0.45 0.1]};

% States is a vector indicating the number of possible states for each uncertainty factor
states=ones(1,length(nodes))*3;

% Initialize the joint probability distribution
jointpd=[];

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
% genie_parser(nodes,states,"printing_network.xdsl");