% (i) Function for retrieving the joint probability distribution for the
% parent nodes

function [CI] = parentCIs(ni, nodes, crossimpacts, states)

    % Input parameters
    % ni: child node index
    % nodes: All node (struct) elements in the Bayesian network
    % crossimpacts: Initial cross-impacts converted using formula (7) on
    % page (14)
    % states: A vector storing the number of states for each node
    
    % Output parameters
    % CI: Cross-impacts matrix that contains the impacts between node ni and
    % its parents
    
    CI = [];
    i = nodes(ni).CIindex;
    for p = nodes(ni).parents
        pi = nodes(p).CIindex;
        CI = [CI, crossimpacts(sum(states(1:i-1))+1:sum(states(1:i)), ...
            sum(states(1:pi-1))+1:sum(states(1:pi)))]; %#ok<AGROW>
    end
end