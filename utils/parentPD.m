% (ii) Function for calculating the joint probability distribution for the
% parent nodes

function pr = parentPD(ni, nodes, pdist, states)

    % Input parameters
    % ni: child node index
    % nodes: All node (struct) elements in the Bayesian network
    % pdist: Joint probability density distribution
    % states: A vector storing the number of states for each node

    % Output parameters
    % pr: Joint probability distribution of the parent nodes

    pr = pdist;
    states = states(1:ni-1);
    
    for i = ni-1:-1:1
        if ~any(i == nodes(ni).parents)
            states(i) = [];
            if length(states) > 1
                pr = reshape(sum(pr, i), states);
            else
                pr = sum(pr, i)';
            end
        end
    end
end