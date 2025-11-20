% (iii) Function for updating the joint probability distribution with the
% updated conditional distribution

function pr = updatePD(ni, nodes, pdist, states, cdist)

    % Input parameters
    % ni: child node index
    % nodes: All node (struct) elements in the Bayesian network
    % pdist: Joint probability density distribution
    % states: A vector storing the number of states for each node
    % cdist: Conditional distribution solved using LSQ

    % Ouput parameters
    % pr: Updated joint probability distribution

    if ni == 1
        pr = cdist;
        return
    end
    parents = nodes(ni).parents;
    % Reshape the conditional distribution to include the states of the current node and its parents
    cdist = reshape(cdist, states([ni, parents]));

    % Add independent uncertainty factors to the conditional distribution
    for i = 1:ni-1
        if ~any(i == nodes(ni).parents)
            np = length(parents);
            % Create a permutation vector to rearrange dimensions
            perm = [1:i, np + 2, i + 1:np + 1];
            % Expand the conditional distribution to include the new dimension for the independent node
            temp = cdist(:) * ones(1, states(i));
            % Reshape and permute the updated conditional distribution
            cdist = permute(reshape(temp, states([ni, parents, i])), perm);
            % Update the parents list to include the new node
            parents = sort([parents, i]);
        end
    end

    % Calculate the updated joint probability distribution
    skprob = ones(states(ni), 1) * pdist(:)';
    newScenarioP = (skprob .* reshape(cdist, states(ni), []))';
    pr = reshape(newScenarioP, states([parents, ni]));

end