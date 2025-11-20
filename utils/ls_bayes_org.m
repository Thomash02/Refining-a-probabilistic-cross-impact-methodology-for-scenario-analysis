% In this module we define the original function for calculating the
% conditional distribution using the least-squares approach based on
% cross-impact coefficients.
%
% We solve the LSQ problem with Matlab's built-in lsqlin solver.
%
% Compared to ls_bayes_odds, this version uses the original interpretation
% of the cross-impact coefficients, where the derived joint probability is
% modeled as:
%       p_kl_ij = C_kl_ij * p_k_i * p_l_j

function [x, C0, resnorm, residual, exitflag, output, lambda] = ls_bayes(C_I, m, scenp, states)

% Input parameters
% C_I    : Cross-impact matrix (log2-scale, or impact statements) between
%          the current node and its parent nodes.
% m      : Prior marginal probabilities of the states of the current node.
% scenp  : Scenario probabilities of the parent nodes (joint distribution
%          over the parent configurations).
% states : Vector indicating the number of states for each node.
%          - states(1)     = number of states of the current node
%          - states(2:end) = number of states of each parent node

% Output parameters
% x        : Solution to the least-squares problem (estimated conditional
%            probabilities q(k | s), vectorized).
% C0       : Transformed cross-impact matrix derived from C_I.
% resnorm  : Residual norm of the LSQ solution.
% residual : Residual vector of the LSQ solution.
% exitflag : Exit condition flag from lsqlin.
% output   : Structure with information returned by lsqlin.
% lambda   : Structure containing the Lagrange multipliers at the solution.

% Initial checks
% If any of the input matrices or vectors are empty, or if the number of
% nodes (length(states)) is less than 2 (no parents), we simply return the
% prior marginal probabilities m reshaped as a column vector. In this case
% no optimization is needed.
resnorm  = 0;
residual = 0;
if isempty(C_I) || isempty(m) || isempty(scenp) || length(states) < 2
    x   = m(:);
    C0  = [];
    return
end

% Start timer to measure the function's execution time
tic

% Transform cross-impact statements into cross-impact coefficients C0
% The original cross-impact matrix C_I is given in a format where
% cross-impacts are expressed through powers of 2. We convert them using:
%       C0 = sqrt(2 .^ C_I)
C0 = sqrt(2.^C_I);

% Basic dimensions
N        = states(1);         % number of states of the current node
nparents = length(states) - 1; % number of parent nodes

% Initialize linear equality constraint matrices
% constraints * x = equals
constraints = [];
equals      = [];

% Alias for scenario probabilities
skenprob2 = scenp;

% 1) Marginal constraints
% Here we build the linear equality constraints enforcing that the
% estimated conditional distributions reproduce the given marginal
% probabilities of the current node:
%
%   sum_{s in S_{D_i}} q(k | s) * q(s) = m(k)
%
% where:
%   - q(k | s) are the unknown conditional distributions to be estimated
%   - q(s) = scenp(s) is the probability of the parent configuration s
%   - m(k) is the prior marginal probability of the k-th state of the
%     current node
%
% We construct one equality constraint per state of the current node.
for mi = 1:N
    temp      = zeros(1, N);
    temp(mi)  = 1;                         % select state mi
    temp2     = (skenprob2(:) * temp)';    % combines temp with scenario probs
    constraints = [constraints; temp2(:)']; % append as a new row
    equals      = [equals; m(mi)];         % corresponding marginal prob
end

% 2) Normalization constraints for conditional distributions
% The conditional distributions must sum to one for each configuration of
% the parents:
%
%   sum_{k = 1}^{n_i} q(k | s_{D_i}) = 1
%
% The helper function conditional_distributions(states) constructs these
% constraints in a scalable way, taking into account the number of states
% for each node.
[Aeq, beq] = conditional_distributions(states);
constraints = [constraints; Aeq];
equals      = [equals;      beq];

% 3) Constructing the objective function
% We now construct the least-squares objective based on the cross-impact
% structure:
%
%   min_x || C * x - d ||_2^2
%
% In this formulation, C and d embody the relationship
%   C_kl_ij * m(s) * p_l_j ≈ q(k | s) * q(s)
%
% where:
%   - C_kl_ij comes from the transformed cross-impact matrix C0
%   - m(s) is the marginal probability of state s of the current node
%   - p_l_j is the marginal probability of state l of parent j
%
% The matrix C and vector d are constructed by looping over parent nodes and
% their states, and then over the states of the current node.

C = [];
d = [];

% Extend scenario probabilities to match dimensions for element-wise
% multiplication when building the rows of C
sceExt = ones(states(1), 1) * scenp(:)';

for pi = 1:nparents                     % loop over parent nodes
    for sp = 1:states(pi+1)             % loop over states of parent pi
        % Projection of the scenario probabilities onto the current state
        % of the parent node pi:
        %   psp = p_l^j (marginal of parent j in state l)
        psp = indexmagic2(pi, sp, states(2:end))' * skenprob2(:);

        % For each state of the current node
        for s = 1:N
            % Each row of C corresponds to the partial scenario probability
            % extended with outcome k, combined with a particular parent state.
            % indexmagic2(1, s, states) selects the s-th state of the current
            % node; indexmagic2(1+pi, sp, states) selects the sp-th state of
            % the pi-th parent. sceExt(:) provides the scenario probabilities.
            C = [C; (indexmagic2(1, s, states) .* ...
                     indexmagic2(1 + pi, sp, states) .* ...
                     sceExt(:))'];

            % Corresponding element of d is the product of:
            %   - cross-impact coefficient from C0
            %   - marginal probability of state s, m(s)
            %   - parent marginal probability psp
            d = [d; C0(s, sum(states(1:pi-1)) + sp) * m(s) * psp]; %#ok<*AGROW>
        end
    end
end

% 4) Regularization for numerical stability
% Append a small regularization row to C and d to improve numerical
% stability of the least-squares problem.
C(end+1, :) = 1E-3;
d(end+1)    = 0;

% 5) Solve the constrained least-squares problem
% We use lsqlin to solve:
%
%   min_x || C * x - d ||_2^2
%   subject to:
%       constraints * x  = equals
%       0 <= x <= 1      (probabilities)
%
[x, resnorm, residual, exitflag, output, lambda] = lsqlin( ...
    C, d, [], [], ...
    constraints, equals, ...
    zeros(1, prod(states)), ... % lower bounds
    ones(1,  prod(states)));    % upper bounds

% Example reshaping if needed to interpret as a multidimensional conditional
% probability table:
rs = reshape(x, states); %#ok<NASGU> % kept for possible debugging / inspection

% Stop timer
toc
end

% Helper function
% conditional_distributions(states)
%
% Constructs the normalization constraints ensuring that for each
% combination of parent states, the conditional distribution of the
% current node sums to one.
%
% Specifically, it builds Aeq and beq such that:
%   Aeq * x = beq
% encodes:
%   sum_{k = 1}^{states(1)} q(k | s_{D_i}) = 1
% for all npstates = prod(states(2:end)) possible parent configurations.

function [Aeq, beq] = conditional_distributions(states)

% npstates = number of possible combinations of all parent states
npstates = prod(states(2:end));

% Right-hand side of the constraints: each conditional distribution sums to 1
beq = ones(npstates, 1);

% Aeq is initialized as a zero matrix, and for each parent configuration we
% set a consecutive block of states(1) entries to 1, corresponding to all
% states of the current node for that parent configuration.
Aeq = zeros(npstates, prod(states));
for i = 1:npstates
    % For row i, columns from (i-1)*states(1)+1 to i*states(1) are set to 1.
    % This ensures that:
    %   q(1 | s_i) + q(2 | s_i) + ... + q(states(1) | s_i) = 1
    Aeq(i, ((i-1)*states(1) + 1):(i*states(1))) = 1;
end

end
