% In this module we define the modified function for
% calculating the conditional distribution using the new approach.
%
% We solve the LSQ problem with Matlab's built in lsqlin solver

function [x, C0, resnorm,residual,exitflag,output,lambda]=ls_bayes_odds(C_I, m, scenp, states)

% Input parameters
% C_I: Cross-impact matrix
% m: Prior marginal probabilities
% scenp: Scenario probabilities
% states: Vector indicating the number of states for each node

% Example input for the function: 
% nodes(i).cdist = ls_bayes_odds(CI, probs{i}, pr, states([i, nodes(i).parents]));

% Output parameters
% x: Solution to the least-squares problem (estimated probabilities)
% CO: initial transformed cross-impact matrix
% resnorm: Residual norm of the solution
% exitflag: Exit condition of the optimization solver
% output: Output information from the optimization solver
% lambda: Lagrange multipliers from the optimization solver

% Initial checks:
% If any of the input matrices or vectors are empty, or if the length 
% of states is less than 2, return x as m reshaped into a column vector, 
% and C0 as an empty array.
resnorm=0;
residual=0;
if isempty(C_I) || isempty(m) || isempty(scenp) || length(states)<2
    x=m(:);
    C0=[];
    return
end

% start timer to measure the function's execution time
tic

% Convert the cross-impact statements into cross-impacts according to the 
% conversion formula 
C_I(isnan(C_I))=-4;
C0=sqrt(2.^C_I);
C0(C0==0.25)=0;
N=states(1); % current node's states
nparents=length(states)-1; % current node's parent nodes

constraints=[]; %Constraints array initialized
equals=[]; % Equals array initialized

skenprob2=scenp; % Alias for scenp

% Marginal constraints
%
% Here we build the linear equality constraint (1) of the objective function: 
% the sum of the estimated partial scenario probabilities conditioned on 
% adding a new uncertainty factor multiplied by the previous 
% estimates for the partial scenario probabilities must equal the marginal 
% probabilities first estimated 
% (i.e. \sum_{s \in S_{D_i}} q(k|s)q(s) = \hat{p}_k^i)
% 
% we append the constraint matrices with constraints corresponding to the 
% state at hand. These matrices are used in the linear equality 
% constraints in the following least squares problem.  
for mi=1:N 
    temp=zeros(1,N); 
    temp(mi)=1; 
    temp2=(skenprob2(:)*temp)';
    constraints=[constraints;temp2(:)']; % Append the resulting row vector 
    % into the constraints matrix as a new row below previous rows
    equals=[equals;m(mi)]; % Append the equals matrix with a new row 
    % corresponding to the prior marginal probabilities of the mi^th state
end

% Additional linear equality constraints built with a helpper function that
% adjusts for varying number of states in the states vector
% 
% these are used for the constraint (2) of the objective function: 
% the partial scenario probabilities conditioned on the uncertainty factor 
% k, must sum to one.
% (i.e. \sum{k=1}^{n_i}q(k|s_{D_i})=1)
[Aeq, beq] = conditional_distributions(states);
constraints=[constraints;Aeq];
equals=[equals;beq];

% Constructing the objective function
%
% Here, we construct the variables C and d corresponding to C*x=d linear
% system in the least squares. In our case, sum(q(k|s)q(s)-C*p_k*p_l, where
% q(s) corresponds to x. 
%
% Initialize c and d as empty arrays
C=[];
d=[];
% Extend scenp to match dimensions for element-wise multiplication
sceExt = ones(states(1),1)*scenp(:)';
for pi=1:nparents % looping over the parent nodes 
    for sp=1:states(pi+1) % and their states

        % projection of scenprob onto the current state of the parent node
        % (i.e. getting the probability p_l^j)
        psp = indexmagic2(pi,sp,states(2:end))'*skenprob2(:);

        % for each state s, construct rows C and d of the current uncertainty
        % factor.
        for s=1:N
            % The probability of the partial scenario extended with
            % outcome k (i.e. q(k|s) ) is 
            C=[C;(indexmagic2(1,s,states).*indexmagic2(1+pi,sp,states).*sceExt(:))'];
                 
            % since, C0 contains the cross-impacts and
            % m(s) is the probability of the state s of the uncertainty
            % factor pi representing \hat{p}_k^i. psp contains the 
            % probability \hat{p}_l^j and represents this.
        
            % Introducing additional variables resembling the notation in the paper 
            %C_kl_ij=C0(s,sum(states(1:pi-1))+sp); %  represents \hat{C}_{kl}^{ij}
            
            % Workaround to fix the issue of the index exeeding the proper
            % value
            col_index = sum(states(1:pi-1)) + sp;
            if col_index <= size(C0, 2)
                C_kl_ij = C0(s, col_index);     
            end
            p_k_i=m(s); % m(s) represent \hat{p}_k^i
            p_l_j=psp; % psp represent \hat{p}_l^j
           
            % New interpretation: the derived joint probability
            p_kl_ij = C_kl_ij*p_k_i*p_l_j / (1-p_k_i+C_kl_ij*p_k_i);

            % Compare: the derived joint probability
            %p_kl_ij = C_kl_ij*p_k_i*p_l_j;
            %d=[d; C0(s,sum(states(1:pi-1))+sp)*m(s)*psp]; %#ok<*AGROW> 

            d=[d; p_kl_ij];
        end
    end
end

%Append a small regularization term to C and d to ensure numerical stability
C(end+1,:)=1E-3;
d(end+1)=0;
%C=sparse(C);
%constraints=sparse(constraints);

% Solving the LSQ problem
% Where C and d define the linear equations
% constraints and equals define the equality constraints 
[x,resnorm,residual,exitflag,output,lambda]=lsqlin(C,d,[],[],constraints,equals,zeros(1,prod(states)),ones(1,prod(states)));

%sum(C*x)
%sum(d)
%residual
%x=reshape(x,states);

% End the timer
toc
end


% Helper function

% function aiding in the construction of the second constraint
% in the objective function. This adjusts the constraints correctly based on 
% the number of states for each uncertainty factor 
function [Aeq, beq] = conditional_distributions(states)

% npstates is the number of possible combinations for the different states 
% of all uncertainty factors excluding the first uncertainty factor
npstates=prod(states(2:end));

% beq is a vector of ones on the right hand side to enforce the sum-to-one constraint
beq=ones(npstates,1);

% Aeq is one only for those  zero matrix to ensure that a conditional probability of a partial
% scenario conditioned on the realization of a state of a new uncertainty
% factor is greater than or equal to zero
Aeq=zeros(npstates,prod(states));
for i=1:npstates
    % for each row, a block of states(1) consecutive columns are set to 1 
    % The matrix resembles a diagonal matrix whose diagonal values are 1.
    % This way each column has only one value that is 1 and all other
    % values are zero. Therefore the sum 
    % previous comb times first
    % uncertainty factors states + 1 to this combs times first uncertainty
    % factors states 
    Aeq(i,((i-1)*states(1)+1):(i*states(1)))=1;
end
%disp(Aeq)
end
