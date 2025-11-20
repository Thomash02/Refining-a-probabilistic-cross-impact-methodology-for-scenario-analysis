% (iv) Function that outputs a multidimensional logical array that sets the
% positions corresponding to the node and state for a specific state of a
% given node

function [out] = indexmagic2(node,state,states)

    % Input parameters
    % node: the index of the node
    % state: the specific state of the node
    % states: vector containing the number of states for each node
    
    % Output
    % out: A multidimensional logical array
    
    % Error checking
    if(states(node)<state)
        error("Invalid state input. State exceeds the number of states for that node.")
    end

    % initialize output: all elements set to false in the same shape as states
    % multidimensional if the length of states is greater than one, else column
    % vector
    if length(states)>1
        out = false(states);
    else
        out = false(states,1);
    end
    
    % Below is a command string that when executed, converts all of the
    % elements in the output vector to false, except the index corresponding 
    % to the given node, which is converted into the given state specified by
    % the state input. This happens, because the array is first
    % initialized to false by the if statements, and the for loops below sets 
    % ":" for all of other elements except the element of the input node 
    % index, which is set to be the state by "call=call+state"
    call="out(";

    for i=1:(node-1)
        call=call+":,";
    end

    call=call+state;

    for i=(node+1):length(states)
        call=call+",:";
    end

    call=call+")=true;";
    eval(call); % execute command string: specified positions become true
    out=out(:); % convert into a column vector 

end
