% In this module we implement a function for producing Bayes-network using
% th GeNIe-software https://www.bayesfusion.com/genie/

function succesful = genie_parser(nodes,states,fname)
    
    % Input parameters
    % nodes: All node (struct) elements in the Bayes network
    % states: A vector storing the number of states for each node
    
    % Output parameters
    % Unit

    load("geniefilestructs.mat","smile","genienode");
    if isempty(fname)
        fname="testi2.xdsl";
    end
    
    for i=1:length(nodes)
       
        smile.nodes.cpt(i).idAttribute="Node"+i;

        for j=1:states(i)
            smile.nodes.cpt(i).state(j).idAttribute="State"+j;
        end

        if isempty(nodes(i).parents)
            cdist=nodes(i).cdist;
        else
            cdist=reshape(nodes(i).cdist,states([i, nodes(i).parents]));
            cdist=permute(cdist,[1,(length(nodes(i).parents):-1:1)+1]);
        end

        smile.nodes.cpt(i).probabilities = strjoin(string(cdist(:)));

        if ~isempty(nodes(i).parents)
            smile.nodes.cpt(i).parents="Node"+strjoin(string(nodes(i).parents)," Node");
        end
        
        smile.extensions.genie.node(i)=genienode;
        smile.extensions.genie.node(i).idAttribute = "Node"+i;
        smile.extensions.genie.node(i).name=nodes(i).name;
        smile.extensions.genie.node(i).position="100 "+ (200*i-100) + " 200 " + (200*i-40);
        
    end
    
    writestruct(smile,fname,"FileType","xml","StructNodeName","smile")
    succesful=true;

end
