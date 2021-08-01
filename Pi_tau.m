function nextAction = Pi_tau(muGraph,sigmaGraph,traversalGraph,currentState,destinationNode,synchronous)
    persistent functionPath;
    
%     persistent outputPath;

if synchronous ~= 1
    nextAction = currentState;
    return;
end


    tau = 0.7;
    [nodeNum,~] = size(muGraph);
    if synchronous == 0     %0
        nextAction = currentState;
    elseif synchronous == 1     %1
        [~,n] = find(functionPath == currentState);
        if traversalGraph(currentState,functionPath(n + 1)) == 1
            nextAction = functionPath(n + 1);
            return;
        end
        for i = 1 : nodeNum
            for j = 1 : nodeNum
                if traversalGraph(i,j) < tau
                    muGraph(i,j) = 0;
                end
            end
        end
        c=sparse(muGraph);
%         currentState
        [~,functionPath,~] = graphshortestpath(c,currentState,destinationNode);
        [~,length] = size(functionPath);
        if length < 2
            nextAction = destinationNode;
        else
            nextAction = functionPath(2);
        end
    else
        [~,length] = size(functionPath);
        if length < 2
            nextAction = destinationNode;
        else
            nextAction = functionPath(2);
        end
    end
end