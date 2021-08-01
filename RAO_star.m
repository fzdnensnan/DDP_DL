function nextAction = RAO_star(muGraph,sigmaGraph,traversalGraph,currentState,destinationNode,synchronous)
%UNTITLED6 此处显示有关此函数的摘要
%   此处显示详细说明
persistent path;
if currentState == destinationNode
    nextAction = destinationNode;
    return;
end

if synchronous == 2
    nextAction = currentState;
    return;
end
if synchronous == 1 && ~isempty(path)
    id = find(path == currentState);
%     [m,~] = max(id);
    if traversalGraph(currentState, path(id + 1)) == 1
        nextAction = path(id + 1);
        return;
    end
end
    


    omega = 2;
    root.index = 1;            %   
    root.parent = [];        %   parent
    root.nodeId = currentState;
    root.prob = 1;
    root.children = []; %   successor set     
    root.type = 0;          %   0-or  1-and
    root.state.nodeId = currentState;   %   nodeId in graph
    root.state.actionStatusList = updateActionList(traversalGraph,traversalGraph);%  uncertain actions status 0-U,1-T,2-A
    root.f = heuristicFunction(root.state.nodeId,muGraph,destinationNode);             %   state-value function
    root.status = false;    %   0-unsolved, 1-solved
    AOtree = root;
    node = root;
    while AOtree(1).status == false && AOtree(1).f ~= inf
        [~,sizeAOtree] = size(AOtree);
        if sizeAOtree > 1000
            nextAction = currentState;
            return;
        end
        node = selectMinCostNode(AOtree);
        [assumptMuGraph,assumptTraversalGraph] = calculateAssumptGraphFromActionStatusList(node.state.actionStatusList,muGraph,traversalGraph);
        [childrenList,parentType] = expandChildrenActionStatusList(node,assumptMuGraph,assumptTraversalGraph);
        for childNode = childrenList
            [assumptMuGraph,assumptTraversalGraph] = calculateAssumptGraphFromActionStatusList(childNode.state.actionStatusList,muGraph,traversalGraph); 
            childNode.f =  heuristicFunction(childNode.state.nodeId,assumptMuGraph,destinationNode);
            if childNode.state.nodeId == destinationNode
                childNode.status = true;
            else
                childNode.status = false;
            end
            AOtree = addChildNodeIntoAOTree(childNode,node.index,parentType,AOtree);
        end
        node = AOtree(node.index);
        [assumptMuGraph,assumptTraversalGraph] = calculateAssumptGraphFromActionStatusList(node.state.actionStatusList,muGraph,traversalGraph);
        AOtree = backpropagate(node,AOtree,omega,assumptMuGraph);    
    end
    path = outputPath(AOtree,root,destinationNode);
    if ~isempty(path)
        nextAction = path(2);
    else
        nextAtion = currenState;
    end


    function path1 = outputPath(AOtree,root,destinationNode)
        [~,treeSize] = size(root);
        path = [root.state.nodeId];
        path1 = [];
        nodeIndex = root.index;
        
        while AOtree(nodeIndex).state.nodeId ~= destinationNode
%             path = [path,AOtree(nodeIndex).state.nodeId];
            for childNodeIndex = AOtree(nodeIndex).children
                if AOtree(childNodeIndex).status == true
                    path = [path,AOtree(childNodeIndex).state.nodeId];
                    nodeIndex = childNodeIndex;
                    break;
                end
            end
        end
        path1 = [path1,path(1)];
        [~,pathSize] = size(path);
        for i = 2 : pathSize
            [~,m] = size(path1);
            if path(i) ~= path1(m)
                path1 = [path1,path(i)];
            end
        end
            
    end

    function AOtree = backpropagate(node,AOtree,omega,assumptMuGraph)
        if isempty(node.children)
            return;
        end
        nodeIndex = node.index;
        while ~isempty(nodeIndex)       % node ~= root
            if AOtree(nodeIndex).type == 0       % OR node
                minf = inf;
%                 minNodeIndex = nodeIndex;
                for childNodeIndex = AOtree(nodeIndex).children
                    AOtree(nodeIndex).status = AOtree(childNodeIndex).status | AOtree(nodeIndex).status;
                    if  AOtree(childNodeIndex).f < minf
                        minf = AOtree(childNodeIndex).f + assumptMuGraph(AOtree(nodeIndex).state.nodeId,AOtree(childNodeIndex).state.nodeId);
%                         minNodeIndex = childNodeIndex;
                    end
                end
                AOtree(nodeIndex).f = minf;
            elseif AOtree(nodeIndex).type == 1      %AND node
                valueExp = 0;
                status = true;
                for childNodeIndex = AOtree(nodeIndex).children
                    status = status & AOtree(childNodeIndex).status;
%                     valuef = valuef + AOtree(childNodeIndex).prob * (AOtree(childNodeIndex).f + assumptMuGraph(AOtree(nodeIndex).state.nodeId,AOtree(childNodeIndex).state.nodeId));
                    valueExp = valueExp + AOtree(childNodeIndex).prob * exp(omega * AOtree(childNodeIndex).f + assumptMuGraph(AOtree(nodeIndex).state.nodeId,AOtree(childNodeIndex).state.nodeId));
                end
                AOtree(nodeIndex).status = status;
                AOtree(nodeIndex).f = (1/omega) * log(valueExp);
            end
            nodeIndex = AOtree(nodeIndex).parent;
       
        end
    end

    function node = selectMinCostNode(AOtree)
        [~,treeSize] = size(AOtree);
        node = AOtree(1);
        if treeSize == 1
            return;
        end
        while ~isempty(node.children)
            minf = inf;
            nextNodeIndex = node.index;
            for childNodeIndex = node.children
                if AOtree(childNodeIndex).status == false
                    if AOtree(childNodeIndex).f < minf
                        minf = AOtree(childNodeIndex).f;
                        nextNodeIndex = AOtree(childNodeIndex).index;
                    end
                end
            end
            node = AOtree(nextNodeIndex);
        end
%         if node.status == true
%             nodeIndex = node.index;
%             while AOtree(nodeIndex).status == true
%                 solvedChildNodeIndex = nodeIndex;
%                 nodeIndex = AOtree(nodeIndex).parent;
%             end
%             minf = inf;
%             for childNodeIndex = AOtree(nodeIndex).children
%                 if childNodeIndex ~= solvedChildNodeIndex
%                     if AOtree(childNodeIndex).f < minf
%                         minf = AOtree(childNodeIndex).f;
%                         nextNodeIndex = AOtree(childNodeIndex).index;
%                     end
%                 end
%             end
%             node = AOtree(nextNodeIndex);
%         end
    end


    function AOtree = addChildNodeIntoAOTree(childNode,parentIndex,parentType,AOtree)
        [~,treeSize] = size(AOtree);
        childNode.index = treeSize + 1;
        childNode.children = [];
        childNode.type = [];
        childNode.nodeId = childNode.state.nodeId;
        AOtree = [AOtree,childNode];
        AOtree(parentIndex).children = [AOtree(parentIndex).children,childNode.index];
        AOtree(parentIndex).type = parentType;
    end

    function [childrenNodes,parentType] = expandChildrenActionStatusList(node,assumptMuGraph,assumptTraversalGraph)
        [~,nodeNum] = size(assumptTraversalGraph);
        parentActionStatusList = node.state.actionStatusList;
        childrenNodes = [];
        parentType = 0;
        for i = 1 : nodeNum
            if assumptTraversalGraph(node.state.nodeId,i) < 1 && assumptTraversalGraph(node.state.nodeId,i) > 0
                childrenNodes = [];
                childNode.parent = node.index;
                childNode.state.nodeId = node.state.nodeId;
                childNode.prob = 1 - assumptTraversalGraph(node.state.nodeId,i);
                [~,listLength] = size(parentActionStatusList);
                for j = 1 : listLength              %   先U后T
                    if parentActionStatusList(j).start == node.state.nodeId && parentActionStatusList(j).end == i
                        parentActionStatusList(j).status = 0;
                        childNode.state.actionStatusList = parentActionStatusList;
                        childrenNodes = [childrenNodes,childNode];
                        parentActionStatusList(j).status = 1;
                        childNode.state.actionStatusList = parentActionStatusList;
                        childNode.prob = assumptTraversalGraph(node.state.nodeId,i);
                        childrenNodes = [childrenNodes,childNode];
                        parentType = 1;         %   parent AND node
                        break;
                    end
                end
            end
        end
            if parentType == 0     
                for i = 1 : nodeNum
                    if assumptTraversalGraph(node.state.nodeId,i) == 1
                        childNode.parent = node.index;
                        childNode.state.nodeId = i;
                        childNode.prob = 1;
                        childNode.state.actionStatusList = parentActionStatusList;
                        childrenNodes = [childrenNodes,childNode];
                        parentType = 0;         %   parent OR node
                    end
                end
            end
        
                        
    end

    function actionStatusList = updateActionList(traversalGraph,assumptTraversalGraph)
        [~,nodeNum] = size(assumptTraversalGraph);
        actionStatusList = [];
        for i = 1 : nodeNum
            for j = 1 : nodeNum
                if traversalGraph(i,j) ~= 0
                    actionStatus.start = i;
                    actionStatus.end = j;
                    if assumptTraversalGraph(i,j) == 0
                        actionStatus.status = 0;
                    elseif assumptTraversalGraph(i,j) == 1
                        actionStatus.status = 1;
                    else
                        actionStatus.status = 2;
                    end
                    actionStatusList = [actionStatusList,actionStatus];
                end
            end
        end             
    end
    
    function [assumptMuGraph,assumptTraversalGraph] = calculateAssumptGraphFromActionStatusList(actionStatusList,muGraph,traversalGraph)
        [~,actionNum] = size(actionStatusList);
        assumptMuGraph = muGraph;
        assumptTraversalGraph = traversalGraph;
        for i = 1 : actionNum
            if actionStatusList(i).status == 0
                assumptMuGraph(actionStatusList(i).start,actionStatusList(i).end) = 0;
                assumptTraversalGraph(actionStatusList(i).start,actionStatusList(i).end) = 0;
            elseif actionStatusList(i).status == 1
                assumptTraversalGraph(actionStatusList(i).start,actionStatusList(i).end) = 1;
            end
        end
    end

    function heuristic = heuristicFunction(nodeId,assumptMuGraph,destinationNode)
        c = sparse(assumptMuGraph);
        [heuristic,~,~] = graphshortestpath(c,nodeId,destinationNode);
    end

end

