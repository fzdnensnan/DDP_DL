function nextAction = Kahni_RSP(muGraph,sigmaGraph,traversalGraph,currentState,destinationNode,synchronous)
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明
persistent outputPath;
% 
if synchronous ~= 1
    nextAction = currentState;
    return;
end

[~,n] = find(outputPath == currentState);
if traversalGraph(currentState,outputPath(n + 1)) == 1
    nextAction = outputPath(n + 1);
    return;
end
zeta = 10;
params = [];
count = 0;
maxIter = 100;
for i = 1 : 3
    param.alpha = -0.5 + i * 0.5;
    param.path = [];
    param.Z1 = 0;
    param.Z2 = 0;
    param.gama = 0;
    params = [params,param];
end
params(1).alpha = 0;
params(2).alpha = 0.5;
params(3).alpha = 1;
% params(4).alpha = (params(1).alpha + params(2).alpha) / 2;
% params(5).alpha = (params(2).alpha + params(3).alpha) / 2;
params = initZWithPath(params,1,muGraph,sigmaGraph,zeta,currentState,destinationNode);
params = initZWithPath(params,2,muGraph,sigmaGraph,zeta,currentState,destinationNode);
params = initZWithPath(params,3,muGraph,sigmaGraph,zeta,currentState,destinationNode);
while count < maxIter
    params(4).alpha = (params(1).alpha + params(2).alpha) / 2;
    params(5).alpha = (params(2).alpha + params(3).alpha) / 2;
    params = initZWithPath(params,4,muGraph,sigmaGraph,zeta,currentState,destinationNode);
    params = initZWithPath(params,5,muGraph,sigmaGraph,zeta,currentState,destinationNode);
    if params(2).Z1 < params(4).Z1 && params(2).Z1 < params(5).Z1
        params(1) = params(4);
        params(3) = params(5);
    else
        if params(4).Z1 < params(5).Z1 
            params(3) = params(2);
            params(2) = params(4);
        elseif params(4).Z1 > params(5).Z1
            params(1) = params(2);
            params(2) = params(5);
        else
            break;
        end
    end
    count = count + 1;
end
minId = 1;
minZ = Inf;
for i = 1 : 3
    if params(i).Z1 < minZ
        minId = i;
        minZ = params(i).Z1;
    end
end
outputPath = params(minId).path;
[~,lenth] = size(outputPath);
if lenth < 2
    nextAction = currentState;
else
    nextAction = outputPath(2);
end

    function params = initZWithPath(params,id,muGraph,sigmaGraph,zeta,currentState,destinationNode)
        params(id).gama = (1 - params(id).alpha) / params(id).alpha;
        c=sparse(params(id).alpha * muGraph + (1 - params(id).alpha) * (sigmaGraph.^2));
%         currentState
        [cost,path,~] = graphshortestpath(c,currentState,destinationNode);
        params(id).Z2 = cost / params(id).alpha;
        params(id).path = path;
        [~,pathLength] = size(path);
        if pathLength == 0
            cost
        end
        mu = 0;
        var = 0;
        for j = 1 : pathLength - 1
            mu = mu + muGraph(j,j+1);
            var = var + sigmaGraph(j,j+1)^2;
        end
        params(id).Z1 = mu + zeta * sqrt(var);

    end
    
end

