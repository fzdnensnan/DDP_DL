function nextAction = Zyf_RSP(muGraph,sigmaGraph,traversalGraph,currentState,destinationNode,synchronous)
persistent path_rsp_zyf;
persistent path;
% persistent A;
% persistent b;
% persistent mu;
% persistent Sigma;

if synchronous ~= 1
    nextAction = currentState;
    return;
end

if synchronous == 1
%     nextAction = getNextActionFromX(A,path_rsp_zyf,currentState);
%     path = convertX2Path(A,path_rsp_zyf,currentState,destinationNode);
    id = find(path == currentState);
    if traversalGraph(currentState,path(id + 1)) == 1
        nextAction = path(id + 1);
        return;
    end
end

[A,b,mu,Sigma] = convertGraph2Ab(muGraph,sigmaGraph,currentState,destinationNode);
%relative duality gap threshold epsilon,num is maximum number of iteration 
%% step 1,variance-covariance maxtrix decomposition
[~,num_edges]=size(A); %num_edges is the number of edges
[eig_Sigma,eig_Lambda] = eig(Sigma);% 对角化分解
%% step2 initialization
zeta = 1;
num = 100;
gap = 0.01;
iter_zyf=0;%set iteraton number
lower_bound_Z=-1000; %initial lower bounder of Z
Lagrangian_multiplier(:,iter_zyf+1)=zeros(1,num_edges);%initial lagrangian multiplier
x_let_dijkstra=func_dijkstra(A,b,mu);%a standard shorest path
x_let=zeros(num_edges,1);
length_x_let_dijkstra=length(x_let_dijkstra);
for i=1:length_x_let_dijkstra
    x_let(x_let_dijkstra(:,i))=1;
end

upper_bound_Z=mu'*x_let+zeta*sqrt(x_let'*Sigma*x_let);%compute Z as the upper bound

gap_zyf=(upper_bound_Z-lower_bound_Z)/upper_bound_Z;
while iter_zyf<num && (upper_bound_Z-lower_bound_Z)>gap
    %% step 3 solve lagtangian telaxation problems

    dynamic_mu=mu+eig_Sigma*Lagrangian_multiplier(:,iter_zyf+1);
    [x_let_dijkstra]=func_dijkstra(A,b,dynamic_mu);%a standard shorest path
    x_let=zeros(num_edges,1);
    length_x_let_dijkstra=length(x_let_dijkstra);
    for i=1:length_x_let_dijkstra
        x_let(x_let_dijkstra(:,i))=1;
    end
    optimal_upper_bound_Z=mu'*x_let+zeta*sqrt(x_let'*Sigma*x_let);% update upper bound of Z
    Lagrangian_function_optimal_value=(mu+eig_Sigma*Lagrangian_multiplier(:,iter_zyf+1))'*x_let;
    %% update upper_bound_Z and save x_let
    if upper_bound_Z>=optimal_upper_bound_Z
        upper_bound_Z=optimal_upper_bound_Z;
        x_zyf_total(:,iter_zyf+1)=x_let';
    else 
        if iter_zyf==0
            x_zyf_total(:,iter_zyf+1)=x_let';
        else
           x_zyf_total(:,iter_zyf+1)=x_zyf_total(:,iter_zyf);
        end
    end
    Lagrangian_function_value=Lagrangian_function_optimal_value;
    %% update low_bound_Z
    if Lagrangian_function_value>=lower_bound_Z
        lower_bound_Z=Lagrangian_function_value;
    end
    % update lagrangian multiplier
    partial_coefficient_lagrangian_multiplier=eig_Sigma'*x_zyf_total(:,iter_zyf+1);
    coefficient_lagrangian_multiplier=(upper_bound_Z-lower_bound_Z)/(partial_coefficient_lagrangian_multiplier'*partial_coefficient_lagrangian_multiplier);
    beta=0.1;
    Lagrangian_multiplier(:,iter_zyf+2)=Lagrangian_multiplier(:,iter_zyf+1)+beta*coefficient_lagrangian_multiplier*partial_coefficient_lagrangian_multiplier;
    for i=1:num_edges
        if Lagrangian_multiplier(i,iter_zyf+2)>zeta*sqrt(eig_Lambda(i,i))
            Lagrangian_multiplier(i,iter_zyf+2)=Lagrangian_multiplier(i,iter_zyf+1);
        elseif Lagrangian_multiplier(i,iter_zyf+2)<-zeta*sqrt(eig_Lambda(i,i))
            Lagrangian_multiplier(i,iter_zyf+2)=Lagrangian_multiplier(i,iter_zyf+1);
        end
    
    end
    iter_zyf=iter_zyf+1;
    path_rsp_zyf=x_zyf_total(:,end);
    gap_zyf=(upper_bound_Z-lower_bound_Z)/upper_bound_Z;
    upper_bound_total_zyf(iter_zyf)=upper_bound_Z;
    lower_bound_total_zyf(iter_zyf)=lower_bound_Z;
    gap_total_zyf(iter_zyf)=(upper_bound_Z-lower_bound_Z);
    gap_real_total_zyf(iter_zyf)=upper_bound_Z-lower_bound_Z;
    iteration_zyf(iter_zyf)=iter_zyf;
    
%     path = convertX2Path(A,path_rsp_zyf,currentState,destinationNode);

end
    path = convertX2Path(A,path_rsp_zyf,currentState,destinationNode);
    nextAction = getNextActionFromX(A,path_rsp_zyf,currentState);


    function nextAction = getNextActionFromX(A,X,currentState)
        [nodeNum,edgeNum] = size(A);
        nextAction = currentState;
        for i = 1 : edgeNum
            if X(i) == 1 && A(currentState,i) == 1
                for j = 1 : nodeNum
                    if A(j,i) == -1
                        nextAction = j;
                        break;
                    end
                end
                break;
            end
        end
    end

    function path = convertX2Path(A,X,currentState,destinationNode)
        [nodeNum,edgeNum] = size(A);
        nextAction = currentState;
        if sum(X) == 0
            path = [];
            return;
        end
        path = currentState;
        while nextAction ~= destinationNode
            for i = 1 : edgeNum
                if X(i) == 1 && A(nextAction,i) == 1
                    for j = 1 : nodeNum
                        if A(j,i) == -1
                            if ismember(j,path)
                                break;
                            end
                            nextAction = j;
                            path = [path,nextAction];
                            break;
                        end
                    end
                    break;
                end
            end
        end
    end

    function [x_let]=func_dijkstra(A,b,mu)
        %% calculate #nodes and #edges out of the map
        [~,num_edges]=size(A); % num_nodes is the number of nodes.num_edges is the number of edges

        %% construct the Dijkstra understandable graph
        start_node = zeros(num_edges,1);
        end_node = zeros(num_edges,1);

        for i=1:num_edges
            start_node(i)=find(A(:,i)==1,1); % find the first element whose value is 1
            end_node(i)=find(A(:,i)==-1,1); % find the first element whose value is -1
        end
        r=find(b==1); % r is origin
        s=find(b==-1);% s is the destinatioin
        start_node_real=start_node';
        end_node_real=end_node';
        mu_real=mu';
        graph_dijkstra = digraph(start_node_real,end_node_real,mu_real);

        %% use dijkstra's algorithm to find the shortest path
        [~,~,x_let] =shortestpath(graph_dijkstra,r,s);
    end
    function [A,b,mu,sigma] = convertGraph2Ab(muGraph,sigmaGraph,currentState,destinationNode)
        [~,nodeNum] = size(muGraph);
        edgeNum = 0;
        for val = muGraph(:)
            if val ~= 0
                edgeNum  = edgeNum + 1;
            end
        end
        A = zeros(nodeNum,edgeNum);
        b = zeros(nodeNum,1);
        mu = zeros(edgeNum,1);
        sigma = zeros(edgeNum,edgeNum);
        count = 1;
        for i = 1 : nodeNum
            for j = 1 : nodeNum
                if muGraph(i,j) ~= 0
                    A(i,count) = 1;
                    A(j,count) = -1;
                    mu(count,1) = muGraph(i,j);
                    sigma(count,count) = sigmaGraph(i,j);
                    count = count + 1;
                end
            end
        end
        b(currentState,1) = 1;
        b(destinationNode,1) = -1;
    end
end


