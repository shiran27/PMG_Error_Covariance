classdef Graph < handle
    % Class for generic graphs (also for sub-graphs, perhaps, in the
    % future)
    
    properties
        index
        targets = [] % list of target objects
        edges = []   % list of edge objects
        agents = [] % list of agents
        distanceMatrix %d_ij = direct edge distance
        angleMatrix %\theta_ij = edge angle
        shortestDistanceMatrix % based on distance
        shortestDistanceMatrix2 % based on the lower bound
        shortestPathMatrix = {}
        neighbors = {}
    end
    
    methods
        
        function obj = Graph(index)
            %GRAPH Construct an instance of this class
            %   Detailed explanation goes here
            obj.index = index;
        end
        
        
        function obj = loadARandomGraph(obj, numOfTargets, numOfAgents, dimentionOfSpace, sizeOfSpace, communicationRadius, targetControllersEnabled)
            
            % generate target List
            targetLocations = sizeOfSpace*rand(numOfTargets,dimentionOfSpace);
            
            % Target Parameters
            a = 0.01+0.4*rand(numOfTargets,1);
%             a = -0.01-0.4*rand(numOfTargets,1);
            q = 0.1+2*rand(numOfTargets,1);
            h = ones(numOfTargets,1);
            r = 10-8*rand(numOfTargets,1);
%             g = 1./r;
            
            obj.targets = [];
            for i = 1:1:numOfTargets
%                 newTarget = Target(i, targetLocations(i,:));
                newTarget = Target(i, a(i), q(i), h(i), r(i), targetLocations(i,:));
                obj.targets = [obj.targets, newTarget];
            end
            
        
            % generate edge List
            edgeCount = 1;
            connected = zeros(1,numOfTargets);
            obj.edges = [];
            for i = 1:1:(numOfTargets-1)
                for j = (i+1):1:numOfTargets
                    newEdge = Edge(edgeCount,[obj.targets(i).index, obj.targets(j).index],[obj.targets(i).position; obj.targets(j).position]);
                    if norm(obj.targets(i).position - obj.targets(j).position, 2) < communicationRadius
                        newEdge.enabled = true;
                        connected(i) = 1; 
                        connected(j) = 1;
                    else
                        newEdge.enabled = false;
                    end
                    obj.edges = [obj.edges, newEdge];
                    edgeCount = edgeCount + 1;
                end
            end
            
            if sum(connected==0)>0 % not fully connected!
                disp('Not fully connected, use multiple agents!');
                % obj = loadARandomGraph(obj, numOfTargets, dimentionOfSpace, sizeOfSpace, communicationRadius);
            end
            
            % omitted for now
            obj = obj.loadDistanceMatrix();
% %             obj = obj.loadShortestPaths();
            % end
            obj.loadNeighbors();
            
            % loading agents
            agentSpeed = 1;
            obj.agents = [];
            for i = 1:1:numOfAgents
                tarId = round(numOfTargets/numOfAgents)*(i-1)+1;
                obj.targets(tarId).residingAgents = [i];
                agentPosition = obj.targets(tarId).position;
                newAgent = Agent(i,agentPosition,agentSpeed,tarId,numOfTargets);
                obj.agents = [obj.agents, newAgent];
            end
            
            % loading initial target state
            for i = 1:1:numOfTargets
                obj.targets(i).Omega = obj.targets(i).Q; % initial error covariance
                obj.targets(i).phi = 0.5; % initial target state
                obj.targets(i).phiHat = 0.5; % estimate of the initial target state
                obj.targets(i).timeSinceLastEvent = 0; %time counter
                obj.targets(i).OmegaAtLastEvent = obj.targets(i).Omega;
            end
            
            
            % target control related information
            b = 0.01+0.4*rand(numOfTargets,1);
% %             r = 5+round(10*rand(numOfTargets,1)); 
% %             alpha = 1; % error dynamics \dot{e} = -alpha e
            for i = 1:1:numOfTargets
                obj.targets(i).B = b(i); % input matrix
                obj.targets(i).K = 2;%(alpha + a(i))/b(i); % initial target state
                obj.targets(i).r = 10*sin(i); % reference input
                if targetControllersEnabled
                    obj.targets(i).controllerEnabled = true;
% %                     obj.targets(i).u = -obj.targets(i).K*obj.targets(i).phiHat + (obj.targets(i).K - a(i)/b(i))*r(i); % initial input based on state estimate
                    rDot = 20*cos(i);
                    obj.targets(i).u = -(1/obj.targets(i).B)*((obj.targets(i).K+obj.targets(i).A)*obj.targets(i).phiHat - (rDot+obj.targets(i).K*obj.targets(i).r));
                else
                    obj.targets(i).controllerEnabled = false;
                    obj.targets(i).u = 0;
                end
            end
            
            
        end
        
        function result = getTarget(obj, targetIndex, requirement)
            count1 = 1;
            for target = obj.targets % search throguh all the targets
                if target.index == targetIndex 
                   if isequal(requirement,'id')
                       result = count1; break;
                   elseif isequal(requirement,'target')
                       result = target; break;
                   end
                end
                count1 = count1 + 1;
            end
        end
        
        function outputArg = getEdgeIndex(obj,targetPair)
            
              % smart way : only works for the main graph
              if obj.index == 0
                i = min(targetPair);
                j = max(targetPair);
                N = length(obj.targets);
                ind = N*(N-1)/2 -((N-i)*(N-i-1)/2 + N - j);
                outputArg = ind; 
                return;
              end
              
              % old way - if this is a subgraph!
            reversedTargetPair = [targetPair(2), targetPair(1)];
            for i = 1:1:length(obj.edges) % search throguh all the edges
                if isequal(obj.edges(i).targets, targetPair) || isequal(obj.edges(i).targets, reversedTargetPair)
                    outputArg = i; % obj.edges(i).index; 
                    return;
                end
            end
            
            outputArg = -1;
            return % not a path on the sub-graph
        end
        
        
        function outputArg = drawGraph(obj, figNum)
            figure(figNum)
            
            grid on
            axis equal
            sizeOfSpace = 1;
            axis([0,sizeOfSpace,0,sizeOfSpace]);
            
            % draw edges
            for i = 1:1:length(obj.edges)
                obj.edges(i).drawEdge();
            end
            
            % draw targets
            for i = 1:1:length(obj.targets)
                obj.targets(i).drawTarget();
            end
            
            xlabel('X')
            ylabel('Y')
%             title('Mission Space')
            axis equal
            axis([0,sizeOfSpace,0,sizeOfSpace]);
        end
        
        
        function outputArg = updateGraph(obj, figNum, timeVal)
            figure(figNum)
            
            % update targets
            thickness = 0.01; 
            scalingFactor = 500; %1000
            for i = 1:1:length(obj.targets)
                posX = obj.targets(i).position(1);
                posY = obj.targets(i).position(2);
                normOmega_i = trace(obj.targets(i).Omega);
                if timeVal ~= 0
                    delete(obj.targets(i).graphicHandles(1));
                    delete(obj.targets(i).graphicHandles(2));
                end
                height1 = normOmega_i/scalingFactor;
                r = rectangle('Position',[posX-thickness, posY, 2*thickness, height1],...
                    'FaceColor',[0.8 0.8 0 0.5],'EdgeColor','k','LineWidth',0.01);
                t = text(posX+0.02, posY+0.02,num2str(normOmega_i,4),'Color','b','FontSize',10);
                obj.targets(i).graphicHandles(1) = r;
                obj.targets(i).graphicHandles(2) = t;
                
                if obj.targets(i).controllerEnabled
                    if timeVal~=0
                        delete(obj.targets(i).graphicHandles(3));
                    end
                    scalingFactor2 = 1;
                    height2 = scalingFactor2*abs(obj.targets(i).phi - obj.targets(i).r); % phiHat-r
                    r2 = rectangle('Position',[posX-thickness, posY+height1, 2*thickness, height2],...
                    'FaceColor',[0.8 0.1 0 0.1],'EdgeColor','k','LineWidth',0.01);
                    obj.targets(i).graphicHandles(3) = r2; 
                end
            end
            
            
            % update agents
            thickness = 0.04; 
            triangle = thickness*[-0.5 -0.5 1; 0.5 -0.5 0];
            for i = 1:1:length(obj.agents)
                pos = transpose(obj.agents(i).position);
                theta = obj.agents(i).orientation;
                rotMatrix = [cos(theta), -sin(theta); sin(theta) cos(theta)];
                triangleT = rotMatrix*triangle + pos;
                if timeVal ~= 0
                    delete(obj.agents(i).graphicHandles(1));
                    delete(obj.agents(i).graphicHandles(2));
                end
                pgon = polyshape(triangleT(1,:),triangleT(2,:));
                r = plot(pgon,'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor','k','LineWidth',0.01);
                t = text(pos(1)-0.04,pos(2)+0.02,num2str(i),'Color','r','FontSize',10);
                obj.agents(i).graphicHandles(1) = r;
                obj.agents(i).graphicHandles(2) = t;
            end
%             pgon = polyshape([0 0 1.5],[0.5 -0.5 0])
         end
        
        function output = triggerCoverednessEvents(obj, startNode, endNode)
            for targetId = 1:1:length(obj.targets)
                if obj.distanceMatrix(startNode,targetId) < 1.5 &  startNode ~= targetId & ~isempty(obj.targets(targetId).residingAgents)
                    agentId = obj.targets(targetId).residingAgents(1);
                    obj.agents(agentId).coverednessEventTriggered = true;
                elseif obj.distanceMatrix(endNode,targetId) < 1.5 &  endNode ~= targetId & ~isempty(obj.targets(targetId).residingAgents)
                    agentId = obj.targets(targetId).residingAgents(1);
                    obj.agents(agentId).coverednessEventTriggered = true;
                end
            end
        end
         
         
        function output = triggerExternalStatePerturbation(obj,activate)
        
            magnitudeQ = 2; % increase the process noise covariance at each target : Q
            magnitudeR = 1; % increase the measurement noise covariance at each agent : R
            for i = 1:1:length(obj.targets)
                if activate
                    obj.targets(i).Q = obj.targets(i).Q*magnitudeQ;
                    obj.targets(i).R = obj.targets(i).R*magnitudeR;
                else
                    obj.targets(i).Q = obj.targets(i).Q/magnitudeQ;
                    obj.targets(i).R = obj.targets(i).R/magnitudeR;
                end
            end
            
            % trigger coverdness event at all the agents!
            for a = 1:1:length(obj.agents)
                obj.agents(a).coverednessEventTriggered = true;
            end
            
        end
        
        
        function outputArg = drawAll(~, subGraphs, cycles)
            figNum = get(gcf,'Number')+1;
            
            for i = 1:1:length(cycles)
                subGraphs(i).drawGraph(figNum); 
                cycles(i).drawCycle(subGraphs(i),'y'); 
            end
            
            
%             for i = 1:1:length(cycles)
%                 cycles(i).drawFullCycle(subGraphs(i)); 
%             end
        end
            
            
            
            
        function obj = loadDistanceMatrix(obj) 
            numOfTargets = length(obj.targets);
            
            dmat = ones(numOfTargets)-eye(numOfTargets);
            dmat(dmat==1)=1000;
            amat = zeros(numOfTargets);
            for i = 1:1:length(obj.edges)
                iInd = obj.getTarget(obj.edges(i).targets(1),'id');
                jInd = obj.getTarget(obj.edges(i).targets(2),'id');
                distVal = obj.edges(i).getLength();
                dmat(iInd,jInd) = distVal;
                dmat(jInd,iInd) = distVal;
                
                posI = obj.getTarget(obj.edges(i).targets(1),'target').position;
                posJ = obj.getTarget(obj.edges(i).targets(2),'target').position;
                deltaPos = posJ-posI;
                angleVal = atan2(deltaPos(2),deltaPos(1));
                amat(iInd,jInd) = angleVal;
                angleVal = atan2(-deltaPos(2),-deltaPos(1));
                amat(jInd,iInd) = angleVal;
            end
            obj.distanceMatrix = dmat;
            obj.angleMatrix = amat;
        end
        
        function TSPCycleFound = findTSPCycle(obj)
            
            numOfTargets = length(obj.targets);
            xy = [];
            for i = 1:1:numOfTargets
                xy = [xy; obj.targets(i).position];
            end
            
            dmat = obj.distanceMatrix;
            
            userConfig = struct('xy',xy,'dmat',dmat,'showProg',false,'showResult',false,'showWaitbar',false);
            resultStruct = tsp_ga(userConfig);
            TSPCycleFound = resultStruct.optRoute;
            
            count = 1;
            for i = TSPCycleFound
                TSPCycleFound(count) = obj.targets(i).index;
                count = count + 1;
            end
            
        end
        
        
        function subGraphs = generateSubGraphs(obj, partitions)
            
            subGraphs = []; count = 1;
            for targetIDList = partitions
                
                subGraph = Graph(count); % create the sub-graph object
                
                subGraph.targets = []; % loading targets
                targetIDList = cell2mat(targetIDList);
                for targetID = targetIDList
                    subGraph.targets = [subGraph.targets, obj.targets(targetID)];
                end
                
                numOfTargets = length(targetIDList); % loading edges
                connected = zeros(1,numOfTargets);
                subGraph.edges = [];
                for i = 1:1:(numOfTargets-1)
                    for j = (i+1):1:numOfTargets
                        edgeID = obj.getEdgeIndex([targetIDList(i),targetIDList(j)]);
                        subGraph.edges = [subGraph.edges, obj.edges(edgeID)];
                        if obj.edges(edgeID).enabled
                            connected(i) = 1; 
                            connected(j) = 1;
                        end
                    end
                end

                if sum(connected==0)>0 % not fully connected!
                    disp(['Sub graph ',num2str(count),' is not fully connected.']);
                end
            
                subGraph = subGraph.loadDistanceMatrix();
                subGraph = subGraph.loadShortestPaths();
                
                subGraphs = [subGraphs, subGraph]; % lead sub-graph
                count = count + 1;
                
            end
             
        end
        
        function obj = loadShortestPaths(obj)
            
            N = length(obj.targets);
            countingVector = 1:1:N;
            
            for startNode = 1:1:N
                
                verticesToCover = ones(1,N); % vertices to cover
                distances = inf*ones(1,N); % shortest distances found so far
                priors = zeros(1,N); % previous target to visit

                distances(startNode) = 0;                        

                while sum(verticesToCover)>0
                    
                    [~,uInd] = min(distances(verticesToCover==1)); % u = vertex in Q with min dist[u]
                    counts = countingVector(verticesToCover==1);
                    u = counts(uInd);
                    
                    verticesToCover(u) = 0; % remove u from Q

                    for v = 1:1:N
                        if verticesToCover(v)>0 % for each neighbor v of u that are still in Q
                            alt = distances(u) + obj.distanceMatrix(u,v); % alt ? dist[u] + length(u, v)
                            if alt < distances(v)               
                                distances(v) = alt;
                                priors(v) = u;
                            end
                        end
                    end
                end

                % unvinding paths
                shortestPaths = cell(1,N);
                for destination = 1:1:N % for different destinations 
                    u = destination;
                    uReal = obj.targets(u).index;
                    path = [];
                    while priors(u)>0
                        path = [uReal ,path];
                        u = priors(u);
                        uReal = obj.targets(u).index;
                    end
                    shortestPaths{destination} = path;
                end  
                
                % loading data
                obj.shortestDistanceMatrix(startNode,:) = distances;
                obj.shortestPathMatrix(startNode,:) = shortestPaths;
            
            end
            
        end
        
        
        
        function path = getHandicappedShortestPath(obj,startNode,endNode)
            
            startNode2 = obj.getTarget(startNode,'id');
            endNode2 = obj.getTarget(endNode,'id');
            
            N = length(obj.targets);
            countingVector = 1:1:N;
            
            verticesToCover = ones(1,N); % vertices to cover
            distances = inf*ones(1,N); % shortest distances found so far
            priors = zeros(1,N); % previous target to visit

            distances(startNode2) = 0;                        

            while sum(verticesToCover)>0

                [~,uInd] = min(distances(verticesToCover==1)); % u = vertex in Q with min dist[u]
                counts = countingVector(verticesToCover==1);
                u = counts(uInd);
                if u == endNode2
                    break; % arrived at the destination
                else
                    verticesToCover(u) = 0; % remove u from Q, and continue the search
                end
                
                
                for v = 1:1:N % update shortest distances to the neighbors
                    if verticesToCover(v)>0 % for each neighbor v of u that are still in Q
                        
                        if u == startNode2 && v == endNode2 % directlyConnectedScenario!
                            distuv = 1000;
                        else
                            distuv = obj.distanceMatrix(u,v);
                        end
                        
                        alt = distances(u) + distuv; % alt ? dist[u] + length(u, v)
                        if alt < distances(v)               
                            distances(v) = alt;
                            priors(v) = u;
                        end
                    end
                end
                
            end

            % unvinding the path
            u = endNode2;
            uReal = obj.targets(u).index;
            path = [];
            while priors(u)>0
                path = [uReal ,path];
                u = priors(u);
                uReal = obj.targets(u).index;
            end
            
        end
        
        
        function path = getHandicappedShortestPath2(obj,startNode,endNode,avoidNode)
            
            startNode2 = obj.getTarget(startNode,'id');
            endNode2 = obj.getTarget(endNode,'id');
            avoidNode2 = obj.getTarget(avoidNode,'id');
            
            N = length(obj.targets);
            countingVector = 1:1:N;
            
            verticesToCover = ones(1,N); % vertices to cover
            distances = inf*ones(1,N); % shortest distances found so far
            priors = zeros(1,N); % previous target to visit

            distances(startNode2) = 0;                        

            while sum(verticesToCover)>0

                [~,uInd] = min(distances(verticesToCover==1)); % u = vertex in Q with min dist[u]
                counts = countingVector(verticesToCover==1);
                u = counts(uInd);
                if u == endNode2
                    break; % arrived at the destination
                else
                    verticesToCover(u) = 0; % remove u from Q, and continue the search
                end
                
                
                for v = 1:1:N % update shortest distances to the neighbors
                    if verticesToCover(v)>0 % for each neighbor v of u that are still in Q
                        
                        if u == avoidNode2 || v == avoidNode2 % arriving or going to go from avoided node!
                            distuv = 1000;
                        else
                            distuv = obj.distanceMatrix(u,v);
                        end
                        
                        alt = distances(u) + distuv; % alt ? dist[u] + length(u, v)
                        if alt < distances(v)               
                            distances(v) = alt;
                            priors(v) = u;
                        end
                    end
                end
                
            end

            % unvinding the path
            u = endNode2;
            uReal = obj.targets(u).index;
            path = [];
            while priors(u)>0
                path = [uReal ,path];
                u = priors(u);
                uReal = obj.targets(u).index;
            end
            
        end
        
        function disparityMatrix = getDisparityMatrix(obj, method)
            
            if isequal(method,'shortest path')
                disparityMatrix =  obj.shortestDistanceMatrix;
                return
            else
                disparityMatrix =  obj.shortestDistanceMatrix2;
            end
            
            % The long way! - disparity based on intermediate target properties as well
            % as edge properties
            
        end
        
        function partitions = partitionTheGraph(obj, numOfAgents, method)
%             partitions = {[1,8,10,13,15],[3,4,5,7,9],[2,6,11,12,14]}; 
            disparityMatrix = obj.getDisparityMatrix(method); % disparity Matrix
            
            if isequal(method,'shortest path')
                sigma = 1/sqrt(numOfAgents);
            else
                d_max = max(max(disparityMatrix));
                d_min = min(min(disparityMatrix(disparityMatrix>0)));
                if d_max > 100000000
                    sigma = (3*d_min);
                else
                    sigma = (0.8*d_min+0.2*d_max); % neighborhood width
                end
            end
            disparityMatrix;
            similarityMatrix = exp(-disparityMatrix.^2/(2*sigma^2)) % similarity matrix
            
            W = similarityMatrix; % weighted adjacency matrix
            
            d_i = sum(W,2);
            D = diag(d_i); % degree matrix
            
            L = D - W; % unnormalized laplacian
            L_rw = D\L; % normalized laplacian based on random walks
            
            [eigVec,eigVal] = eig(L_rw); % need to pick only the first non zero minimum amplitude eigenvalues!
            [sortedEigs,sortedIndices] = sort(diag(eigVal)); % minimum eigenvalue should be zero
            if eigVal(sortedIndices,sortedIndices)~=0
                disp('Error');
            end
            requiredInds = sortedIndices(1:(numOfAgents)); % (2:(numOfAgents+1)) %
            U = eigVec(:,requiredInds'); % soretd eigenvectors: \in R^{M \times N} % N=numOfagents, M=numOfTargets
            
            % M - rows of U needs to be clustered using k-means 
            idx = kmeans(U,numOfAgents); % \in \R^{M \times 1} - cluster index of each element
            
            partitions = cell(1,numOfAgents);
            targetPool = transpose(1:1:length(obj.targets));
            for i = 1:1:numOfAgents
                partition = targetPool(idx==i);
                partitions{1,i} = partition';
            end
            
        end
        
        function [minDistances, minCycles, currentCycle] = constructACycleStartingFromANode(obj,startNode, method)
            currentCycle = Cycle([startNode]);
            currentCycle = currentCycle.updateEdgeList(obj);
            minDistances = inf*ones(1,length(obj.targets));
            minCycles = repmat(currentCycle,1,length(obj.targets));
            
            i = obj.getTarget(startNode,'id');
            minDistances(i) = currentCycle.getLowerBound(obj);
            minCycles(i) = currentCycle;
            if isequal(method,'slow')
                numOfExternalTargets = length(obj.targets)-1;
            else
%                 numOfExternalTargets = round(length(obj.targets)); % For Case 7
                numOfExternalTargets = round(0.7*length(obj.targets));
            end
            
            for count = 1:1:numOfExternalTargets
            
                minlb = inf;
                minExpandedCycle = currentCycle;
                minExternalTarget = startNode;
                for target = obj.targets
                    if sum(target.index == currentCycle.targetList)>0
                        continue % already in the current cycle
                    end
                    expandedCycle = currentCycle.expandCycleToInclude(obj,target.index, method);
                    [lb,~] = expandedCycle.getLowerBound(obj);
                    if lb < minlb
                        minlb = lb;
                        minExpandedCycle = expandedCycle;
                        minExternalTarget = target.index;
                    end
                end

                currentCycle = minExpandedCycle;
                
                % Case 7
%                 obj.drawGraph(count+1);
%                 currentCycle.drawCycle(obj,'r');
%                 count
                % End Case 7
                
%                 currentCycle.drawFullCycle(obj);
                i = obj.getTarget(minExternalTarget,'id');
                minDistances(i) = minlb;
                minCycles(i) = currentCycle;
                
                
%                 minExternalTarget
%                 currentCycle.targetList
            
            end

        end
        
        function obj = loadShortestPaths2(obj,method)
            
            N = length(obj.targets);
            countingVector = 1:1:N;
            
            for startNode = 1:1:N
                disp(['Similarities w.r.t. Node ',num2str(startNode),' loaded.'])
                
                [dis,~,~] = obj.constructACycleStartingFromANode(startNode,method) 
%                 obj.drawGraph(startNode+2);
%                 fCyc.drawCycle(obj,'y');
%                 fCyc.drawFullCycle(obj);
                % loading data
                obj.shortestDistanceMatrix2(startNode,:) = dis;
            
            end
            obj.shortestDistanceMatrix2 = min(obj.shortestDistanceMatrix2,obj.shortestDistanceMatrix2')
            
        end
        
        
        
        
        function [subGraphs,cycles] = balanceThePartitions(obj,subGraphs,cycles) % Not used!
            % find the sub-graph, cycle and the set of targets that can be sold
            
            % finding the worst target in the worst cycle/subGraph
            maxlb = 0;
            maxSubGraph = 0; % dummy
            maxTarget = 0; % dummy
            subGraphIndex = 1;
            for cycle = cycles
                [lb,ij] =  cycle.getLowerBound(obj);
                if lb > maxlb
                    maxlb = lb;
                    maxSubGraph = subGraphIndex; 
                    maxTarget = ij(1);
                end
                subGraphIndex = subGraphIndex + 1;
            end
            
            % finding the gain in the global cost if the max target was
            % covered by some other
            
            
            
            % Remobing the max target: COntraction
            contractedTargetList = [];
            for target = subGraphs(maxSubGraph).targets
                if target.index~= maxTarget
                    contractedTargetList = [contractedTargetList, target.index];
                end
            end
            contractedSubGraph = obj.generateSubGraphs({contractedTargetList});
            contractedSubGraph = contractedSubGraph(1);
            cycleTSP = Cycle(contractedSubGraph.findTSPCycle()); % Cycle([8,1,10,8,1,10,8,15,13]);
            cycleTSP = cycleTSP.updateEdgeList(contractedSubGraph);
            contractedCycle = cycleTSP.executeGreedyExpansions(contractedSubGraph);
            [lbCon,~] = contractedCycle.getLowerBound(contractedSubGraph);
            
            bestExpandedCycleIndex = 0;
            bestExpandedCycle = [];
            bestExpandedSubGraph = [];
            maxGain = 0;
            for i = 1:1:length(cycles)
                if i == maxSubGraph
                    continue;
                end
                
                expandedTargetList = [maxTarget];
                for target = subGraphs(i).targets
                    expandedTargetList = [expandedTargetList, target.index];
                end
                
                expandedSubGraph = obj.generateSubGraphs({expandedTargetList});
                expandedSubGraph = expandedSubGraph(1);
                cycleTSP = Cycle(expandedSubGraph.findTSPCycle()); % Cycle([8,1,10,8,1,10,8,15,13]);
                cycleTSP = cycleTSP.updateEdgeList(expandedSubGraph);
                expandedCycle = cycleTSP.executeGreedyExpansions(expandedSubGraph);
                [lbExp,~] = expandedCycle.getLowerBound(expandedSubGraph);
               
                gain = maxlb - max(lbCon,lbExp);
                if gain > maxGain
                    disp('Interation Helped!')
                    maxGain = gain
                    bestExpandedCycleIndex = i;
                    bestExpandedSubGraph = expandedSubGraph;
                    bestExpandedCycle = expandedCycle;
                    % Great 1
                end
            end
            
            % answer
            if bestExpandedCycleIndex~=0
                subGraphs(maxSubGraph) = contractedSubGraph;
                subGraphs(bestExpandedCycleIndex) = expandedSubGraph;
                cycles(maxSubGraph) = contractedCycle;
                cycles(bestExpandedCycleIndex) = expandedCycle;
            end
            
        end
        
        
        % method 2
        function [subGraphs,cycles,gainVal] = balanceThePartitions2(obj,subGraphs,cycles)
            gainVal = 0;
            % try to sell one of the targets in the critical subcycle
            
            % finding the worst target in the worst cycle/subGraph
            maxlb = 0;
            maxInd = 0; % dummy
            maxTarget = 0; % dummy
            subGraphIndex = 1;
            for cycle = cycles
                [lb,ij] =  cycle.getLowerBound(obj);
                if lb > maxlb
                    maxlb = lb;
                    maxInd = subGraphIndex; 
                    maxTarget = ij;
                end
                subGraphIndex = subGraphIndex + 1;
            end
            
            % finding indivdual targets to sell
           indexInCycle = maxTarget(3);
           
           subCycleVector = cycles(maxInd).subCycleMatrix(:,indexInCycle); % we need find the individual targets in this subcycle
           subCycleIndices = find(subCycleVector==1); 
           targetsToRemove = [];
           
           for ind = subCycleIndices'
               target = cycles(maxInd).targetList(ind);
               if sum(targetsToRemove==target)==0
                    targetsToRemove = [targetsToRemove, target];
               end
           end
           targetsToRemove;
            
           
           % finding the target TO remove and the neighbor to expand
           bestTargetToSell = 0;
           bestNeighborToBuy = 0;
           maxGain = 0;
           
           bestConExpCycles = [];
           bestlbs = [];
           for target = targetsToRemove
               
                % finding the gain in the global cost if the target was
                % covered by some other
                
                % Remobing the max target: COntraction
%                 cycles(maxInd).targetList
%                 target
                contractedCycle = cycles(maxInd).contractCycleToExclude(subGraphs(maxInd),target);
                [lbCon,~] = contractedCycle.getLowerBound(subGraphs(maxInd));
                
                for i = 1:1:length(cycles)
                    if i == maxInd
                        continue;
                    end
                    
                    % need to expand
                    expandedCycle = cycles(i).expandCycleToInclude(obj,target,'slow');
                    [lbExp,~] = expandedCycle.getLowerBound(obj);

                    gain = maxlb - max(lbCon,lbExp);
                    if gain > maxGain
%                         disp('Interation Helped!')
                        maxGain = gain;
                        bestTargetToSell = target;
                        bestNeighborToBuy = i;
                        contractedCycle;
                        bestConExpCycles = [contractedCycle, expandedCycle];
                        bestlbs = [lbCon, lbExp];
                    end
                end
           end
           
            % final answer
            if bestNeighborToBuy~=0
                
                % contracted
                contractedTargetList = [];
                for target = subGraphs(maxInd).targets
                    if target.index ~= bestTargetToSell
                        contractedTargetList = [contractedTargetList, target.index];
                    end
                end
                contractedSubGraph = obj.generateSubGraphs({contractedTargetList});
                contractedSubGraph = contractedSubGraph(1);
                cycleTSP = Cycle(contractedSubGraph.findTSPCycle()); 
                cycleTSP = cycleTSP.updateEdgeList(contractedSubGraph);
                cycleTSP = cycleTSP.createCycleEvaluator(contractedSubGraph);
                cycleTSP.searchWithLB = 1;
                contractedCycle = cycleTSP.executeGreedyExpansions(contractedSubGraph);
                [lbCon,~] = contractedCycle.getLowerBound(subGraphs(maxInd));
                
                                
                % Expanded
                expandedTargetList = [bestTargetToSell];
                for target = subGraphs(bestNeighborToBuy).targets
                    expandedTargetList = [expandedTargetList, target.index];
                end
                expandedSubGraph = obj.generateSubGraphs({expandedTargetList});
                expandedSubGraph = expandedSubGraph(1);
                cycleTSP = Cycle(expandedSubGraph.findTSPCycle()); % Cycle([8,1,10,8,1,10,8,15,13]);
                cycleTSP = cycleTSP.updateEdgeList(expandedSubGraph);
                cycleTSP = cycleTSP.createCycleEvaluator(expandedSubGraph);
                cycleTSP.searchWithLB = 1;
                expandedCycle = cycleTSP.executeGreedyExpansions(expandedSubGraph);
                [lbExp,~] = expandedCycle.getLowerBound(obj);
                
                % gain after obtainng and refining TSP solutions may be
                % different
                gain = maxlb - max(lbCon,lbExp); 
                if max(bestlbs) < max(lbCon,lbExp) % Regular contracted and expanded cycles are better than TSP
                    disp('TSP solution improvement is worse!')
%                     bestConExpCycles(1,1)
%                     contractedSubGraph
                    subGraphs(maxInd) = contractedSubGraph;
                    cycles(maxInd) = bestConExpCycles(1,1).updateEdgeList(contractedSubGraph);
                    subGraphs(bestNeighborToBuy) = expandedSubGraph;
                    cycles(bestNeighborToBuy) = bestConExpCycles(1,2).updateEdgeList(expandedSubGraph);
                
                    targetSold = bestTargetToSell
                    gainVal = maxGain
                elseif gain>0   % TSP cycles are better
                    disp('TSP solution improvement is better!')
                    subGraphs(maxInd) = contractedSubGraph;
                    cycles(maxInd) = contractedCycle;
                    subGraphs(bestNeighborToBuy) = expandedSubGraph;
                    cycles(bestNeighborToBuy) = expandedCycle;
                
                    targetSold = bestTargetToSell
                    gainVal = maxGain
                else
                    disp('TSP solution improvement is neagative!')
                end
                
            end
            
            
        end
        
        function output = assignCyclesToAgents(obj, controlMethod)
            
            numOfAgents = length(obj.agents);
            
            obj.loadShortestPaths2('fast'); % try 'slow' if have patience 
            partitions = obj.partitionTheGraph(numOfAgents,'smallest lower bound'); % spectral clustering approach 2

            subGraphs = obj.generateSubGraphs(partitions);
            cycles = [];
            for subGraph = subGraphs

                cycleTSP = Cycle(subGraph.findTSPCycle()); % Cycle([8,1,10,8,1,10,8,15,13]);
                cycleTSP = cycleTSP.updateEdgeList(subGraph);
                cycleTSP = cycleTSP.createCycleEvaluator(subGraph);
                cycleTSP.searchWithLB = 1;
                
                greedyExpandedCycle = cycleTSP.executeGreedyExpansions(subGraph);
                cycles = [cycles, greedyExpandedCycle];  

            end
            %drawing the initial subgraphs and cycles
            obj.drawAll(subGraphs,cycles);
            for i = 1:1:length(cycles)
                cycles(i).drawFullCycle(subGraphs(i)); 
            end
            % % end

            % % interchange algorithm
            [subGraphs,cycles,gainVal] = obj.balanceThePartitions2(subGraphs,cycles);
            while gainVal > 0
                obj.drawAll(subGraphs,cycles); % drawing intermediate subgraphs and cycles
                [subGraphs,cycles,gainVal] = obj.balanceThePartitions2(subGraphs,cycles);
            end
            % drawing the final cycles
            for i = 1:1:length(cycles)
                if isequal(controlMethod,'Periodic')
                    cost = cycles(i).getCost()
                end
                cycles(i).drawFullCycle(subGraphs(i)); 
            end
            
            
            % for each agent, we need to find the closest cycle and the
            % entry path
            obj.loadShortestPaths(); % load dijkstras shortest paths and distances
            assignmentCosts = zeros(numOfAgents,numOfAgents);
            assignmentTarget = zeros(numOfAgents,numOfAgents); % entry point target
            iCount = 1;
            for agent = obj.agents
                startTarget = agent.residingTarget;
                jCount = 1;
                for cycle = cycles
                    shortestDist = inf;
                    shortestDistTarget = 0;
                    for endTarget = cycle.targetList
                        distVal = obj.shortestDistanceMatrix(startTarget,endTarget);
                        if distVal<shortestDist
                            shortestDist = distVal;
                            shortestDistTarget = endTarget;
                        end
                    end
                    assignmentCosts(iCount,jCount) = shortestDist; % best cost
                    assignmentTarget(iCount,jCount) = shortestDistTarget; % best entry point target
                    jCount = jCount + 1;
                end
                iCount = iCount + 1;
            end
            % solving the assignment problem
            matchings = matchpairs(assignmentCosts,1000000);% list of [agent,cycle] pairs
            
            for match = 1:1:size(matchings,1)
                agentId = matchings(match,1);
                cycleId = matchings(match,2);
                startTarget = obj.agents(agentId).residingTarget;
                endTarget = assignmentTarget(agentId,cycleId);
                shortestPath = obj.shortestPathMatrix(startTarget,endTarget);
                obj.agents(agentId).assignedCycle = cycles(cycleId);
                obj.agents(agentId).entryPathToCycle = [startTarget, shortestPath{1}];
                
                entryPathState = -length(shortestPath{1});
                cycleState = find(cycles(cycleId).targetList==endTarget);
                cycleState = cycleState(1);
                obj.agents(agentId).cycleState = [entryPathState, cycleState];
                
                % dwell times
                if isequal(controlMethod,'Periodic')
%                     cost = cycle(cycleId).getCost(); % prerequisite for dwell time computation
                    obj.agents(agentId).dwellTimes = cycles(cycleId).cycleEvaluator.getDwellTimes();
                end
            end
            % entry paths and cycles are loaded to each agent!
            
            
        end
        
        
        function output = loadNeighbors(obj)
            
            for currentTargetId = 1:1:length(obj.targets) 
                neighborSet = [currentTargetId];
                for targetId = 1:1:length(obj.targets)
                    if obj.distanceMatrix(currentTargetId,targetId) < 1.5 & currentTargetId ~= targetId 
                        neighborSet = [neighborSet, targetId];
                    end
                end
                obj.neighbors{currentTargetId} = neighborSet;
            end
        
        end
        
        
        
        
    end
end

