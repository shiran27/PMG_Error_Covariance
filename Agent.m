classdef Agent < handle
    % Class for a generic agent
    
    properties
        index
        position
        orientation % for graphics
        speed
        residingTarget
        
        agentMode % dwelling or travelling
        nextEventTime % time to next mode change
        coverednessEventTriggered
        
        assignedCycle
        entryPathToCycle
        cycleState % [entryPathState,assignedCycleState]
        dwellTimes
        
        learnMode % used in RHC mode 
        pastData % used in RHC mode just tot store data
        
        slidingWindowSize
        pastInputOutputData % struct with fields: inputs, outputs, data count
        neuralNetworks
        executionTimes % [time, RHCPId, executionTime] array
        
        graphicHandles
        
    end
    
    methods
        
        function obj = Agent(index,position,speed,residingTarget,numOfTargets)
            %UNTITLED2 Construct an instance of this class
            obj.index = index;
            obj.position = position;
            obj.speed = speed;
            obj.residingTarget = residingTarget;
            obj.orientation = 0;
            obj.coverednessEventTriggered = false;
            
            obj.learnMode = false; % for RHC
             
            % for RHC-Learn
            obj.slidingWindowSize = 75; %25
            obj.pastInputOutputData.inputs = cell(numOfTargets,2);
            obj.pastInputOutputData.outputs = cell(numOfTargets,2);
            obj.pastInputOutputData.dataCounts = zeros(numOfTargets,2);

            hiddenLayerSize = [10];%[10]; % one layer with 10 neurons
            net = patternnet(hiddenLayerSize);
            nets = cell(numOfTargets,2);
            [nets{:}] = deal(net); % repmat version
            obj.neuralNetworks = nets;

            
        end
        
        
        function distanceTraveled = update(obj,deltaT,currentTime,graph,horizon,controlMethod)
            distanceTraveled = 0;
            if currentTime ~= 0 & ~obj.coverednessEventTriggered
                if currentTime >  obj.nextEventTime & obj.agentMode == 1 % travel time exceeded
                    % pick a dwell time and start dwelling 
                    disp(['Agent ',num2str(obj.index),' finished travelling at t = ',num2str(currentTime)])
                    obj.agentMode = 0;
                    obj.residingTarget = obj.residingTarget(2);
                    graph.targets(obj.residingTarget).residingAgents(1) = obj.index;
                    graph.targets(obj.residingTarget).timeSinceLastEvent = 0;
                    graph.targets(obj.residingTarget).OmegaAtLastEvent = graph.targets(obj.residingTarget).Omega;
                    
                    % % Important:
%                     tic;
                    t_0 = cputime;
                    dwellTimeAhead = obj.dwellTimeAhead(graph,horizon,controlMethod);
                    executionTime = cputime - t_0;
%                     executionTime = toc;
                    if dwellTimeAhead>0.001
                        obj.executionTimes(end,end) = obj.executionTimes(end,end) + 1;
                        count = obj.executionTimes(end,end);
                        obj.executionTimes(count,:) = [currentTime, 1, executionTime];
% %                         obj.executionTimes = [obj.executionTimes; currentTime, 1, executionTime];
                    end
                    % % End Important.
                    
                    
                    obj.nextEventTime = currentTime + dwellTimeAhead;
                    distanceTraveled = norm(obj.position - graph.targets(obj.residingTarget).position);
                    obj.position = graph.targets(obj.residingTarget).position;
                    
                elseif currentTime >  obj.nextEventTime & obj.agentMode == 0 % dwell time exceeded
                    % pick a neighbor to visit and start travelling
                    disp(['Agent ',num2str(obj.index),' finished dwelling at t = ',num2str(currentTime),' at taregt ',num2str(obj.residingTarget)])
                    
                    % % Important:
%                     tic;
                    t_0 = cputime;
                    [nextTargetId, travelTimeAhead] = obj.nextVisitTarget(graph,horizon,controlMethod); 
                    executionTime = cputime - t_0;
%                     executionTime = toc;
                    
                    
                    if nextTargetId ~= obj.residingTarget
                        if isempty(obj.executionTimes)
                            obj.executionTimes = zeros(3000,3);
                        end
                        obj.executionTimes(end,end) = obj.executionTimes(end,end) + 1;
                        count = obj.executionTimes(end,end);
                        obj.executionTimes(count,:) = [currentTime, 2, executionTime];
% %                         obj.executionTimes = [obj.executionTimes; currentTime, 2, executionTime];
                        % % End Important
                        
                        obj.agentMode = 1;
                        graph.targets(obj.residingTarget).residingAgents = [];
                        graph.targets(obj.residingTarget).timeSinceLastEvent = 0;
                        graph.targets(obj.residingTarget).OmegaAtLastEvent = graph.targets(obj.residingTarget).Omega;
                        
                        obj.residingTarget = [obj.residingTarget, nextTargetId];
                        obj.nextEventTime = currentTime + travelTimeAhead;
                    
                        theta = graph.angleMatrix(obj.residingTarget(1),obj.residingTarget(2));
                        obj.position = obj.position + obj.speed*deltaT*[cos(theta), sin(theta)];
                        obj.orientation = theta;
                        
                        % trigger an event at neighbor targets of i telling i is now uncovered!
                        % also trigger an event telling at neighbors of j telling it is now covered!
                        graph.triggerCoverednessEvents(obj.residingTarget(1),obj.residingTarget(2));
                    else
                        disp('Error: No neighbor!')
                    end
                elseif obj.agentMode == 1 % travelling
                    % keep travelling
                    theta = obj.orientation; %graph.angleMatrix(obj.residingTarget(1),obj.residingTarget(2));
                    distanceTraveled = obj.speed*deltaT;
                    obj.position = obj.position + distanceTraveled*[cos(theta), sin(theta)];
                else
                    % keep dwelling
                end
            else % only at begining: t = 0 or when coveredness related event is triggered
                obj.agentMode = 0;  % in dwell mode 
                
                % % Important:
%                 tic;
                t_0 = cputime;
                dwellTimeAhead = obj.dwellTimeAhead(graph,horizon,controlMethod);
                executionTime = cputime - t_0;
%                 executionTime = toc;
                if dwellTimeAhead > 0.001
                    if isempty(obj.executionTimes)
                        obj.executionTimes = zeros(3000,3);
                        obj.executionTimes(end,end) = obj.executionTimes(end,end) + 1;
                        count = obj.executionTimes(end,end);
                        obj.executionTimes(count,:) = [currentTime, 1, executionTime];
                    else
                        obj.executionTimes(end,end) = obj.executionTimes(end,end) + 1;
                        count = obj.executionTimes(end,end);
                        obj.executionTimes(count,:) = [currentTime, 1, executionTime];
% %                         obj.executionTimes = [obj.executionTimes; currentTime, 1, executionTime];
                    end
                end
                % % End Important
                
                obj.nextEventTime = currentTime + dwellTimeAhead;
                if obj.coverednessEventTriggered
                    disp(['Agent ',num2str(obj.index),' recomputed controls!']);
                    obj.coverednessEventTriggered = false;
                end
                
            end
        end
        
        
        
        function dwellTime = dwellTimeAhead(obj, graph, horizon, controlMethod)
            if isequal(controlMethod,'random')
                dwellTime = 5*rand(1,1);
                
                
            elseif isequal(controlMethod,'fixed')
                dwellTime = 2.5;
                
                
            elseif isequal(controlMethod,'BDC') | isequal(controlMethod,'BDC-Periodic')
                currentTargetId = obj.residingTarget;
                epsilon = graph.targets(currentTargetId).BDCThreshold; % if epsilon = 0.1, get Omega/Omega_ss to 1.1 (or 0.9 if Omgea_0 is less than Omega_ss)
                dwellTime = graph.targets(currentTargetId).dwellTimeToReduceCovarianceUptoFraction(epsilon);
                
                
            elseif isequal(controlMethod,'RHC')
                currentTargetId = obj.residingTarget;
                % finding feasible neighbors
                feasibleTargets = [];
                for targetId = 1:1:length(graph.targets)
                    if graph.distanceMatrix(currentTargetId,targetId) < 1.5 & isempty(graph.targets(targetId).residingAgents) & currentTargetId ~= targetId 
                        covered = false;
                        for agentId = 1:1:length(graph.agents)
                            if length(graph.agents(agentId).residingTarget)>1 & obj.index ~= agentId
                                if graph.agents(agentId).residingTarget(2) == targetId
                                    covered = true;
                                    break
                                end
                            end
                        end
                        if ~covered
                            feasibleTargets = [feasibleTargets, targetId];
                        end
                    end
                end
                
                nextTargetId = currentTargetId;
                
                if obj.learnMode 
                    dataFile = zeros(length(feasibleTargets)+1,4);
                    dataFile(1,:) = [currentTargetId, graph.targets(currentTargetId).Omega, 0, 0]; %targetId_i, Omega_i, ~, j^*
                    dataCount = 2;
                end
                
                dwellTime = 0.001;
                minCost = inf;
                for targetId = feasibleTargets
                    % solve for the optimum dwell time to be spent at target i=currentTargetId  and target j=targetId upun visiting.
                    % associated cost is used to find the best dwell time at target i  
                    [cost, u_i, u_j] = obj.solveRHCP1(targetId,feasibleTargets,graph,horizon); 
                    
                    if obj.learnMode
                        dataFile(dataCount,:) = [targetId, graph.targets(targetId).Omega, u_i, u_j]; % j, Omega_j, u_i
                        dataCount = dataCount + 1;
                    end
                    
                    if cost < minCost 
                        minCost = cost;
                        nextTargetId = targetId;
                        dwellTime = u_i;
                    end
                end
                
                if obj.learnMode & currentTargetId ~= nextTargetId
                    dataFile(1,4) = nextTargetId;
                    obj.pastData{length(obj.pastData)+1} = dataFile;
                end
                
            elseif isequal(controlMethod,'RHC-Learn') | isequal(controlMethod,'RHC-Re-Learn')  
                
                currentTargetId = obj.residingTarget;
                % finding feasible neighbors
                feasibleTargets = [];
                for targetId = 1:1:length(graph.targets)
                    if graph.distanceMatrix(currentTargetId,targetId) < 1.5 & isempty(graph.targets(targetId).residingAgents) & currentTargetId ~= targetId 
                        covered = false;
                        for agentId = 1:1:length(graph.agents)
                            if length(graph.agents(agentId).residingTarget)>1 & obj.index ~= agentId
                                if graph.agents(agentId).residingTarget(2) == targetId
                                    covered = true;
                                    break
                                end
                            end
                        end
                        if ~covered
                            feasibleTargets = [feasibleTargets, targetId];
                        end
                    end
                end
                
                
                stillLearning = obj.pastInputOutputData.dataCounts(currentTargetId, 1) <= obj.slidingWindowSize;
                needsReLearning = false;
                
                if ~stillLearning % learning has finished !
                    
                    % collect input data
                    input = zeros(length(graph.neighbors{currentTargetId}),1);
                    input(1) =  graph.targets(currentTargetId).Omega; % current target covariance
                    for j = feasibleTargets
                        Omega_j = graph.targets(j).Omega;
                        input(graph.neighbors{currentTargetId}==j) = Omega_j;
                    end
                    
                    % apply that to the nueral network
                    net = obj.neuralNetworks{currentTargetId, 1};
                    result = net(input);
                    jStar = graph.neighbors{currentTargetId}(vec2ind(result)+1);
                    
                    % if some ambiguity left: evaluate RHCP1 for all the ones and Re-learn the neural net
                    if (1-max(result))>0.2 & isequal(controlMethod,'RHC-Re-Learn') % 0.1
                        needsReLearning = true;
                    else
                        % if suceess e > 0.1 evaluate only the RHCP1 for the j* found
                        [~, dwellTime, ~] = obj.solveRHCP1(jStar,feasibleTargets,graph,horizon);     
                    end
                    
                end
                
                if stillLearning | needsReLearning
                    
                    % evaluate RHCP1 for all the neighbors
                    dataFile = zeros(length(feasibleTargets)+1,4);
                    dataFile(1,:) = [currentTargetId, graph.targets(currentTargetId).Omega, 0, 0]; %targetId_i, Omega_i, ~, j^*
                    dataCount = 2;

                    nextTargetId = currentTargetId;
                    dwellTime = 0.001;
                    minCost = inf;
                    for targetId = feasibleTargets
                        [cost, u_i, u_j] = obj.solveRHCP1(targetId,feasibleTargets,graph,horizon); 
                        
                        dataFile(dataCount,:) = [targetId, graph.targets(targetId).Omega, u_i, u_j]; % j, Omega_j, u_i
                        dataCount = dataCount + 1;
                    
                        if cost < minCost 
                            minCost = cost;
                            nextTargetId = targetId;
                            dwellTime = u_i;
                        end
                    end
                
                    if currentTargetId ~= nextTargetId
                        dataFile(1,4) = nextTargetId;
                        obj.updatePastInputOutputData(graph, dataFile); % learning happaens inside this if data filled 
                    end

                end
          
                
            
            elseif isequal(controlMethod,'RHC-Periodic')
                currentTargetId = obj.residingTarget;
                % finding the next target to visit
                if obj.cycleState(1) < 0 % have to enter the cycle yet
                    nextTarget = obj.entryPathToCycle(end+obj.cycleState(1)+1);
                    feasibleTargets = [nextTarget];
                else % already in the cycle
                    if obj.cycleState(2)==length(obj.assignedCycle.targetList)
                        nextTragetIndex = 1;
                    else
                        nextTragetIndex = obj.cycleState(2) + 1;
                    end
                    nextTarget = obj.assignedCycle.targetList(nextTragetIndex);
                    
                    feasibleTargets = [];
                    for targetId = 1:1:length(graph.targets)
                        if graph.distanceMatrix(currentTargetId,targetId) < 1.5 & isempty(graph.targets(targetId).residingAgents) & currentTargetId ~= targetId & sum(obj.assignedCycle.targetList==targetId)>0
                            covered = false;
                            for agentId = 1:1:length(graph.agents)
                                if length(graph.agents(agentId).residingTarget)>1 & obj.index ~= agentId
                                    if graph.agents(agentId).residingTarget(2) == targetId
                                        covered = true;
                                        break
                                    end
                                end
                            end
                            if ~covered
                                feasibleTargets = [feasibleTargets, targetId];
                            end
                        end
                    end
                end
                
                if isempty(feasibleTargets)
                    dwellTime = 1;
                else
                    [~, dwellTime, ~] = obj.solveRHCP1(nextTarget,feasibleTargets,graph,horizon); 
                end
                
                
            elseif isequal(controlMethod,'Periodic')
                if obj.cycleState(1) < 0 % have to enter the cycle yet
                    dwellTime = 0;
                else % already in the cycle
                    currentTargetIndex = obj.cycleState(2);
                    dwellTime = obj.dwellTimes(currentTargetIndex);
                end
                
                
            end
        end
        
        
        
        function [nextTargetId, travelTimeAhead] = nextVisitTarget(obj, graph, horizon, controlMethod)
            currentTargetId = obj.residingTarget;
            
            if isequal(controlMethod,'RHC-Periodic') | isequal(controlMethod,'BDC-Periodic') | isequal(controlMethod,'Periodic')
                % finding the next target to visit
                if obj.cycleState(1) < 0 % have to enter the cycle
                    nextTargetId = obj.entryPathToCycle(end+obj.cycleState(1)+1);
                    obj.cycleState(1) = obj.cycleState(1) + 1; 
                else % already in the cycle
                    if obj.cycleState(2)==length(obj.assignedCycle.targetList)
                        obj.cycleState(2) = 1;
                    else
                        obj.cycleState(2) = obj.cycleState(2) + 1;
                    end
                    nextTargetId = obj.assignedCycle.targetList(obj.cycleState(2));
                end
                travelTimeAhead = graph.distanceMatrix(currentTargetId,nextTargetId)/obj.speed;
                return;
                
                
            end
            % otherwise:
            
            
            % feasible neighbors
            feasibleTargets = [];
            for targetId = 1:1:length(graph.targets)
                if graph.distanceMatrix(currentTargetId,targetId) < 1.5 & isempty(graph.targets(targetId).residingAgents) & currentTargetId ~= targetId 
                    covered = false;
                    for agentId = 1:1:length(graph.agents)
                        if length(graph.agents(agentId).residingTarget)>1 & obj.index ~= agentId
                            if graph.agents(agentId).residingTarget(2) == targetId
                                covered = true;
                                break
                            end
                        end
                    end
                    
                    if ~covered
                        feasibleTargets = [feasibleTargets, targetId];
                    end
                end
            end
%             feasibleTargets
            

            if isequal(controlMethod,'RHC-Learn') | isequal(controlMethod,'RHC-Re-Learn')
                
                stillLearning = obj.pastInputOutputData.dataCounts(currentTargetId, 2) <= obj.slidingWindowSize;
                needsReLearning = false;
                
                if ~stillLearning % learning has finished !
                    
                    % collect input data
                    input = zeros(length(graph.neighbors{currentTargetId}),1);
                    input(1) =  graph.targets(currentTargetId).Omega; % current target covariance
                    for j = feasibleTargets
                        Omega_j = graph.targets(j).Omega;
                        input(graph.neighbors{currentTargetId}==j) = Omega_j;
                    end
                    
                    % apply that to the nueral network
                    net = obj.neuralNetworks{currentTargetId, 2};
                    result = net(input);
                    
                    % if some ambiguity left: evaluate RHCP1 for all the ones and Re-learn the neural net
                    if (1-max(result))>0.2 & isequal(controlMethod,'RHC-Re-Learn') % 0.1, 0.2
                        needsReLearning = true;
                    else
                        % if suceess e > 0.1 we know the RHCP2 solution is directly j* found
                        nextTargetId = graph.neighbors{currentTargetId}(vec2ind(result)+1);
                    end
                    
                end
                
                if stillLearning | needsReLearning
                    
                    % evaluate RHCP2 for all the neighbors
                    dataFile = zeros(length(feasibleTargets)+1,3);
                    dataFile(1,:) = [currentTargetId, graph.targets(currentTargetId).Omega, 0]; %targetId_i, Omega_i, j^*
                    dataCount = 2;

                    nextTargetId = currentTargetId;
                    minCost = inf;
                    for targetId = feasibleTargets
                        [cost, u_j] = obj.solveRHCP2(targetId,feasibleTargets,graph,horizon); 
                    
                        dataFile(dataCount,:) = [targetId, graph.targets(targetId).Omega, u_j]; % j, Omega_j, u_j
                        dataCount = dataCount + 1;
                                        
                        if cost < minCost 
                            minCost = cost;
                            nextTargetId = targetId;
                        end
                    end
                
                    if nextTargetId ~= currentTargetId 
                        dataFile(1,3) = nextTargetId;
                        obj.updatePastInputOutputData(graph,dataFile);
                    end

                end
                
                travelTimeAhead = graph.distanceMatrix(currentTargetId,nextTargetId)/obj.speed;
                return
                
            end
            % otherwise



            if obj.learnMode 
                dataFile = zeros(length(feasibleTargets)+1,3);
                dataFile(1,:) = [currentTargetId, graph.targets(currentTargetId).Omega, 0]; %targetId_i, Omega_i, j^*
                dataCount = 2;
            end
            
            nextTargetId = currentTargetId;
            minCost = inf;
            for targetId = feasibleTargets
                if  isequal(controlMethod,'BDC') | isequal(controlMethod,'fixed') | isequal(controlMethod,'random') %isequal(controlMethod,'highestCovariance')
                    cost = -1*graph.targets(targetId).Omega; % select the target with the highest covariance out of feasible targets
                
                elseif isequal(controlMethod,'RHC')
                    % solve for the optimum dwell time to be spent at target j=targetId upun visiting.
                    % associated cost is used to fined th ebest taget to visit next 
                    [cost, u_j] = obj.solveRHCP2(targetId,feasibleTargets,graph,horizon); 
                    
                    if obj.learnMode
                        dataFile(dataCount,:) = [targetId, graph.targets(targetId).Omega, u_j]; % j, Omega_j, u_j
                        dataCount = dataCount + 1;
                    end    
                end
                
                if cost < minCost 
                    minCost = cost;
                    nextTargetId = targetId;
                end
            end
            
            if obj.learnMode & nextTargetId ~= currentTargetId 
                dataFile(1,3) = nextTargetId;
                obj.pastData{length(obj.pastData)+1} = dataFile;
            end
            
            travelTimeAhead = graph.distanceMatrix(currentTargetId,nextTargetId)/obj.speed;  
        end
        
        
        
        function outputArg = updatePastInputOutputData(obj, graph, dataBlock)
            i = dataBlock(1,1); % current target
            jStar = dataBlock(1,size(dataBlock,2)); % chosen next target
            RHCPId = 1*(size(dataBlock,2)==4) + 2*(size(dataBlock,2)==3);%1 or 2
            
            obj.pastInputOutputData.dataCounts(i,RHCPId) = obj.pastInputOutputData.dataCounts(i,RHCPId) + 1; 

            % output
            neighborIndex = find(graph.neighbors{i}==jStar)-1;
            output = zeros(length(graph.neighbors{i})-1,1);
            output(neighborIndex) = 1;
            
            
            % input
            input = zeros(length(graph.neighbors{i}),1);
            for jId = 1:1:size(dataBlock,1)
                j = dataBlock(jId,1);
                Omega_j = dataBlock(jId,2);
                input(graph.neighbors{i}==j) = Omega_j;
            end
            
            
            % loading input and output data
            inputsMatrix = obj.pastInputOutputData.inputs{i,RHCPId};
            outputsMatrix = obj.pastInputOutputData.outputs{i,RHCPId};
            if size(inputsMatrix,2) == obj.slidingWindowSize
                inputsMatrix = [inputsMatrix(:,2:end), input];
                outputsMatrix = [outputsMatrix(:,2:end), output];
            else
                inputsMatrix = [inputsMatrix, input];
                outputsMatrix = [outputsMatrix, output];
            end
            obj.pastInputOutputData.inputs{i,RHCPId} = inputsMatrix;
            obj.pastInputOutputData.outputs{i,RHCPId} = outputsMatrix;

            
            % learning the net if data is filled
            if obj.pastInputOutputData.dataCounts(i,RHCPId) > obj.slidingWindowSize
                net = obj.neuralNetworks{i,RHCPId};
                [flag,~,~] = isconfigured(net);
                if ~flag
                    disp(['Target ',num2str(i),' RHCP ',num2str(RHCPId),' learning began']);
                    net = configure(net,inputsMatrix,outputsMatrix);
                    net.trainParam.showWindow = 0;
                    net = train(net,inputsMatrix,outputsMatrix);
                    obj.neuralNetworks{i,RHCPId} = net;
                end
            end
            
            
        end
        
        
        function [minCost, uSol] = solveRHCP2(obj, nextTarget, neighbors, graph, horizon)
            currentTarget = obj.residingTarget;
            rho = graph.distanceMatrix(currentTarget,nextTarget);
            H = horizon; % horizon length
            if H<rho % do not waste time!
                minCost = inf;
                uSol = 0;
                return;
            end
                
            
            syms u; % Dwell time to be spent at nextTarget 0\leq u \leq H-rho 
            
            %  This should include the cost of currenttarget
%             H_tilde = 0; % examle: exp(3*u) + exp(2*u);
            active = 0;
            inactive = 0;
            for l = [currentTarget, neighbors]
                if l  ~= nextTarget % currentTarget %
                    H_l = graph.targets(l).localObjectiveIdlePerod(0,rho+u); % startGap, duration 
%                     H_tilde = H_tilde + H_l;
                    inactive = inactive + H_l;
                elseif l == nextTarget % l = j case
                    H_j1 = graph.targets(l).localObjectiveIdlePerod(0,rho); % startGap, duration;
                    H_j2 = graph.targets(l).localObjectiveDwellPerod(rho,u); % startGap, duration;
%                     H_tilde = H_tilde + (H_j1 + H_j2);
                    active = active + H_j2;
                    inactive = inactive + H_j1;
                end
            end

            A = vpa(simplify(active),6); 
            B = vpa(simplify(inactive),6);
            J = vpa(-A/(A+B),6);
            
% %             Gradient Descent:
% %             tic
            uSol = obj.solveRHCP2UsingGD(J,u, H-rho);
% %             toc
% %             End gradient descent

            
% %             Alternative method:
%             tic
% % % % %             Adot = diff(A,u); Bdot = diff(B,u);
% % % % %             uSol = vpasolve(A*Bdot==B*Adot,u,[0,H-rho],'Random',true);
% % % % %             if isempty(uSol)
% % % % %                 syms lambda_1 lambda_2;
% % % % %                 Hamil = J + lambda_1*(-u) + lambda_2*(u-(H-rho));
% % % % %                 eqn1 = diff(Hamil,u)==0;
% % % % %                 eqn2 = lambda_1*(-u)==0;
% % % % %                 eqn3 = lambda_2*(u-(H-rho))==0;
% % % % %                 [uSol,~,~] = vpasolve([eqn1,eqn2,eqn3],[u,lambda_1,lambda_2],[0 H-rho;0 inf;0 inf],'Random',true);
% % % % %             end
%             toc
% %             End Alternative    

            if ~isempty(uSol)
                uSol = double(uSol);
                minCost = double(subs(J,u,uSol));
                %disp(['Solution Found: i=',num2str(currentTarget),', j=',num2str(nextTarget),', cost=',num2str(minCost),', u=',num2str(uSol)])
%                 active;
%                 inactive;
            else
                uSol = 0;
                minCost = double(subs(J,u,0));
                %disp(['No Solution! i=',num2str(currentTarget),', j=',num2str(nextTarget),', cost=',num2str(minCost)])
                
%                 vpasolve(A*Bdot==B*Adot,u,[0,10])
%                 active
%                 inactive              
            end
        end
        
        
        
        function [minCost, u_iSol, u_jSol] = solveRHCP1(obj, nextTarget, neighbors, graph, horizon)
            currentTarget = obj.residingTarget;
            rho = graph.distanceMatrix(currentTarget,nextTarget);
            H = horizon; % horizon length
            if H<rho % do not waste time!
                u_iSol = 1;
                u_jSol = 0;
                minCost = inf;
                return;
            end
            
            syms u_i u_j; % Dwell time to be spent at currentTarget
            
            %  This should include the cost of currenttarget
%             H_tilde = 0; % examle: exp(3*u) + exp(2*u);
            active = 0;
            inactive = 0;
            for l = [currentTarget, neighbors]
                if l == currentTarget
                    H_i1 = graph.targets(l).localObjectiveDwellPerod(0,u_i); % startGap, duration;
                    H_i2 = graph.targets(l).localObjectiveIdlePerod(u_i,u_j+rho); % startGap, duration;
%                     H_tilde = H_tilde + H_i1 + H_i2;
                    active = active + H_i1;
                    inactive = inactive + H_i2;
                elseif l ~= nextTarget
                    H_l = graph.targets(l).localObjectiveIdlePerod(0,u_i+rho+u_j); % startGap, duration 
%                     H_tilde = H_tilde + H_l;
                    inactive = inactive + H_l;
                else % l = j case
                    H_j1 = graph.targets(l).localObjectiveIdlePerod(0,u_i+rho); % startGap, duration;
                    H_j2 = graph.targets(l).localObjectiveDwellPerod(u_i+rho,u_j); % startGap, duration;
%                     H_tilde = H_tilde + (H_j1 + H_j2);
                    active = active + H_j2;
                    inactive = inactive + H_j1;
                end
            end
            
            A = vpa(simplify(active),6);
            B = vpa(simplify(inactive),6);
            J = vpa(simplify(-A/(B+A)),6);
            
% %             Gradient Descent:
%             tic
            [u_iSol, u_jSol] = obj.solveRHCP1UsingGD(J, u_i, u_j, H-rho);
%             toc
% %             End gradient descent
            
            
% %             % Alternative method
% % %             syms lambda_1 lambda_2 lambda_3; % Dwell time to be spent at currentTarget
% % %             Hamil = J + lambda_1*(-u_i) + lambda_2*(-u_j) + lambda_3*(u_i+u_j-(H-rho));
% % %             eqn1 = diff(Hamil,u_i)==0;
% % %             eqn2 = diff(Hamil,u_j)==0;
% % %             eqn3 = lambda_1*(-u_i)==0;
% % %             eqn4 = lambda_2*(-u_j)==0;
% % %             eqn5 = lambda_3*(u_i+u_j-(H-rho))==0;
% % %             [u_iSol,u_jSol,~,~,~] = vpasolve([eqn1,eqn2,eqn3,eqn4,eqn5],[u_i,u_j,lambda_1,lambda_2,lambda_3],[0 H-rho;0 H-rho;0 inf;0 inf;0 inf],'Random',true);
% % %             
% % %             if isempty(u_iSol) | isempty(u_jSol) 
% % %                 syms dummy;
% % %                 Adot_i = diff(A,u_i); Bdot_i = diff(B,u_i);
% % %                 Adot_j = diff(A,u_j); Bdot_j = diff(B,u_j);
% % %                 eqn1 = Adot_i*B == Bdot_i*A;
% % %                 eqn2 = Adot_j*B == Bdot_j*A;
% % %                 eqn3 = u_i+u_j+dummy == H-rho;
% % %                 [u_iSol,u_jSol,~] = vpasolve([eqn1,eqn2,eqn3],[u_i,u_j,dummy],[0,H-rho;0,H-rho;0,H-rho],'Random',true);
% % %             end
% %             % end alternate method
            
            
            if ~isempty(u_iSol) & ~isempty(u_jSol) 
                u_iSol = double(u_iSol);
                u_jSol = double(u_jSol);
                minCost = double(subs(J,[u_i,u_j],[u_iSol,u_jSol]));
%                 disp(['Solution Found: i=',num2str(currentTarget),', j=',num2str(nextTarget),', cost=',num2str(minCost),', u_i=',num2str(u_iSol),', u_j=',num2str(u_jSol)])
%                 A
%                 B
%                 H
%                 J
            else
                u_iSol = 1;
                u_jSol = 0;
                minCost = double(subs(J,[u_i,u_j],[u_iSol,u_jSol]));
                %disp(['No Solution! i=',num2str(currentTarget),', j=',num2str(nextTarget),', cost=',num2str(minCost)])
%                 H-rho<0
%                 A
%                 B
%                 H
%                 J
                
            end
        end
        
        
        function uSol = solveRHCP2UsingGD(~, J, u, uLim)
%             JFun = matlabFunction(J);
            DJFun = matlabFunction(diff(J,u));
            
            u_k = 0;
            beta = 0.1; % step size
            for k = 1:1:100000
                
                grad = DJFun(u_k);
                u_kNew = u_k - beta*grad;
                if u_kNew < 0
                    u_kNew = 0;
                elseif u_kNew > uLim
                    u_kNew = uLim;
                end

                if abs(u_k-u_kNew)<0.0001
                    uSol = u_kNew;
                    break
                end
                
                u_k = u_kNew;
            end
            
            uSol = u_kNew;
            
        end
        
        
        
        function [u_iSol, u_jSol] = solveRHCP1UsingGD(~, J, u_i, u_j, uLim)
%             JFun = matlabFunction(J);
            DJ_iFun = matlabFunction(diff(J,u_i));
            DJ_jFun = matlabFunction(diff(J,u_j));
            
            u_ik = 0; 
            u_jk = 0;
            beta = 0.1; % step size
            for k = 1:1:100000
                
                grad_i = DJ_iFun(u_ik, u_jk);
                grad_j = DJ_jFun(u_ik, u_jk);
                u_ikNew = u_ik - beta*grad_i;
                u_jkNew = u_jk - beta*grad_j;
                
                % projections
                if u_ikNew+u_jkNew > uLim
                    u_ikNew = 0.5*(u_ikNew - u_jkNew + uLim);
                    u_jkNew = 0.5*(u_jkNew - u_ikNew + uLim);
                end
                
                if u_ikNew < 0
                    u_ikNew = 0;
                end
                
                if u_jkNew < 0
                    u_jkNew = 0;
                end
                
                
                
                % checking convergence:
                if abs(u_ikNew-u_ik)<0.0001 & abs(u_jkNew-u_jk)<0.0001
                    u_iSol = u_ikNew;
                    u_jSol = u_jkNew;
                    break
                end
                
                u_ik = u_ikNew;
                u_jk = u_jkNew;
            end
            
            % did not converge:
            u_iSol = u_ikNew;
            u_jSol = u_jkNew;
                    
        end
        
    end
    
end