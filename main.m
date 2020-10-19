clear all
close all
clc


% create a random graph
graph = Graph(0); % create the graph object : index 0 for main graph 

%Case 1: 7, 7, 0.7, 2;   Case 2: 250, 7, 0.45, 2;   Case 3: 500, 10, 0.45, 4;    Case 4: 77, 10, 0.45, 4 
rand('seed', 77); %% 7, 250, 500
numOfTargets = 10; % 7, 10
dimentionOfSpace = 2;
sizeOfSpace = 1;  % 1x1 square
communicationRadius = 0.45; %0.7,  0.45

% load the agents
numOfAgents = 4; %2, 4

% Control Options: RHC-Re-Learn, Periodic, RHC-Periodic, RHC-Learn, RHC, BDC-Periodic, BDC, fixed, random
controlMethod = 'BDC-Periodic';
targetControllersEnabled = true;

% create the graph and draw
graph = graph.loadARandomGraph(numOfTargets, numOfAgents, dimentionOfSpace, sizeOfSpace, communicationRadius, targetControllersEnabled);
graph.drawGraph(1);

% if we try the periodic controller
if isequal(controlMethod,'RHC-Periodic') | isequal(controlMethod,'BDC-Periodic') | isequal(controlMethod,'Periodic')
    graph.assignCyclesToAgents(controlMethod);
end

% BDC learning Thresholds
if isequal(controlMethod,'BDC') | isequal(controlMethod,'BDC-Periodic')
    graph.loadBDCThresholds('results/Case3/ST/Case3_RHC_Periodic.mat');
end




timeResolution = 0.001;
periodT = 50; % 50, 200, 500, 
fixedHorizonValue = 2.5; %2.50727
cost = 0;
controlCost = 0;
textHandle1 = text();  textHandle2 = text(); textHandle3 = text(); 
plotMode = true; 
videoMode = true;
perturbMode = []; %[1 40; 0 2.1; 1 7; 0 7.1; inf 1]; % [1  2; 0  2.1; 1  7; 0  7.1]: activate at 2 and 7  for a duration of 0.1 s 
data = zeros(periodT/timeResolution+1, numOfTargets, 5);
dataAgent = zeros(periodT/timeResolution+1, numOfAgents, 2);
timeSteps = 1;
frameCount = 1;
totalDistance = 0;
maxCovarianceRecorded = zeros(2,numOfTargets); % maxOmega, time 
for t = 0:timeResolution:periodT
    
    % target update
    for i = 1:1:numOfTargets
        [costVal, dataVal, costVal2] = graph.targets(i).update(timeResolution,plotMode,t); 
        cost = cost + costVal;
        if targetControllersEnabled
            controlCost = controlCost + costVal2;
        end
        
        if plotMode
            data(timeSteps,i,:) = dataVal; %[phi ,phiHat ,Omega ,r ,eta]
            if maxCovarianceRecorded(1,i) < graph.targets(i).Omega
                maxCovarianceRecorded(1,i) = graph.targets(i).Omega; % re-assign
                maxCovarianceRecorded(2,i) = t; % current time
            end
        end
        
        
    end
    
    % agent update
    horizon = min(fixedHorizonValue,periodT-t);
    for i = 1:1:numOfAgents
        distanceTraveled = graph.agents(i).update(timeResolution,t,graph,horizon,controlMethod); 
        if plotMode
            dataAgent(timeSteps,i,:) = [distanceTraveled,0]; %[distance ,nullForNow]
            totalDistance = totalDistance + distanceTraveled;
        end
    end
    
    % display & video
    if rem(t,0.1)==0 
        time = t
        if plotMode 
            graph.updateGraph(1,t);
            if t~=0
                delete(textHandle1);
                textHandle1 = text(0.05,0.96,['Cost: ',num2str(cost/t,5)],'Color','k','FontSize',10);
            end
            delete(textHandle2);
            textHandle2 = text(0.05,0.92,['Time: ',num2str(t,5)],'Color','k','FontSize',10);
            if targetControllersEnabled
                delete(textHandle3);
                textHandle3 = text(0.05,0.88,['Ctrl. Cost: ',num2str(periodT*controlCost/t,5)],'Color','k','FontSize',10);
            end
        end
        
        if videoMode
            frameArray(frameCount) = getframe(gcf);
            frameCount = frameCount + 1;
        end
        
        if length(perturbMode)>0
            count = perturbMode(end,end);
            if perturbMode(count,2)==t
                action = perturbMode(count,1);
                graph.triggerExternalStatePerturbation(action);
                perturbMode(end,end) = perturbMode(end,end) + 1; % increasing the count
            end
        end
    end
    
    timeSteps = timeSteps + 1;
end

disp(['Total Distance = ',num2str(totalDistance)]);
[maxOmega, maxTarget] = max(maxCovarianceRecorded(1,:));
maxTime = maxCovarianceRecorded(2,maxTarget);
disp(['Max Omega: ',num2str(maxOmega),'; at Target: ',num2str(maxTarget),'; at Time: ',num2str(maxTime)]);

run('graphicsOfSimulation.m')

