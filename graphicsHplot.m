clear all
close all
clc

costArray = [];
frames = [];
for fixedHorizonValue = 10:5:25 
    fixedHorizonValue 
    
    % create a random graph
    graph = Graph(0); % create the graph object : index 0 for main graph 

    %Case 1: 7, 7, 0.7, 2;   Case 2: 250, 7, 0.45, 2;   Case 3: 500, 10, 0.45, 4
    rand('seed', 250); %% 7, 250, 500
    numOfTargets = 7; % 7, 10
    dimentionOfSpace = 2;
    sizeOfSpace = 1;  % 1x1 square
    communicationRadius = 0.45; %0.7,  0.45

    % load the agents
    numOfAgents = 2; %2, 4

    % create the graph and draw
    graph = graph.loadARandomGraph(numOfTargets, numOfAgents, dimentionOfSpace, sizeOfSpace, communicationRadius)
    graph.drawGraph(1);

    timeResolution = 0.001;
    periodT = 50; % 50

    cost = 0;  
    textHandle1 = text();  textHandle2 = text(); 
    plotMode = false; 
    videoMode = false;
    data = zeros(periodT/timeResolution+1, numOfTargets, 3);
    timeSteps = 1;
    frameCount = 1;
    for t = 0:timeResolution:periodT

        % target update
        for i = 1:1:numOfTargets
            [costVal, dataVal] = graph.targets(i).update(timeResolution,plotMode); 
            cost = cost + costVal;
            if plotMode
                data(timeSteps,i,:) = dataVal; %{phi,phiHat,Omega]
            end
        end
 
        % agent display
        horizon = min(fixedHorizonValue,periodT-t);
        for i = 1:1:numOfAgents
            graph.agents(i).update(timeResolution,t,graph,horizon);     
        end

        % display
        if rem(t,0.1)==0
            time = t
%             graph.updateGraph(1,t);
%             delete(textHandle1);
%             textHandle1 = text(0.05,0.96,['Cost: ',num2str(cost,5)],'Color','k','FontSize',10);
%             delete(textHandle2);
%             textHandle2 = text(0.05,0.92,['Time: ',num2str(t,5)],'Color','k','FontSize',10);
% 
%             frameArray(frameCount) = getframe(gcf);
%             frameCount = frameCount + 1;
        end

        timeSteps = timeSteps + 1;
    end
    
    costArray = [costArray; fixedHorizonValue, cost];
%     frames = [frames; frameArray];
end



% plotting
load('Hdata1.mat')
c1 = costArray;
load('Hplot2.mat')
c2 = costArray;
c3 = [2.5, 3577.3]
load('Hdata3.mat')
c4 = costArray;
load('Hplot4.mat')
c5 = costArray;
load('Hplot5.mat')
c6 = costArray;
load('Hplot6.mat')
c7 = costArray;
c = [c1;c2;c3;c4;c5;c6;c7]
figure
plot(c(2:end,1),c(2:end,2))
grid on
xlabel('Planning Horizon Limit: $H$','Interpreter','latex')
ylabel('Global Objective: $J_T$','Interpreter','latex')

    
 