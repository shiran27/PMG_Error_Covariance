


%% Case 1: OPT1: CMO Method I
% % clear all
% % close all
% % clc
% % graph = Graph(0); % create the graph object : index 0 for main graph 
% % rand('seed', 4); 
% % numOfTargets = 5;
% % dimentionOfSpace = 2;
% % sizeOfSpace = 1;  
% % communicationRadius = 0.6; 
% % 
% % graph = graph.loadARandomGraph(numOfTargets, dimentionOfSpace, sizeOfSpace, communicationRadius)
% % graph.drawGraph(1);
% % 
% % cycleTemp = Cycle([5 2 4 1 3]) %[5 2 4 1 3 2 4 3]
% % cycleTemp = cycleTemp.updateEdgeList(graph);
% % cycleTemp.drawFullCycle(graph);
% % 
% % [gainValue, expandedCycle] = cycleTemp.greedy1OPT(graph);
% % expandedCycle.drawFullCycle(graph);
% End Case 1


%% Case 2: OPT1: CMO Method II
% % clear all
% % close all
% % clc
% % graph = Graph(0); % create the graph object : index 0 for main graph 
% % rand('seed', 4); 
% % numOfTargets = 5; 
% % dimentionOfSpace = 2;
% % sizeOfSpace = 1;  
% % communicationRadius = 0.7; 
% % 
% % graph = graph.loadARandomGraph(numOfTargets, dimentionOfSpace, sizeOfSpace, communicationRadius)
% % graph.drawGraph(1);
% % 
% % cycleTemp = Cycle([5 2 4 1 3 2 4 3])
% % cycleTemp = cycleTemp.updateEdgeList(graph);
% % cycleTemp.drawFullCycle(graph);
% % 
% % [gainValue, expandedCycle_11] = cycleTemp.greedyExpandCycle1(graph)
% % expandedCycle_11.drawFullCycle(graph)
% End Case 2


%% Case 3: OPT1: CMO Method III
% % clear all
% % close all
% % clc
% % graph = Graph(0); % create the graph object : index 0 for main graph 
% % rand('seed', 4); 
% % numOfTargets = 5; 
% % dimentionOfSpace = 2;
% % sizeOfSpace = 1;  
% % communicationRadius = 0.7; 
% % 
% % graph = graph.loadARandomGraph(numOfTargets, dimentionOfSpace, sizeOfSpace, communicationRadius)
% % graph.drawGraph(1);
% % 
% % cycleTemp = Cycle([5 2 4 1 3 2 4 3])
% % cycleTemp = cycleTemp.updateEdgeList(graph);
% % cycleTemp.drawFullCycle(graph);
% % 
% % [gainValue, expandedCycle_21] = cycleTemp.greedyExpandCycle2(graph)
% % expandedCycle_21.drawFullCycle(graph)
% End Case 3


%% Case 4: Greedy Alg 1
% % clear all
% % close all
% % clc
% % graph = Graph(0); % create the graph object : index 0 for main graph 
% % rand('seed', 270); %% 4,5,7 seed
% % numOfTargets = 7; % 15 % 12
% % dimentionOfSpace = 2;
% % sizeOfSpace = 1;  % 1x1 square
% % communicationRadius = 0.6; %0.7 %0.5 % 0.7
% % 
% % % create the graph and draw
% % graph = graph.loadARandomGraph(numOfTargets, dimentionOfSpace, sizeOfSpace, communicationRadius)
% % graph.drawGraph(1);
% % 
% % cycleTSP = Cycle(graph.findTSPCycle());
% % cycleTSP = cycleTSP.updateEdgeList(graph);
% % % cycleTSP.drawCycle(graph,'y');
% % cycleTSP.drawFullCycle(graph);
% % cycleTemp = cycleTSP;
% % 
% % currentTargetCycle = cycleTemp.executeGreedyExpansions(graph)
% End Case 4


%% Case 5: Greedy Alg 1
% % clear all
% % close all
% % clc
% % graph = Graph(0); % create the graph object : index 0 for main graph 
% % rand('seed', 100); %% 4,5,7 seed
% % numOfTargets = 7; % 15 % 12
% % dimentionOfSpace = 2;
% % sizeOfSpace = 1;  % 1x1 square
% % communicationRadius = 0.7; %0.7 %0.5 % 0.7% create the graph and draw
% % 
% % % create the graph and draw
% % graph = graph.loadARandomGraph(numOfTargets, dimentionOfSpace, sizeOfSpace, communicationRadius)
% % graph.drawGraph(1);
% % 
% % cycleTSP = Cycle(graph.findTSPCycle());
% % cycleTSP = cycleTSP.updateEdgeList(graph);
% % % cycleTSP.drawCycle(graph,'y');
% % cycleTSP.drawFullCycle(graph);
% % cycleTemp = cycleTSP;
% % 
% % currentTargetCycle = cycleTemp.executeGreedyExpansions(graph)
% End Case 5


%% Case 6: Greedy Alg 1
% % clear all
% % close all
% % clc
% % graph = Graph(0); % create the graph object : index 0 for main graph 
% % rand('seed', 50); %% 4,5,7 seed
% % numOfTargets = 7; % 15 % 12
% % dimentionOfSpace = 2;
% % sizeOfSpace = 1;  % 1x1 square
% % communicationRadius = 0.5; %0.7 %0.5 % 0.7
% % 
% % % create the graph and draw
% % graph = graph.loadARandomGraph(numOfTargets, dimentionOfSpace, sizeOfSpace, communicationRadius)
% % graph.drawGraph(1);
% % 
% % cycleTSP = Cycle(graph.findTSPCycle());
% % cycleTSP = cycleTSP.updateEdgeList(graph);
% % % cycleTSP.drawCycle(graph,'y');
% % cycleTSP.drawFullCycle(graph);
% % cycleTemp = cycleTSP;
% % 
% % currentTargetCycle = cycleTemp.executeGreedyExpansions(graph)
% End Case 6




%% Case 7
% % clear all
% % close all
% % clc
% % graph = Graph(0); % create the graph object : index 0 for main graph 
% % rand('seed', 27); %% 4,5,7 seed
% % numOfTargets = 7; % 15 % 12
% % dimentionOfSpace = 2;
% % sizeOfSpace = 1;  % 1x1 square
% % communicationRadius = 0.6; %0.7 %0.5 % 0.7
% % 
% % graph = graph.loadARandomGraph(numOfTargets, dimentionOfSpace, sizeOfSpace, communicationRadius)
% % graph.drawGraph(1);
% % 
% % [dis, cyc, endCyc] = graph.constructACycleStartingFromANode(1,'fast');
% % % minDistances = [Inf   49.1277   47.5390   11.2918   46.5593   42.1195   19.5325];
% End Case 7




%% Case 8
% % clear all
% % close all
% % clc
% % graph = Graph(0); % create the graph object : index 0 for main graph 
% % rand('seed', 4); %% 4,5,7 seed
% % numOfTargets = 15; % 15 % 12
% % dimentionOfSpace = 2;
% % sizeOfSpace = 1;  % 1x1 square
% % communicationRadius = 0.7; %0.7 %0.5 % 0.7
% % 
% % graph = graph.loadARandomGraph(numOfTargets, dimentionOfSpace, sizeOfSpace, communicationRadius)
% % graph.drawGraph(1);
% % numOfAgents = 3;
% % 
% % % Option 1: Clustering based on 'smart' similarity values
% % graph = graph.loadShortestPaths2('fast'); % try 'slow' if have patience 
% % partitions = graph.partitionTheGraph(numOfAgents,'smallest lower bound'); % spectral clustering approach 2
% % 
% % % Option 2: Clustering based on shortest distance based similarity values
% % % partitions = graph.partitionTheGraph(numOfAgents,'shortest path'); % spectral clustering approach 1
% % 
% % subGraphs = graph.generateSubGraphs(partitions);
% % cycles = [];
% % for subGraph = subGraphs
% %     
% %     cycleTSP = Cycle(subGraph.findTSPCycle()); % Cycle([8,1,10,8,1,10,8,15,13]);
% %     cycleTSP = cycleTSP.updateEdgeList(subGraph);
% %     
% %     greedyExpandedCycle = cycleTSP.executeGreedyExpansions(subGraph);
% %     cycles = [cycles, greedyExpandedCycle];  
% % 
% % end
% % %drawing the initial subgraphs and cycles
% % graph.drawAll(subGraphs,cycles);
% % for i = 1:1:length(cycles)
% %     cycles(i).drawFullCycle(subGraphs(i)); 
% % end
% End Case 8




%% Case 9
% % clear all
% % close all
% % clc
% % graph = Graph(0); % create the graph object : index 0 for main graph 
% % rand('seed', 7); %% 4,5,7 seed
% % numOfTargets = 12; % 15 % 12
% % dimentionOfSpace = 2;
% % sizeOfSpace = 1;  % 1x1 square
% % communicationRadius = 0.5; %0.7 %0.5 % 0.7
% % 
% % graph = graph.loadARandomGraph(numOfTargets, dimentionOfSpace, sizeOfSpace, communicationRadius)
% % graph.drawGraph(1);
% % numOfAgents = 3;
% % graph = graph.loadShortestPaths2('fast'); % try 'slow' if have patience 
% % partitions = graph.partitionTheGraph(numOfAgents,'smallest lower bound'); % spectral clustering approach 2
% % 
% % subGraphs = graph.generateSubGraphs(partitions);
% % cycles = [];
% % for subGraph = subGraphs
% %     
% %     cycleTSP = Cycle(subGraph.findTSPCycle()); % Cycle([8,1,10,8,1,10,8,15,13]);
% %     cycleTSP = cycleTSP.updateEdgeList(subGraph);
% %     
% %     greedyExpandedCycle = cycleTSP.executeGreedyExpansions(subGraph);
% %     cycles = [cycles, greedyExpandedCycle];  
% % 
% % end
% % 
% % graph.drawAll(subGraphs,cycles);
% % for i = 1:1:length(cycles)
% %     cycles(i).drawFullCycle(subGraphs(i)); 
% % end
% % 
% % [subGraphs,cycles,gainVal] = graph.balanceThePartitions2(subGraphs,cycles);
% % while gainVal > 0
% %     graph.drawAll(subGraphs,cycles); % drawing intermediate subgraphs and cycles
% %     [subGraphs,cycles,gainVal] = graph.balanceThePartitions2(subGraphs,cycles);
% % end
% % 
% % for i = 1:1:length(cycles)
% %     cycles(i).drawFullCycle(subGraphs(i)); 
% % end
% [16.3, 6.9, 9.7] ---> [11.2,11.6,10.9]
% End Case 9




%% Case 9
% % clear all
% % close all
% % clc
% % graph = Graph(0); % create the graph object : index 0 for main graph 
% % rand('seed', 27); %% 4,5,7 seed
% % numOfTargets = 12; % 15 % 12
% % dimentionOfSpace = 2;
% % sizeOfSpace = 1;  % 1x1 square
% % communicationRadius = 0.5; %0.7 %0.5 % 0.7
% % 
% % graph = graph.loadARandomGraph(numOfTargets, dimentionOfSpace, sizeOfSpace, communicationRadius)
% % graph.drawGraph(1);
% % numOfAgents = 3;
% % graph = graph.loadShortestPaths2('fast'); % try 'slow' if have patience 
% % partitions = graph.partitionTheGraph(numOfAgents,'smallest lower bound'); % spectral clustering approach 2
% % 
% % subGraphs = graph.generateSubGraphs(partitions);
% % cycles = [];
% % for subGraph = subGraphs
% %     
% %     cycleTSP = Cycle(subGraph.findTSPCycle()); % Cycle([8,1,10,8,1,10,8,15,13]);
% %     cycleTSP = cycleTSP.updateEdgeList(subGraph);
% %     
% %     greedyExpandedCycle = cycleTSP.executeGreedyExpansions(subGraph);
% %     cycles = [cycles, greedyExpandedCycle];  
% % 
% % end
% % 
% % graph.drawAll(subGraphs,cycles);
% % for i = 1:1:length(cycles)
% %     cycles(i).drawFullCycle(subGraphs(i)); 
% % end
% % 
% % [subGraphs,cycles,gainVal] = graph.balanceThePartitions2(subGraphs,cycles);
% % while gainVal > 0
% %     graph.drawAll(subGraphs,cycles); % drawing intermediate subgraphs and cycles
% %     [subGraphs,cycles,gainVal] = graph.balanceThePartitions2(subGraphs,cycles);
% % end
% % 
% % for i = 1:1:length(cycles)
% %     cycles(i).drawFullCycle(subGraphs(i)); 
% % end
% [13.0, 4.4, 5.0] ---> [8.6, 8.0, 7.6]
% End Case 9




%% Case 10
% % clear all
% % close all
% % clc
% % graph = Graph(0); % create the graph object : index 0 for main graph 
% % rand('seed', 27); %% 4,5,7 seed
% % numOfTargets = 12; % 15 % 12
% % dimentionOfSpace = 2;
% % sizeOfSpace = 1;  % 1x1 square
% % communicationRadius = 0.5; %0.7 %0.5 % 0.7
% % 
% % graph = graph.loadARandomGraph(numOfTargets, dimentionOfSpace, sizeOfSpace, communicationRadius)
% % graph.drawGraph(1);
% % numOfAgents = 3;
% % 
% % % Option 1: Clustering based on 'smart' similarity values
% % graph = graph.loadShortestPaths2('fast'); % try 'slow' if have patience 
% % partitions = graph.partitionTheGraph(numOfAgents,'smallest lower bound'); % spectral clustering approach 2
% % 
% % subGraphs = graph.generateSubGraphs(partitions);
% % cycles = [];
% % for subGraph = subGraphs
% %     
% %     cycleTSP = Cycle(subGraph.findTSPCycle()); % Cycle([8,1,10,8,1,10,8,15,13]);
% %     cycleTSP = cycleTSP.updateEdgeList(subGraph);
% %     
% %     greedyExpandedCycle = cycleTSP.executeGreedyExpansions(subGraph);
% %     cycles = [cycles, greedyExpandedCycle];  
% % 
% % end
% % %drawing the initial subgraphs and cycles
% % graph.drawAll(subGraphs,cycles);
% % for i = 1:1:length(cycles)
% %     cycles(i).drawFullCycle(subGraphs(i)); 
% % end
% % [subGraphs,cycles,gainVal] = graph.balanceThePartitions2(subGraphs,cycles);
% % while gainVal > 0
% %     graph.drawAll(subGraphs,cycles); % drawing intermediate subgraphs and cycles
% %     [subGraphs,cycles,gainVal] = graph.balanceThePartitions2(subGraphs,cycles);
% % end
% % 
% % for i = 1:1:length(cycles)
% %     cycles(i).drawFullCycle(subGraphs(i)); 
% % end
% End Case 10





%% Sam's main code : Very slow
% % clear all
% % close all
% % clc
% % graph = Graph(0); % create the graph object : index 0 for main graph 
% % rand('seed', 7); %% 4,5,7 seed
% % numOfTargets = 7; % 15 % 12
% % dimentionOfSpace = 2;
% % sizeOfSpace = 1;  % 1x1 square
% % communicationRadius = 0.7; %0.7 %0.5 % 0.7
% % 
% % % % partition the graph
% % numOfAgents = 2;
% % 
% % graph = graph.loadARandomGraph(numOfTargets, numOfAgents, dimentionOfSpace, sizeOfSpace, communicationRadius);
% % graph.drawGraph(1);
% % 
% % 
% % graph = graph.loadShortestPaths2('fast'); % try 'slow' if have patience 
% % partitions = graph.partitionTheGraph(numOfAgents,'smallest lower bound'); % spectral clustering approach 2
% % 
% % subGraphs = graph.generateSubGraphs(partitions);
% % cycles = [];
% % for subGraph = subGraphs
% %     
% %     cycleTSP = Cycle(subGraph.findTSPCycle()); % Cycle([8,1,10,8,1,10,8,15,13]);
% %     cycleTSP = cycleTSP.updateEdgeList(subGraph);
% %     cycleTSP = cycleTSP.createCycleEvaluator(subGraph);
% %     cycleTSP.searchWithLB = 0;
% %     
% %     greedyExpandedCycle = cycleTSP.executeGreedyExpansions(subGraph);
% %     cycles = [cycles, greedyExpandedCycle];  
% % 
% % end
% % %drawing the initial subgraphs and cycles
% % graph.drawAll(subGraphs,cycles);
% % for i = 1:1:length(cycles)
% %     cycles(i).drawFullCycle(subGraphs(i)); 
% % end
% % % end Sam's main code: Very slow




%% Sam's main code : Fast, with target exchanges
clear all
close all
clc
graph = Graph(0); % create the graph object : index 0 for main graph 
rand('seed', 7); %% try this with 2 to see the error
numOfTargets = 7; % 
dimentionOfSpace = 2;
sizeOfSpace = 1;  % 1x1 square
communicationRadius = 0.7; 

% % partition the graph
numOfAgents = 2;

graph = graph.loadARandomGraph(numOfTargets, numOfAgents, dimentionOfSpace, sizeOfSpace, communicationRadius);
graph.drawGraph(1);


graph = graph.loadShortestPaths2('fast'); % try 'slow' if have patience 
partitions = graph.partitionTheGraph(numOfAgents,'smallest lower bound'); % spectral clustering approach 2

subGraphs = graph.generateSubGraphs(partitions);
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
graph.drawAll(subGraphs,cycles);
costs1 = [];
lowerBounds1 = [];
for i = 1:1:length(cycles)
    cycles(i).drawFullCycle(subGraphs(i)); 
    cost = cycles(i).getCost();
    costs1 = [costs1, cost];
    lowerBound = cycles(i).getLowerBound(subGraphs(i));
    lowerBounds1 = [lowerBounds1, lowerBound];
    
    dwellTimes{i} = cycles(i).cycleEvaluator.getDwellTimes();
end
%interchange algorithm
[subGraphs,cycles,gainVal] = graph.balanceThePartitions2(subGraphs,cycles);
while gainVal > 0
    graph.drawAll(subGraphs,cycles); % drawing intermediate subgraphs and cycles
    [subGraphs,cycles,gainVal] = graph.balanceThePartitions2(subGraphs,cycles);
end
%drawing the fianl subgraphs and cycles with costs, lowerBounds and dwell times
costs2 = [];
lowerBounds2 = [];
dwellTimes = cell(numOfAgents,1);
for i = 1:1:numOfAgents
    cycles(i).drawFullCycle(subGraphs(i)); 
    
    cost = cycles(i).getCost();
    costs2 = [costs2, cost];
    lowerBound = cycles(i).getLowerBound(subGraphs(i));
    lowerBounds2 = [lowerBounds2, lowerBound];
    
    dwellTimes{i} = cycles(i).cycleEvaluator.getDwellTimes();
end
costs1
lowerBounds1
costs2
lowerBounds2
dwellTimes{:}
% end Sam's main code: Fast, with target exchanges










%% Old Codes
% % graph = Graph(0); % create the graph object : index 0 for main graph 
% % rand('seed', 7); %% 4,5,7 seed
% % numOfTargets = 7; % 15 % 12
% % dimentionOfSpace = 2;
% % sizeOfSpace = 1;  % 1x1 square
% % communicationRadius = 0.7; %0.7 %0.5 % 0.7
% % 
% % % % partition the graph
% % numOfAgents = 2;
% % 
% % graph = graph.loadARandomGraph(numOfTargets, numOfAgents, dimentionOfSpace, sizeOfSpace, communicationRadius)
% % graph.drawGraph(1);
% % 
% % 
% % 
% % 
% % % Option 1: Clustering based on 'smart' similarity values
% % graph = graph.loadShortestPaths2('fast'); % try 'slow' if have patience 
% % partitions = graph.partitionTheGraph(numOfAgents,'smallest lower bound'); % spectral clustering approach 2
% % 
% % % Option 2: Clustering based on shortest distance based similarity values
% % % partitions = graph.partitionTheGraph(numOfAgents,'shortest path'); % spectral clustering approach 1
% % 
% % % Option 3: Clustering manually
% % % partitions = {[1,8,10,13,15],[3,4,5,7,9],[2,6,11,12,14]}; % Manual - visually Balanced
% % % partitions = {[1,8,13,15],[3,4,5,7,9,10],[2,6,11,12,14]}; % Manual - visually unbalanced
% % 
% % subGraphs = graph.generateSubGraphs(partitions);
% % cycles = [];
% % for subGraph = subGraphs
% %     
% %     cycleTSP = Cycle(subGraph.findTSPCycle()); % Cycle([8,1,10,8,1,10,8,15,13]);
% %     cycleTSP = cycleTSP.updateEdgeList(subGraph);
% %     
% %     greedyExpandedCycle = cycleTSP.executeGreedyExpansions(subGraph);
% %     cycles = [cycles, greedyExpandedCycle];  
% % 
% % end
% % %drawing the initial subgraphs and cycles
% % graph.drawAll(subGraphs,cycles);
% % for i = 1:1:length(cycles)
% %     cycles(i).drawFullCycle(subGraphs(i)); 
% % end
% % % % end
% % 
% % % % interchange algorithm
% % % [balancedSubGraphs, balancedCycles] = graph.balanceThePartitions1(subGraphs,cycles); % sell only the max target: Does not seem to work 
% % [subGraphs,cycles,gainVal] = graph.balanceThePartitions2(subGraphs,cycles);
% % while gainVal > 0
% %     graph.drawAll(subGraphs,cycles); % drawing intermediate subgraphs and cycles
% %     [subGraphs,cycles,gainVal] = graph.balanceThePartitions2(subGraphs,cycles);
% % end
% % % drawing the final cycles
% % for i = 1:1:length(cycles)
% %     cycles(i).drawFullCycle(subGraphs(i)); 
% % end
% % % % end iterchange
%% end Old code


% % old codes for cycle expansion
% % Testing TCEOs an distance matrix 2
% cycleTemp = Cycle([5]) % 5,2,4,1,4,2
% cycleTemp = cycleTemp.updateEdgeList(graph);
% cycleTemp.drawCycle(graph,'r');
% cycleTemp.drawFullCycle(graph);
% 
% [gainValue, expandedCycle] = expandCycleToInclude3(obj,graph,externalTarget)

% [dis, cyc, endCyc] = graph.constructACycleStartingFromANode(1,'fast');
% graph.drawGraph(count+1);
% endCyc.drawCycle(graph,'r');
% endCyc.drawFullCycle(graph)
% % End

% Old code for 1 - agent cases:

% % Uncomment this block to use a temporary cycle (testing: only when numOfTargets = 5) 
% cycleTemp = Cycle([5 2 4 1 3 2 4 3])
% cycleTemp = cycleTemp.updateEdgeList(graph);
% % cycleTemp.drawCycle(graph,'r');
% cycleTemp.drawFullCycle(graph);
% % end

% % uncommment this block to see the TSP cycle
% cycleTSP = Cycle(graph.findTSPCycle());
% cycleTSP = cycleTSP.updateEdgeList(graph);
% % cycleTSP.drawCycle(graph,'y');
% cycleTSP.drawFullCycle(graph);
% cycleTemp = cycleTSP;
% % end


% % Lowebound 
% [lowerBound,i] = cycleTemp.getLowerBound(graph)
% % end


% % expanding method 0: step 1
% [gainValue, expandedCycle] = cycleTemp.greedy1OPT(graph);
% expandedCycle.drawFullCycle(graph);
% cycleTemp = expandedCycle

% [gainValue, expandedCycle] = cycleTemp.greedy1OPT(graph);
% expandedCycle.drawFullCycle(graph);
% cycleTemp = expandedCycle
% % end

% % expanding method 1: step 1
% [gainValue, expandedCycle_11] = cycleTemp.greedyExpandCycle1(graph)
% expandedCycle_11.drawFullCycle(graph)
% 
% % expanding method 1: step 2
% [gainValue, expandedCycle_12] = expandedCycle_11.greedyExpandCycle1(graph)
% expandedCycle_12.drawFullCycle(graph)
% % end


% % expanding method 2: step 1
% cycleTemp = expandedCycle_11;
% [gainValue, expandedCycle_21] = cycleTemp.greedyExpandCycle2(graph)
% expandedCycle_21.drawFullCycle(graph)
% 
% % expanding method 2: step 2
% [gainValue, expandedCycle_22] = expandedCycle_21.greedyExpandCycle2(graph)
% expandedCycle_22.drawFullCycle(graph)
% % end


% % Automated greedy expansion
% currentTargetCycle = cycleTemp.executeGreedyExpansions(graph)
% % end