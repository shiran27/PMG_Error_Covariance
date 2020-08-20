classdef Cycle
    % Class for a generic cycle
    
    properties
        targetList
        edgeList
        transitTimeList
        
        auxTargetIndexList
        subCycleMatrix
        
        searchWithLB
        cycleEvaluator
        conversionGraphIdToCycleTargetId
    end
    
    methods
        
        function obj = Cycle(targetList)
            % Constructor
            obj.targetList = targetList;
            obj.searchWithLB = 1;
        end
        
        function obj = createCycleEvaluator(obj,graph)
            indexTarget = 1;
            obj.conversionGraphIdToCycleTargetId = zeros(max(obj.targetList),1);
            targetsSS = cell(length(obj.targetList),1);
            
            for indexTarget = 1:length(graph.targets)
                %if(obj.auxTargetIndexList(indexList) == 1)
                    targetsSS{indexTarget} = TargetSS(graph.targets(indexTarget).A,graph.targets(indexTarget).G,graph.targets(indexTarget).G,[]);
                    obj.conversionGraphIdToCycleTargetId(graph.targets(indexTarget).index) = indexTarget;
                    %indexTarget = indexTarget+1;
                %end
            end
            obj.cycleEvaluator = TrajectoryHandler(targetsSS,1,1);
        end
        
        function obj = updateEdgeList(obj, graph)
            
            % Update the list of edges based on the list of targets on the
            % given graph
            N = length(obj.targetList);
            obj.edgeList = zeros(1,N);
            obj.transitTimeList = zeros(1,N);
            if N>1
                for i = 1:1:(N-1)
                    edgeIndex = graph.getEdgeIndex(obj.targetList(i:i+1)); % edge index relative to the "graph"
                    obj.edgeList(i) = edgeIndex; 
                    obj.transitTimeList(i) = graph.edges(edgeIndex).getLength();
                end
                edgeIndex = graph.getEdgeIndex([obj.targetList(end), obj.targetList(1)]);
                obj.edgeList(end) = edgeIndex; % edge index relative to the "graph"
                obj.transitTimeList(end) = graph.edges(edgeIndex).getLength();
            else
                obj.edgeList = [];
                obj.transitTimeList = 0;
            end
            
            
            % update the sub-cycle data
            obj.auxTargetIndexList = zeros(1,N);
            obj.subCycleMatrix = zeros(N,N);
            for relativeTargetID = 1:1:length(graph.targets) % search for each target occurance
                i = graph.targets(relativeTargetID).index;
                occurances = find(obj.targetList == i);
                
                if isempty(occurances) % neglected target
                    continue;
                end
                
                auxTargetCount = 1;
                subCycleVectorSum = zeros(N,1);
                for ind = occurances
                    obj.auxTargetIndexList(ind) = auxTargetCount;
                    if auxTargetCount > 1
                        startInd = occurances(auxTargetCount-1)+1;
                        endInd = ind;
                        subCycleVector = zeros(N,1);
                        subCycleVector(startInd:endInd) = ones(endInd-startInd+1,1);
                                                
                        obj.subCycleMatrix(:,ind) = subCycleVector;
                        subCycleVectorSum = subCycleVectorSum + subCycleVector;
                    end
                    auxTargetCount = auxTargetCount + 1;  
                end
                obj.subCycleMatrix(:,occurances(1)) = subCycleVectorSum == 0;
            end
            
        end
        
        
        function outputArg = printCycle(obj)
            
            obj.targetList  % list of target indeces
            
            edgeArray = [];
            for i=1:1:length(obj.edgeList)
                edgeArray(i) = obj.edgeList(i).targets;
            end
            edgeArray'      % list of target pairs (indeces) corresponding to edges
            
            obj.edgeList    % list of edge indeces
            
        end
        
        function outputArg = drawCycle(obj, graph, color)
            hold on
            for i = obj.edgeList
                graph.edges(i).highlightEdge(color);
            end
            sizeOfSpace = 1;
            axis([0,sizeOfSpace,0,sizeOfSpace]);
            % highlihht the associated edges
        end
        
        
        
        function outputArg = drawFullCycle(obj, graph)
            r = 0.3;
            center = [0.5,0.5];
            
            ratio = 2*pi/sum(obj.transitTimeList);
            N = length(obj.targetList);
            theta = zeros(N,1);
            for i = 2:1:N
                theta(i) = sum(obj.transitTimeList(1:i-1))*ratio;
            end
            pos = center + [r*cos(theta), r*sin(theta)]; % exploiting dimension fluidity of sum operator
            
            figure
            hold on; grid on;
            viscircles([0.5,0.5],r,'Color','r');
            plot(pos(:,1),pos(:,2),'*b');
            axis equal;  axis([0,1,0,1]);
            
            r= r + 0.05; posOld = pos;
            pos = center + [r*cos(theta), r*sin(theta)];
            l_bData = [];
            textString1 = ['\Xi = \{ '];
            for i = 1:1:N
                textString1i = [num2str(obj.targetList(i)),'^',num2str(obj.auxTargetIndexList(i))];
                textString1 = [textString1, textString1i, ','];
                text(pos(i,1),pos(i,2),textString1i,'Color','blue','FontSize',10);
                % lowebound
                subCycleTime = obj.transitTimeList*obj.subCycleMatrix(:,i);
                
                lb_i = graph.getTarget(obj.targetList(i),'target').computeLowerBound(subCycleTime);
                thickness = 0.01; scalingFactor = 500;
                rectangle('Position',[posOld(i,1)-thickness, posOld(i,2), 2*thickness, lb_i/scalingFactor],...
                    'FaceColor',[0.5 .5 0 0.5],'EdgeColor','b','LineWidth',0.01)
                lbData(i) = lb_i;
            end
            [val,i] = max(lbData); 
            textString1 = [textString1(1:end-1),' \}'];
            textString2 = ['$\hat{J}(\bar{\Xi})=',num2str(round(val,3)),' $ (Crit. Target: $',...
                num2str(obj.targetList(i)),'^',num2str(obj.auxTargetIndexList(i)),' )$'];
            text(0.01,0.05,textString1,'Color','blue','FontSize',10);
            text(0.01,0.95,textString2,'Color','blue','FontSize',10,'Interpreter','Latex');
%             title(['Cycle View: Max Lower Bound: ',num2str(val),' at Target ',...
%                 num2str(obj.targetList(i)),'^',num2str(obj.auxTargetIndexList(i))])
            
            % plot not included targets
            count = 1;
            ratio = 0.5/length(graph.targets);
            for target = graph.targets
                if sum(target.index == obj.targetList)==0
                    plot(0.9,0.1+ratio*count,'*b');
                    text(0.95,0.1+ratio*count,num2str(target.index),'Color','blue','FontSize',10);
                    count = count + 1;
                end
            end
        end
        
        
        % more efficient way to Lb
        function [lb, lb_ij] = getLowerBound(obj, graph)
            lb = 0;
            for relativeTargetID = 1:length(graph.targets)
                i = graph.targets(relativeTargetID).index;
                occurances = find(obj.targetList == i);
                
                if isempty(occurances)
                    continue;
                end
                
                subCycleTimes = obj.transitTimeList*obj.subCycleMatrix(:,occurances);
                lb_i = graph.targets(relativeTargetID).computeLowerBound(subCycleTimes);
                if lb_i >= lb
                    lb = lb_i;
                    [~,ind] = max(subCycleTimes);
                    lb_ij = occurances(ind);
                end
            end
            lb_ij = [obj.targetList(lb_ij), obj.auxTargetIndexList(lb_ij), lb_ij]; % target index, aux target index, plcaemnt in the cycle
        end
        
        %% operations with complete target-cycles
        % 1 - OPT method
        function [gainValue, expandedCycle] = greedy1OPT(obj,graph)
            
           [lb,ij] =  obj.getLowerBound(graph);
           if obj.searchWithLB==0
                lb = obj.getCost();
           end
           indexInCycle = ij(3);
           
           subCycleVector = obj.subCycleMatrix(:,indexInCycle); % we need to break this subcycle into two
           subCycleIndices = find(subCycleVector==1); % ex:: if this is [5,6,7,8] we only can add new aux target of i at 6
           if length(subCycleIndices)>3 % otherwise no way to remove an edge within the subcycle
               maxGain = 0;
               maxGainCycle = obj;
               for ind = subCycleIndices'
                   
                   prevInd = ind - 1;
                   if prevInd < 1
                       prevInd2 = graph.getTarget(obj.targetList(end),'id');
                       ind2 = graph.getTarget(obj.targetList(ind),'id');
                       path = graph.shortestPathMatrix{prevInd2,ind2};
                       if length(path)==1
                            path = graph.getHandicappedShortestPath(obj.targetList(end),obj.targetList(ind));
                       end
                       newTargetList = [path(1:end-1), obj.targetList(1:end)];
                   else
                       prevInd2 = graph.getTarget(obj.targetList(ind-1),'id');
                       ind2 = graph.getTarget(obj.targetList(ind),'id');
                       path = graph.shortestPathMatrix{prevInd2,ind2};
                       if length(path)==1
                            path = graph.getHandicappedShortestPath(obj.targetList(ind-1),obj.targetList(ind));
                       end
                       newTargetList = [obj.targetList(1:prevInd), path(1:end-1), obj.targetList(ind:end)];
                   end
                   
                   
                    cycleTemp = Cycle(newTargetList);
                    cycleTemp = cycleTemp.updateEdgeList(graph);
                    cycleTemp = cycleTemp.createCycleEvaluator(graph);
                    cycleTemp.searchWithLB = obj.searchWithLB;
                    if obj.searchWithLB==0
                        lbTemp = obj.getCost();
                    else
                        [lbTemp,~] = cycleTemp.getLowerBound(graph);
                    end
                    gain = lb - lbTemp; % higher the better
                    if gain > maxGain
                        maxGain = gain;
                        maxGainCycle = cycleTemp;
                    end
                   
               end
               
               gainValue = maxGain;
               expandedCycle = maxGainCycle;
           else
               gainValue = 0;
               expandedCycle = obj;
           end
        
        end
            
            
        % Greedy expansion method 1
        function [gainValue, expandedCycle] = greedyExpandCycle1(obj,graph)
           
            [lb,ij] =  obj.getLowerBound(graph);
            if obj.searchWithLB==0
                lb = obj.getCost();
            end
            indexInCycle = ij(3);
           
           subCycleVector = obj.subCycleMatrix(:,indexInCycle); % we need to break this subcycle into two
           subCycleIndices = find(subCycleVector==1); % ex:: if this is [5,6,7,8] we only can add new aux target of i at 6
           if length(subCycleIndices)>3 % otherwise no way to remove an edge within the subcycle
               maxGain = 0;
               maxGainCycle = obj;
               for ind = subCycleIndices'
                   
                   nextInd = ind + 1;
                   prevInd = ind - 1;
                   if nextInd > length(obj.targetList)
                       nextInd = 1;
                       newTargetList = [obj.targetList, ij(1)];
                   elseif  prevInd < 1
                       prevInd = length(obj.targetList);
                       newTargetList = [obj.targetList(1:ind), ij(1), obj.targetList(nextInd:end)];
                   else
                       newTargetList = [obj.targetList(1:ind), ij(1), obj.targetList(nextInd:end)];
                   end
                   cond = sum(obj.targetList([prevInd,ind,nextInd]) == ij(1)) == 0;
                   
                   if cond
                        cycleTemp = Cycle(newTargetList);
                        cycleTemp = cycleTemp.updateEdgeList(graph);
                        cycleTemp = cycleTemp.createCycleEvaluator(graph);
                        cycleTemp.searchWithLB = obj.searchWithLB;
                        if obj.searchWithLB==0
                            lbTemp = obj.getCost();
                        else
                            [lbTemp,~] = cycleTemp.getLowerBound(graph);
                        end
                        gain = lb - lbTemp; % higher the better
                        if gain > maxGain
                            maxGain = gain;
                            maxGainCycle = cycleTemp;
                        end
                   end
                   
               end
               
               gainValue = maxGain;
               expandedCycle = maxGainCycle;
           else
               gainValue = 0;
               expandedCycle = obj;
           end
           
           
        end
        
        
        % Greedy expansion method 2: 
        function [gainValue, expandedCycle] = greedyExpandCycle2(obj,graph)
           [lb,ij] =  obj.getLowerBound(graph);
           if obj.searchWithLB==0
                lb = obj.getCost();
            end
            
           indexInCycle = ij(3);
           
           subCycleVector = obj.subCycleMatrix(:,indexInCycle); % we need to break this subcycle into two
           subCycleIndices = find(subCycleVector==1); 
           
           maxGain = 0;
           maxGainCycle = obj;
           for ind = subCycleIndices'

               if ind == length(obj.targetList)
                   newTargetList = [obj.targetList, ij(1), obj.targetList(end)];
               else
                   newTargetList = [obj.targetList(1:ind), ij(1), obj.targetList(ind:end)];
               end
               cond = obj.targetList(ind) ~= ij(1);

               if cond
                    cycleTemp = Cycle(newTargetList);
                    cycleTemp = cycleTemp.updateEdgeList(graph);
                    cycleTemp = cycleTemp.createCycleEvaluator(graph);
                    cycleTemp.searchWithLB = obj.searchWithLB;
                    if obj.searchWithLB==0
                        lbTemp = obj.getCost();
                    else
                        [lbTemp,~] = cycleTemp.getLowerBound(graph);
                    end
            
                    gain = lb - lbTemp; % higher the better
                    if gain > maxGain
                        maxGain = gain;
                        maxGainCycle = cycleTemp;
                    end
               end

           end

           gainValue = maxGain;
           expandedCycle = maxGainCycle;
           
        end
        
        
        % Iterative Greedy Expansion Algorithm
        function  currentTargetCycle = executeGreedyExpansions(obj,graph)
            
            gain = 1;
            currentTargetCycle = obj;
            
            % 1-OPT
            while gain > 0
                [gainValue, expandedCycle] = currentTargetCycle.greedy1OPT(graph);
                if gainValue > 0
                    currentTargetCycle = expandedCycle;
                    gain = gainValue;
%                     currentTargetCycle.drawFullCycle(graph);
                else
                    gain = 0;
                end
            end
            
            % Further Expansions
            while gain > 0 
                [gainValue1, expandedCycle_1] = currentTargetCycle.greedyExpandCycle1(graph);
                [gainValue2, expandedCycle_2] = currentTargetCycle.greedyExpandCycle2(graph);
                if gainValue1 > gainValue2 && gainValue1 > 0
                    currentTargetCycle = expandedCycle_1;
                    gain = gainValue1;
%                     currentTargetCycle.drawFullCycle(graph);
                elseif gainValue2 >= gainValue1 && gainValue2 > 0
                    currentTargetCycle = expandedCycle_2;
                    gain = gainValue2;
%                     currentTargetCycle.drawFullCycle(graph);
                else
                    gain = 0;
                end
            end
            
        end
        
        
        function  currentTargetCycle = executeGreedyExpansions2(obj,graph)
            
            gain = 1;
            currentTargetCycle = obj;
            
            % Further Expansions
            while gain > 0 
                [gainValue1, expandedCycle_1] = currentTargetCycle.greedyExpandCycle1(graph);
                [gainValue2, expandedCycle_2] = currentTargetCycle.greedyExpandCycle2(graph);
                if gainValue1 > gainValue2 && gainValue1 > 0
                    currentTargetCycle = expandedCycle_1;
                    gain = gainValue1;
% %                     currentTargetCycle.drawFullCycle(graph);
                elseif gainValue2 >= gainValue1 && gainValue2 > 0
                    currentTargetCycle = expandedCycle_2;
                    gain = gainValue2;
% %                     currentTargetCycle.drawFullCycle(graph);
                else
                    gain = 0;
                end
            end
            
        end
        
        
        %% operations with incomplete target-cycles
        
        % Target Cycle Expansion Operation 1
        function [gainValue, expandedCycle] = expandCycleToInclude1(obj,graph,externalTarget)
           [lb,ij] =  obj.getLowerBound(graph);
           indexInCycle = ij(3);
           
           subCycleVector = obj.subCycleMatrix(:,indexInCycle); % we need to break this subcycle into two
           subCycleIndices = find(subCycleVector==1); % ex:: if this is [5,6,7,8] we only can add new aux target of i at 6
           if length(subCycleIndices)>=2 % otherwise no way to remove an edge within the subcycle
               maxGain = -inf;
               maxGainCycle = obj;
               for ind = subCycleIndices'
                   
                   nextInd = ind + 1;
                   prevInd = ind - 1;
                   if nextInd > length(obj.targetList)
                       nextInd = 1;
                       newTargetList = [obj.targetList, externalTarget];
                   elseif  prevInd < 1
                       prevInd = length(obj.targetList);
                       newTargetList = [obj.targetList(1:ind), externalTarget, obj.targetList(nextInd:end)];
                   else
                       newTargetList = [obj.targetList(1:ind), externalTarget, obj.targetList(nextInd:end)];
                   end
                   
                    cycleTemp = Cycle(newTargetList);
                    cycleTemp = cycleTemp.updateEdgeList(graph);
                    cycleTemp = cycleTemp.createCycleEvaluator(graph);
                    cycleTemp.searchWithLB = obj.searchWithLB;
                    if obj.searchWithLB==0
                        lbTemp = obj.getCost();
                    else
                        [lbTemp,~] = cycleTemp.getLowerBound(graph);
                    end
                    gain = lb - lbTemp; % higher the better
                    if gain > maxGain
                        maxGain = gain;
                        maxGainCycle = cycleTemp;
                    end
                   
               end
               
               gainValue = maxGain;
               expandedCycle = maxGainCycle;
               
           else
               gainValue = -inf;
               expandedCycle = obj;
           end
        end
        
        % Target Cycle Expansion Operation 2
        function [gainValue, expandedCycle] = expandCycleToInclude2(obj,graph,externalTarget)
           [lb,ij] =  obj.getLowerBound(graph);
           if obj.searchWithLB==0
                obj.getCost();
           end
           indexInCycle = ij(3);
           
           subCycleVector = obj.subCycleMatrix(:,indexInCycle); % we need to break this subcycle into two
           subCycleIndices = find(subCycleVector==1); 
           
           maxGain = -inf;
           maxGainCycle = obj;
           for ind = subCycleIndices'

               
               if length(obj.targetList)==1
                    newTargetList = [obj.targetList, externalTarget];
               else
                   if ind == length(obj.targetList)
                       newTargetList = [obj.targetList, externalTarget, obj.targetList(end)];
                   else
                       newTargetList = [obj.targetList(1:ind),  externalTarget, obj.targetList(ind:end)];
                   end
               end
               
               cycleTemp = Cycle(newTargetList);
               cycleTemp = cycleTemp.updateEdgeList(graph);
               cycleTemp = cycleTemp.createCycleEvaluator(graph);
               cycleTemp.searchWithLB = obj.searchWithLB;
               if obj.searchWithLB==0
                   lbTemp = obj.getCost();
               else
                   [lbTemp,~] = cycleTemp.getLowerBound(graph);
               end
               gain = lb - lbTemp; % higher the better
               if gain > maxGain
                   maxGain = gain;
                   maxGainCycle = cycleTemp;
               end
           end

           gainValue = maxGain;
           expandedCycle = maxGainCycle;
        end
        
        
        % Target Cycle Expansion Operation 3
        function [gainValue, expandedCycle] = expandCycleToInclude3(obj,graph,externalTarget)
            if obj.searchWithLB == 0
                lb = obj.getCost();
            else
                [lb,~] =  obj.getLowerBound(graph);
            end
            % connect the external target so that it cancels out some target
            % revisits
            N = length(obj.targetList);
            startEndNodes = [];
            subCycleMatrix2 = [obj.subCycleMatrix, obj.subCycleMatrix];
            targetList2 = [obj.targetList, obj.targetList];
            targetListList = {};
            for i = 1:1:N
                
                searchNode = i+1;
                while sum(subCycleMatrix2(:,searchNode)==0)>0
                    endNode = searchNode + 1;
                    if endNode > (N+1)
                        tempList = [targetList2((endNode-N-1):i), externalTarget];
                    else
                        tempList = [targetList2(1:i), externalTarget, targetList2(endNode:end-N)];
                    end
                    
                    faulty = false;
                    for t = obj.targetList
                        if sum(tempList==t)==0
                            faulty = true;
                        end
                    end
                    if faulty
                        break
                    end
                    targetListList{length(targetListList)+1} = tempList;
                    startEndNodes = [startEndNodes; targetList2(i), targetList2(endNode)];
                    if endNode > 2*N
                        break
                    end
                    searchNode = endNode;
                end
                
                
            end
            
            
           maxGain = -inf;
           maxGainCycle = obj;
           for i = 1:1:length(targetListList)
               newTargetList = targetListList{i};
               cycleTemp = Cycle(newTargetList);
               cycleTemp = cycleTemp.updateEdgeList(graph);
               cycleTemp = cycleTemp.createCycleEvaluator(graph);
               cycleTemp.searchWithLB = obj.searchWithLB;
               if obj.searchWithLB==0
                   lbTemp = obj.getCost();
               else
                   [lbTemp,~] = cycleTemp.getLowerBound(graph);
               end
               gain = lb - lbTemp; % higher the better
               if gain > maxGain
                   maxGain = gain;
                   maxGainCycle = cycleTemp;
               end
           end

           gainValue = maxGain;
           expandedCycle = maxGainCycle;
            
        end
        
        % Target Cycle Expansion Operation Combined
        function expandedCycle = expandCycleToInclude(obj,graph,externalTarget,method)
            % combine the above three functions with the refinements
            % to get the best way to expand the existing target-cycle to inclue the external target 
            slow = isequal(method,'slow');
            if obj.searchWithLB == 0
                lb = obj.getCost();
            else
                [lb, ~] = obj.getLowerBound(graph);
            end
            cycles = [];
            gains = [];
            
            [~, expandedCycle_1] = obj.expandCycleToInclude1(graph,externalTarget);
            if sum(expandedCycle_1.targetList==externalTarget)>0
                if slow
                    expandedCycle_1 = expandedCycle_1.executeGreedyExpansions2(graph);
                end
                if obj.searchWithLB == 0
                    lb1 = expandedCycle_1.getCost();
                else
                    [lb1, ~] = expandedCycle_1.getLowerBound(graph);
                end
                cycles = [cycles, expandedCycle_1];
                gains = [gains, lb-lb1];
            end
            
            [~, expandedCycle_2] = obj.expandCycleToInclude2(graph,externalTarget);
            if slow 
                expandedCycle_2 = expandedCycle_2.executeGreedyExpansions2(graph);
            end
            if obj.searchWithLB == 0
                lb2 = expandedCycle_2.getCost();
            else
                [lb2, ~] = expandedCycle_2.getLowerBound(graph);
            end
            cycles = [cycles, expandedCycle_2];
            gains = [gains, lb-lb2];
            
            [~, expandedCycle_3] = obj.expandCycleToInclude3(graph,externalTarget);
            if sum(expandedCycle_3.targetList==externalTarget)>0
                if slow
                    expandedCycle_3 = expandedCycle_3.executeGreedyExpansions2(graph);
                end
                if obj.searchWithLB == 0
                    lb3 = expandedCycle_3.getCost();
                else
                    [lb3, ~] = expandedCycle_3.getLowerBound(graph);
                end
                cycles = [cycles, expandedCycle_3];
                gains = [gains, lb-lb3];
            end
            
            [~,ind] = max(gains);
            expandedCycle = cycles(ind);
            
        end
        
        
        % contracting a cycle
        function cycle = contractCycleToExclude(obj, graph, internalTarget)
            
            contractedCycle = obj.targetList;
            occurances = 0;
            ind = 1;
            for i = obj.targetList
                if i==internalTarget
%                     contractedCycle(ind-occurances) = [];
                    
                    % replace by the shortest path!
                    current = contractedCycle(ind-occurances); % internalTarget or the sold Traget
                    if (ind-occurances-1) < 1
                        prev = contractedCycle(end);
                        next = contractedCycle(ind-occurances+1); 
                        shortestPath = graph.getHandicappedShortestPath2(prev,next,current);
                        contractedCycle = [shortestPath(1:end-1), contractedCycle((ind-occurances+1):end)];
                    elseif (ind-occurances+1) > length(contractedCycle) 
                        prev = contractedCycle(ind-occurances-1);
                        next = contractedCycle(1);
                        shortestPath = graph.getHandicappedShortestPath2(prev,next,current);
                        contractedCycle = [contractedCycle(1:(ind-occurances-1)), shortestPath(1:end-1)];
                    else
                        prev = contractedCycle(ind-occurances-1);
                        next = contractedCycle(ind-occurances+1); 
                        shortestPath = graph.getHandicappedShortestPath2(prev,next,current);
                        contractedCycle = [contractedCycle(1:(ind-occurances-1)), shortestPath(1:end-1), contractedCycle((ind-occurances+1):end)];
                    end
%                     prev
                    
%                     next
                 
                    
                    
                    % End replace by the shortest path!
                    occurances = occurances + 1;
                end
                ind = ind + 1;
            end
            
            % remove adjacent repetitions
            contractedCycleTemp = [];
            for i = 1:1:length(contractedCycle)
                if length(contractedCycle)==1
                    contractedCycleTemp = contractedCycle;
                    continue
                end
                
                if i==1
                    prev = contractedCycle(end);
                else
                    prev = contractedCycle(i-1);
                end
                
                if contractedCycle(i)~=prev
                    contractedCycleTemp = [contractedCycleTemp, contractedCycle(i)];
                end
            end
            contractedCycle = contractedCycleTemp;
            
            cycle = Cycle(contractedCycle); % Cycle([8,1,10,8,1,10,8,15,13]);
            cycle = cycle.updateEdgeList(graph);
            if length(cycle.targetList)>2
                cycle = cycle.executeGreedyExpansions(graph);
            end
        end
        
        
        %%
        % Sam
        function cost = getCost(obj)
                visitingSequence = obj.conversionGraphIdToCycleTargetId(obj.targetList);
                travelTime = obj.transitTimeList;
                tMin = 0.1*sum(travelTime);
                tMax = 3*sum(travelTime);
                obj.cycleEvaluator.visitingSequence = visitingSequence;
                obj.cycleEvaluator.travelTime = travelTime;
                [cost,~]=obj.cycleEvaluator.goldenRatioSearch(tMin,tMax);
        end
%         function lb = getLowerBound(obj,graph)
%             visitingSequence = obj.targetList;
%             travelTime = obj.transitTimeList;
%             % Sam
%             lb = 0;
%             for indexTarget = 1:length(graph.targets)
%                 currentTargetVisits = find(visitingSequence == indexTarget);
%                 currentTargetTravelTime = zeros(size(currentTargetVisits));
%                 for indexCurrentVisit = 1:length(currentTargetVisits)-1
%                     currentTargetTravelTime(indexCurrentVisit) = sum(travelTime(currentTargetVisits(indexCurrentVisit):currentTargetVisits(indexCurrentVisit+1)-1));
%                 end
%                 currentTargetTravelTime(end) = sum(travelTime(currentTargetVisits(end):end))+...
%                     sum(travelTime(1:currentTargetVisits(1)-1));
%                 currentLb = graph.targets(indexTarget).computeLowerBound(currentTargetTravelTime);
%                 if currentLb > lb
%                     lb = currentLb;
%                 end
%             end
%         end
        
        
    end
end

