classdef TrajectoryHandler < handle
    properties
        targets
        visitingSequence
        travelTime
        totalTime
        cost
        alpha
        beta
        Niter
        logTraces
        logTime
        goldenSearchTol
        %visitTime
        targetVisitTime
        logGoldenRatioPeriod
        logGoldenRatioCost
    end
    methods
        function obj = TrajectoryHandler(targets,visitingSequence,travelTime)
            obj.targets = targets;
            obj.visitingSequence = visitingSequence;
            obj.travelTime = travelTime;
            obj.alpha=1e-3;
            obj.Niter = 200;
            obj.goldenSearchTol = 0.01;
        end
        function [times,peak] = explorePeak(obj,tMin,tMax)
            times = linspace(tMin,tMax,100);
            peak = times;
            for index=1:length(times)
                peak(index) = obj.optimizeVisitingTimes(times(index));
            end
        end
        function avgPeak = optimizeVisitingTimes(obj,totalTime)
            travelTime = obj.travelTime;
            alpha = obj.alpha;
            visits = obj.visitingSequence;
            uniqueVisits = unique(obj.visitingSequence);
            Ntargets = length(uniqueVisits);
            
            obj.logTraces = zeros(max(uniqueVisits),obj.Niter);
            obj.logTime = zeros(max(uniqueVisits),obj.Niter);
            
            if length(obj.targetVisitTime)== 0
                targetVisitTime(uniqueVisits) = ones(1,(Ntargets))*totalTime/Ntargets;
            else
                targetVisitTime(uniqueVisits) = obj.targetVisitTime(uniqueVisits)/sum(obj.targetVisitTime(uniqueVisits))*totalTime;
            end
            
            targets = obj.targets;
            visitTime = zeros(size(travelTime));
                
            for indexAux = 1:Ntargets
                indexTarget = uniqueVisits(indexAux);
                currentTargetVisits = find(visits == indexTarget);
                for indexVisit = 1:length(currentTargetVisits)
                    visitTime(currentTargetVisits(indexVisit)) = targetVisitTime(indexTarget)/length(currentTargetVisits);
                end
            end
            
            newVisitTime = zeros(size(visitTime));
            avgPeak = 0;
            for indexAux = 1:Ntargets
                indexTarget = uniqueVisits(indexAux);
                currentTargetVisits = find(visits == indexTarget);
                if length(currentTargetVisits)==0
                    disp('Hi')
                end
                currentTargetTravelTime = zeros(size(currentTargetVisits));
                for indexCurrentVisit = 1:length(currentTargetVisits)-1
                    currentTargetTravelTime(indexCurrentVisit) =sum(travelTime(currentTargetVisits(indexCurrentVisit):currentTargetVisits(indexCurrentVisit+1)-1))+...
                        sum(visitTime(currentTargetVisits(indexCurrentVisit)+1:currentTargetVisits(indexCurrentVisit+1)-1));
                end
                currentTargetTravelTime(end) = sum(travelTime(currentTargetVisits(end):end))+...
                    sum(travelTime(1:currentTargetVisits(1)-1))+...
                    sum(visitTime(currentTargetVisits(end)+1:end))+...
                    sum(visitTime(1:currentTargetVisits(1)-1));
                targets{indexTarget}.updateTargetTime(currentTargetTravelTime,sum(currentTargetTravelTime)+targetVisitTime(indexTarget));
                newVisitTime(currentTargetVisits) = targets{indexTarget}.tOn;
                
                avgPeak = avgPeak + log(targets{indexTarget}.X)/Ntargets;
            end
            maxPeak = 10;
            minPeak = 1;
            indexIteration = 1;
            while (maxPeak-minPeak)/minPeak > 1e-4
                maxPeak = 0;
                minPeak = inf;
                visitTime = newVisitTime;
                newAvgPeak = 0;
                for indexAux=1:Ntargets
                    indexTarget = uniqueVisits(indexAux);
                    targetVisitTime(indexTarget) = targetVisitTime(indexTarget)+alpha*(log(targets{indexTarget}.X)-avgPeak);
                    currentTargetVisits = find(visits == indexTarget);
                    currentTargetTravelTime = zeros(size(currentTargetVisits));
                    for indexCurrentVisit = 1:length(currentTargetVisits)-1
                        currentTargetTravelTime(indexCurrentVisit) =sum(travelTime(currentTargetVisits(indexCurrentVisit):currentTargetVisits(indexCurrentVisit+1)-1))+...
                            sum(visitTime(currentTargetVisits(indexCurrentVisit)+1:currentTargetVisits(indexCurrentVisit+1)-1));
                    end
                    currentTargetTravelTime(end) = sum(travelTime(currentTargetVisits(end):end))+...
                        sum(travelTime(1:currentTargetVisits(1)-1))+...
                        sum(visitTime(currentTargetVisits(end)+1:end))+...
                        sum(visitTime(1:currentTargetVisits(1)-1));
                    targets{indexTarget}.updateTargetTime(currentTargetTravelTime,sum(currentTargetTravelTime)+targetVisitTime(indexTarget));
                    newVisitTime(currentTargetVisits) = (targets{indexTarget}.tOn);
                    obj.logTraces(indexTarget,indexIteration) = targets{indexTarget}.X;
                    if targets{indexTarget}.X < minPeak
                        minPeak = targets{indexTarget}.X;
                    end
                    if targets{indexTarget}.X > maxPeak
                        maxPeak = targets{indexTarget}.X;
                    end
                    obj.logTime(indexTarget,indexIteration) = sum(targets{indexTarget}.tOn);
                    newAvgPeak = newAvgPeak + log(targets{indexTarget}.X)/Ntargets;
                end
                avgPeak = newAvgPeak;
                indexIteration = indexIteration+1;
            end
            avgPeak = exp(avgPeak);
            obj.targetVisitTime = targetVisitTime;
        end
        function [costOpt,tOpt] = goldenRatioSearch(obj,tMin,tMax)
            goldenRatio = (sqrt(5)+1)/2;
            c = tMax - (tMax-tMin)/goldenRatio;
            d = tMin + (tMax-tMin)/goldenRatio;
            costC = obj.optimizeVisitingTimes(c);
            costD = obj.optimizeVisitingTimes(d);
            obj.logGoldenRatioCost = [costC,costD];
            obj.logGoldenRatioPeriod = [c,d];
            while abs(c-d)>obj.goldenSearchTol       
                if costC<costD
                    tMax = d;
                    costD = costC;
                    d = c;
                    c = tMax - (tMax-tMin)/goldenRatio;
                    costC = obj.optimizeVisitingTimes(c);
                    obj.logGoldenRatioCost = [obj.logGoldenRatioCost,costC];
                    obj.logGoldenRatioPeriod = [obj.logGoldenRatioPeriod,c];
                else
                    tMin = c;
                    costC = costD;
                    c=d;
                    d = tMin + (tMax-tMin)/goldenRatio;
                    costD = obj.optimizeVisitingTimes(d);
                    
                    obj.logGoldenRatioCost = [obj.logGoldenRatioCost,costD];
                    obj.logGoldenRatioPeriod = [obj.logGoldenRatioPeriod,d];
                end
            end
            costOpt = (costC+costD)/2;
            tOpt = (c+d)/2;
        end
        function lb = getLowerBound(obj,visitingSequence,travelTime)
            lb = 0;
            for indexTarget = 1:length(obj.targets)
                currentTargetVisits = find(visitingSequence == indexTarget);
                currentTargetTravelTime = zeros(size(currentTargetVisits));
                for indexCurrentVisit = 1:length(currentTargetVisits)-1
                    currentTargetTravelTime(indexCurrentVisit) =sum(travelTime(currentTargetVisits(indexCurrentVisit):currentTargetVisits(indexCurrentVisit+1)-1));
                end
                currentTargetTravelTime(end) = sum(travelTime(currentTargetVisits(end):end))+...
                    sum(travelTime(1:currentTargetVisits(1)-1));
                currentLb = obj.targets{indexTarget}.computeLowerBound(max(currentTargetTravelTime));
                if currentLb > lb
                    lb = currentLb;
                end
            end
        end
        
        
        function dwellTimes = getDwellTimes(obj)
            if length(obj.visitingSequence)==1
                dwellTimes = inf;
            else
                dwellTimes = zeros(length(obj.visitingSequence),1);
                indexTargetVisit = ones(length(obj.targets),1);
                for indexVisit = 1:length(dwellTimes)
                    currentTarget = obj.visitingSequence(indexVisit);
                    dwellTimes(indexVisit) = obj.targets{currentTarget}.tOn(indexTargetVisit(currentTarget));
                    indexTargetVisit(currentTarget) = indexTargetVisit(currentTarget)+1;
                end
            end
        end
        
    end
end