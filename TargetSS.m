classdef TargetSS < handle
    properties
        A
        G
        Q
        R
        X
        tArrival
        tDeparture
        tTotal
        NdoublingSteps
        dtPlot
        peakTraces
        tOn
        tol
        position
        lambda
        Z
        S
    end
    methods
        function obj = TargetSS(A,G,Q,pos)
            obj.A = A;
            obj.G = G;
            obj.Q = Q;
            obj.NdoublingSteps = 5;     
            %obj.X = obj.updateCov(tArrival,tDeparture,tTotal);
            obj.dtPlot = 1e-3;
            obj.tol = 1e-6;
            obj.position = pos;
            
            if length(A)==1
                obj.R = 1/obj.G;
                obj.lambda = -sqrt(obj.R*obj.A^2+obj.Q)/sqrt(obj.R);
                
                z21 = obj.Q*sqrt(obj.R)/(2*sqrt(obj.R*obj.A^2+obj.Q));
                z22 = 0.5*(sqrt(obj.R)*obj.A/sqrt(obj.R*obj.A^2+obj.Q)+1);
                z11 = -0.5*obj.Q*sqrt(obj.R)/sqrt(obj.R*obj.A^2+obj.Q);
                z12 = 0.5-0.5*obj.A*sqrt(obj.R)/sqrt(obj.R*obj.A^2+obj.Q);
                obj.Z = [z11,z12;z21,z22];
                
                s11 = (-(obj.A*sqrt(obj.R))-sqrt(obj.Q+obj.A^2*obj.R))/(obj.Q*sqrt(obj.R));
                s12 = (-(obj.A*sqrt(obj.R))+sqrt(obj.Q+obj.A^2*obj.R))/(obj.Q*sqrt(obj.R));
                s21 = 1;
                s22 = 1;
                obj.S = [s11,s12;s21,s22];
            end

        end
        function M = get1DOnMatrixExponential(obj,time)
            exponential = exp(obj.lambda*time);
            M = obj.S*diag([exponential,1/exponential])*obj.Z;
        end
        function M = get1DOffMatrixExponential(obj,time)
            exponential = exp(obj.A*time);
            exp21 = (exponential-1/exponential)*obj.Q/(2*obj.A);
            M = [1/exponential,0;exp21,exponential];
        end
        function X = updateCov(obj,tArrival,tDeparture,tTotal)
            obj.tTotal = tTotal;
            Hon = [-obj.A,obj.G;obj.Q,obj.A'];
            Hoff = [-obj.A,zeros(size(obj.G));obj.Q,obj.A'];
            if iscolumn(tArrival)
                tArrival = tArrival';
            end
            if iscolumn(tDeparture)
                tDeparture = tDeparture';
            end
            obj.tArrival = tArrival;
            obj.tDeparture = tDeparture;
            obj.peakTraces = zeros(size(obj.tArrival));
            obj.tOn = obj.tDeparture-obj.tArrival;
            
            %Intercalating arrival and departure times
            eventTimes = [tArrival;tDeparture];
            eventTimes = eventTimes(:)';
            eventTimes = [0,eventTimes,obj.tTotal];
            
            eventDuration = diff(eventTimes);
            
            HOnePeriod = eye(size(Hon));
            H = cell(length(eventDuration),1);       
            
            for index = 1:length(eventDuration)
                if mod(index,2)==1
                    if isscalar(obj.A)
                        H{index} = obj.get1DOffMatrixExponential(eventDuration(index));
                        %HOnePeriod = obj.get1DOffMatrixExponential(eventDuration(index))*HOnePeriod;
                    else
                        H{index} = expm(eventDuration(index)*Hoff);
                        %HOnePeriod = expm(eventDuration(index)*Hoff)*HOnePeriod;
                    end
                else
                    if isscalar(obj.A)
                        H{index} = obj.get1DOnMatrixExponential(eventDuration(index));
                        %HOnePeriod = obj.get1DOnMatrixExponential(eventDuration(index))*HOnePeriod;
                    else
                        H{index} = expm(eventDuration(index)*Hon);
                        %HOnePeriod = expm(eventDuration(index)*Hon)*HOnePeriod;
                    end
                end
                HOnePeriod = H{index}*HOnePeriod;
            end
            Htotal = HOnePeriod;
            n = length(obj.A);
            
            %Doubling algorithm (structure preserving)
            Y11 = HOnePeriod(1:n,1:n);
            Y12 = HOnePeriod(1:n,n+1:2*n);
            Y21 = HOnePeriod(n+1:2*n,1:n);
            Y22 = HOnePeriod(n+1:2*n,n+1:2*n);
            
            %Decomposition in sympletic pairs
            Adoubling = inv(Y11);
            Gdoubling = Adoubling*Y12;
            Hdoubling = Y21*Adoubling;
            
            %Running the doubling algorithm
            for index=1:obj.NdoublingSteps
                Wdoubling = eye(length(obj.A))+Gdoubling*Hdoubling;
                Winv = inv(Wdoubling);
                V1 = Winv*Adoubling;
                V2 = Gdoubling*Winv';
                
                Hdoubling = Hdoubling+V1'*Hdoubling*Adoubling;
                Gdoubling = Gdoubling+Adoubling*V2*Adoubling';
                Adoubling = Adoubling*V1;
            end
            X = Hdoubling;
            
            obj.X = X;
            previousCov = X;
            for indexEvent = 1:length(eventDuration)-1
                if mod(indexEvent,2)==1
                    if isscalar(obj.A)
                        Hcurrent = obj.get1DOffMatrixExponential(eventDuration(indexEvent));
                    else
                        Hcurrent = expm(eventDuration(indexEvent)*Hoff);
                    end
                else
                    if isscalar(obj.A)
                        Hcurrent = obj.get1DOnMatrixExponential(eventDuration(indexEvent));
                    else
                        Hcurrent = expm(eventDuration(indexEvent)*Hon);
                    end
                end
%                 if mod(indexEvent,2)==0
%                     Hcurrent = expm(Hon*eventDuration(indexEvent));
%                 else
%                     Hcurrent = expm(Hoff*eventDuration(indexEvent));
%                 end
                Hcurrent = Hcurrent*[eye(size(obj.A));previousCov];
                V = Hcurrent(1:length(obj.A),:);
                U = Hcurrent(length(obj.A)+1:end,:);
                previousCov = U/V;
                if (mod(indexEvent,2))==1 
                    obj.peakTraces((indexEvent+1)/2)=trace(previousCov);
                end
            end
            
        end
        
        function plotTraceCovariance(obj,delay)
            Hon = [-obj.A,obj.G;obj.Q,obj.A'];
            Hoff = [-obj.A,zeros(size(obj.G));obj.Q,obj.A'];
            
            tArrival = obj.tArrival;
            tDeparture = obj.tDeparture;
            tTotal = obj.tTotal;
            
            %Intercalating arrival and departure times
            eventTimes = [tArrival;tDeparture];
            eventTimes = eventTimes(:)';
            eventTimes = [0,eventTimes,obj.tTotal];
            
            eventDuration = diff(eventTimes);
            
            times = [];
            traces = [];
            
            previousCov = obj.X;
            for indexEvent = 1:length(eventDuration)
                newTimes = (0:obj.dtPlot:eventDuration(indexEvent));
                if length(newTimes)>0
                    newTraces = zeros(size(newTimes));
                    for indexTimeStep = 1:length(newTimes)
                        if mod(indexEvent,2)==0
                            Hcurrent = expm(Hon*newTimes(indexTimeStep));
                        else
                            Hcurrent = expm(Hoff*newTimes(indexTimeStep));
                        end
                        Hcurrent = Hcurrent*[eye(size(obj.A));previousCov];
                        V = Hcurrent(1:length(obj.A),:);
                        U = Hcurrent(length(obj.A)+1:end,:);
                        cov = U*inv(V);
                        newTraces(indexTimeStep) = trace(cov);
                    end
                    if mod(indexEvent,2)==0
                        Hcurrent = expm(Hon*eventDuration(indexEvent));
                    else
                        Hcurrent = expm(Hoff*eventDuration(indexEvent));
                    end
                    Hcurrent = Hcurrent*[eye(size(obj.A));previousCov];
                    V = Hcurrent(1:length(obj.A),:);
                    U = Hcurrent(length(obj.A)+1:end,:);
                    previousCov = U*inv(V);
                    
                    times = [times,newTimes+eventTimes(indexEvent)];
                    traces = [traces,newTraces];
                end
                
            end
            figure(948);
            times = times+delay;
            indexes = find(times>obj.tTotal);
            times(indexes) = times(indexes)-obj.tTotal;
            [~,indexes] = sort(times);
            times = times(indexes);
            traces = traces(indexes);
            plot(times,traces,'LineWidth',3);
            xlabel('Time')
            ylabel('Trace');
            
        end
        
        function time = computeTimeFromPeak1D(obj,peak,tTravel)
            if ~isscalar(obj.A)
                disp('This is not an 1D system');
            end
            time = zeros(size(tTravel));
            
            zPrime = (obj.Z(2,1)+obj.Z(2,2)*peak)/(obj.Z(1,1)+obj.Z(1,2)*peak);
            
            for index=1:length(tTravel)
                auxUnderPeak = exp(-2*obj.A*tTravel(index));
                underPeak = peak*auxUnderPeak-obj.Q*(1-auxUnderPeak)/(2*obj.A);
                time(index) = -log((obj.S(1,1)*underPeak-obj.S(2,1))/(zPrime*obj.S(2,2)-zPrime*obj.S(1,2)*underPeak))/(2*obj.lambda);
            end
        end
        function [peak,time] = searchForOptimalTimeDistribution1D(obj,tOn,tTravel)
            if length(tTravel)>1
                %Computing the algebraic ricatti solution
                xss = (obj.A+sqrt(obj.A^2+obj.G*obj.Q))/obj.G;

                %Computing a lower bound for the peak
                maxTime = max(tTravel);
                peakMin = exp(2*obj.A*maxTime)*(xss+obj.Q*(1-exp(-2*obj.A*maxTime))/(2*obj.A));

                tMax = inf;

                %Computing an upper bound for the peak
                if isscalar(obj.A)
                    y = obj.get1DOnMatrixExponential(tOn);
                    y11 = y(2,1);
                    y12 = y(2,2);
                    y21 = y(1,1);
                    y22 = y(1,2);
                    tOff = sum(tTravel);
                    qtilda = obj.Q/(2*obj.A)*(exp(2*obj.A*tOff)-1);
                    a = y22;
                    b = y21-qtilda*y22-y12*exp(2*obj.A*tOff);
                    c = -y11-qtilda*y21*exp(2*obj.A*tOff);
                    peak = (-b+sqrt(b^2-4*a*c))/(2*a);
                    obj.X = peak;
                    peakMax = obj.X;
                    tMin = tOn;
                else
                    tVisit = ones(size(tTravel))*tOn/length(tTravel);
                    travelSoFar = [0,tTravel];
                    travelSoFar = travelSoFar(1:end-1);
                    tDeparture = cumsum([travelSoFar+tVisit]);
                    tArrival = tDeparture-tVisit;
                    obj.updateCov(tArrival,tDeparture,sum(tTravel+tVisit));
                    peakMax = max(obj.peakTraces);
                    tMin = sum(obj.computeTimeFromPeak1D(peakMax,tTravel));
                end
                while log(tMax/tMin)>obj.tol
                    avgPeak = sqrt(peakMax*peakMin);
                    tCurrentSplit = obj.computeTimeFromPeak1D(avgPeak,tTravel);
                    tCurrent = sum(tCurrentSplit);
                    if tCurrent>tOn
                        tMax = tCurrent;
                        peakMin = avgPeak;
                    else
                        tMin = tCurrent;
                        peakMax = avgPeak;
                    end   
                end
                peak = sqrt(peakMax*peakMin);
                if isempty(tCurrentSplit)
                    time = sqrt(tMin*tMax);
                else
                    time = tCurrentSplit;
                end
            else
                y = obj.get1DOnMatrixExponential(tOn);
                y11 = y(2,1);
                y12 = y(2,2);
                y21 = y(1,1);
                y22 = y(1,2);
                tOff = tTravel;
                qtilda = obj.Q/(2*obj.A)*(exp(2*obj.A*tOff)-1);
                a = y22;
                b = y21-qtilda*y22-y12*exp(2*obj.A*tOff);
                c = -y11-qtilda*y21*exp(2*obj.A*tOff);
                peak = (-b+sqrt(b^2-4*a*c))/(2*a);
                obj.X = peak;
                peak = obj.X;
                time = tOn;
                obj.updateCov(0,tOn,tOn+tTravel);
                peak = obj.X;
            end
            obj.tOn = time;
        end
        
        function updateTargetTime(obj,tTravel,tTotal)
            if iscolumn(tTravel)
                tTravel = tTravel';
            end
            travelSoFar = [0,tTravel];
            travelSoFar = travelSoFar(1:end-1);
            [peak,tVisit] = obj.searchForOptimalTimeDistribution1D(tTotal-sum(tTravel),tTravel);
            obj.tOn = tVisit;
            obj.X = peak;
            obj.tDeparture = cumsum([travelSoFar+tVisit]);
            obj.tArrival = obj.tDeparture-tVisit;
            obj.tTotal = sum(tTravel+tVisit);
        end
        
        function peakMin = computeLowerBound(obj,maxTime)
            maxTime = max(maxTime);
            xss = (obj.A+sqrt(obj.A^2+obj.G*obj.Q))/obj.G;
            peakMin = exp(2*obj.A*maxTime)*(xss+obj.Q*(1-exp(-2*obj.A*maxTime))/(2*obj.A));
        end
    end
end