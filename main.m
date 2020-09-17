clear all
close all
clc


% create a random graph
graph = Graph(0); % create the graph object : index 0 for main graph 

%Case 1: 7, 7, 0.7, 2;   Case 2: 250, 7, 0.45, 2;   Case 3: 500, 10, 0.45, 4;    Case 4: 77, 10, 0.45, 4 
rand('seed', 500); %% 7, 250, 500
numOfTargets = 10; % 7, 10
dimentionOfSpace = 2;
sizeOfSpace = 1;  % 1x1 square
communicationRadius = 0.45; %0.7,  0.45

% load the agents
numOfAgents = 4; %2, 4

% Control Options: RHC-Re-Learn, Periodic, RHC-Periodic, RHC-Learn, RHC, BDC-Periodic, BDC, fixed, random
controlMethod = 'Periodic';
targetControllersEnabled = true;

% create the graph and draw
graph = graph.loadARandomGraph(numOfTargets, numOfAgents, dimentionOfSpace, sizeOfSpace, communicationRadius, targetControllersEnabled);
graph.drawGraph(1);

% if we try the periodic controller
if isequal(controlMethod,'RHC-Periodic') | isequal(controlMethod,'BDC-Periodic') | isequal(controlMethod,'Periodic')
    graph.assignCyclesToAgents(controlMethod);
end


timeResolution = 0.001;
periodT = 50; % 50, 200, 500
fixedHorizonValue = 2.5; %2.50727
cost = 0;
controlCost = 0;
textHandle1 = text();  textHandle2 = text(); textHandle3 = text(); 
plotMode = true; 
videoMode = true;
perturbMode = []; %[1 40; 0 2.1; 1 7; 0 7.1; inf 1]; % [1  2; 0  2.1; 1  7; 0  7.1]: activate at 2 and 7  for a duration of 0.1 s 
data = zeros(periodT/timeResolution+1, numOfTargets, 4);
timeSteps = 1;
frameCount = 1;
for t = 0:timeResolution:periodT
    
    % target update
    for i = 1:1:numOfTargets
        [costVal, dataVal, costVal2] = graph.targets(i).update(timeResolution,plotMode,t); 
        cost = cost + costVal;
        if targetControllersEnabled
            controlCost = controlCost + costVal2;
        end
        
        if plotMode
            data(timeSteps,i,:) = dataVal; %{phi,phiHat,Omega]
        end
    end
    
    % agent update
    horizon = min(fixedHorizonValue,periodT-t);
    for i = 1:1:numOfAgents
        graph.agents(i).update(timeResolution,t,graph,horizon,controlMethod);     
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
                textHandle3 = text(0.05,0.88,['Ctrl. Cost: ',num2str(controlCost,5)],'Color','k','FontSize',10);
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




% plotting data: if plotMode
if plotMode
    timeSeries = 0:timeResolution:periodT;
    sumOmega = zeros(periodT/timeResolution+1,1);
    for i = 1:1:numOfTargets
        figure()
        subplot(1,3,1)
        plot(timeSeries, data(:,i,1),'r','DisplayName','$\phi_i(t)$')
        hold on 
        plot(timeSeries, data(:,i,2),'b','DisplayName','$\hat{\phi}_i(t)$')
        ylabel('$\phi_i(t)$ and $\hat{\phi}_i(t)$','Interpreter','Latex')
        if targetControllersEnabled
            plot(timeSeries, data(:,i,4),'k','DisplayName','$r_i(t)$')
            ylabel('$\phi_i(t)$, $\hat{\phi}_i(t)$ and $r_i(t)$','Interpreter','Latex')
        end
        legend('Interpreter','Latex','Location','S')
        xlabel('Time - $t$','Interpreter','Latex')
        grid on

        subplot(1,3,2)
        plot(timeSeries, data(:,i,1)-data(:,i,2),'b','DisplayName','$(\phi_i(t)-\hat{\phi}_i(t))$')
        if targetControllersEnabled
            hold on 
            plot(timeSeries, data(:,i,2)-data(:,i,4),'r','DisplayName','$(\hat{\phi}_i(t)-r_i(t))$')
            plot(timeSeries, data(:,i,1)-data(:,i,4),'k','DisplayName','$(\phi_i(t)-r_i(t))$')
        end
        title(['Target ',num2str(i)])
        ylabel('Error Metrics')
        legend('Interpreter','Latex','Location','S')
        xlabel('Time - $t$','Interpreter','Latex')
        grid on

        subplot(1,3,3)
        plot(timeSeries, data(:,i,3),'b')
        ylabel('$\Omega_i(t) = E(e_i e_i^T)$','Interpreter','Latex')
        xlabel('Time - $t$','Interpreter','Latex')
        grid on
        
        sumOmega = sumOmega + data(:,i,3);
    end
    
    sumOmega_old = 0;
    sumVal = 0;
    intSumOmega = zeros(periodT/timeResolution+1,1);
    count = 1;
    factor = 0.5*timeResolution;
    timePeriod = 0;
    for sumOmega_k = sumOmega'
        timePeriod = timePeriod + timeResolution;
        sumVal = sumVal + factor*(sumOmega_k + sumOmega_old);
        intSumOmega(count) = sumVal/timePeriod;
        count = count + 1;
        sumOmega_old = sumOmega_k; 
    end
     figure()
     plot(timeSeries, intSumOmega,'b')
     ylabel('Cost: $J_T$','Interpreter','Latex')
     xlabel('Time: $T$','Interpreter','Latex')
     grid on
     
     
     executionTImes = [];
     for a = 1:1:numOfAgents
        executionTImes = [executionTImes; graph.agents(a).executionTimes(:,3)];
     end
     disp(['Mean execution time: ',num2str(mean(executionTImes))]);
end

if videoMode
    % create the video writer with 1 fps
    writerObj = VideoWriter('myVideo.avi');
    writerObj.FrameRate = 10;
    % set the seconds per image
    % open the video writer
    open(writerObj);
    % write the frames to the video
    for i=1:length(frameArray)
        % convert the image to a frame
        frame = frameArray(i) ;    
        writeVideo(writerObj, frame);
    end
    % close the writer object
    close(writerObj);
end
    
   

% implement periodic controls: !  dwell times: ?
% veryfy the steady state Omega_i


% figure(2)
% p = [0.4,0,0.2,5]
% r = rectangle('Position',p)
% axis([0,1,0,100])
% grid on
% for t = 0:timeResolution:100
% delete(r)
% p = [0.4,0,0.2,5+t/2]
% r = rectangle('Position',p)
% drawnow
% % pause(0.05)
% end
