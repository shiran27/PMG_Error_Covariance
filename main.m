clear all
close all
clc


% create a random graph
graph = Graph(0); % create the graph object : index 0 for main graph 

%Case 1: 7, 7, 0.7, 2;   Case 2: 250, 7, 0.45, 2;   Case 3: 500, 10, 0.45, 4
rand('seed', 7); %% 7, 250, 500
numOfTargets = 7; % 7, 10
dimentionOfSpace = 2;
sizeOfSpace = 1;  % 1x1 square
communicationRadius = 0.7; %0.7,  0.45

% load the agents
numOfAgents = 2; %2, 4

% Control Options: Periodic, RHC-Periodic, RHC-Learn, RHC, variable, fixed, random
controlMethod = 'Periodic';

% create the graph and draw
graph = graph.loadARandomGraph(numOfTargets, numOfAgents, dimentionOfSpace, sizeOfSpace, communicationRadius);
graph.drawGraph(1);

% if we try the periodic controller
if isequal(controlMethod,'RHC-Periodic') | isequal(controlMethod,'Periodic')
    graph.assignCyclesToAgents(controlMethod);
end


timeResolution = 0.001;
periodT = 50; % 50
fixedHorizonValue = 2.5; 
cost = 0;  
textHandle1 = text();  textHandle2 = text(); 
plotMode = false; 
videoMode = true;
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
    
    % agent update
    horizon = min(fixedHorizonValue,periodT-t);
    for i = 1:1:numOfAgents
        graph.agents(i).update(timeResolution,t,graph,horizon,controlMethod);     
    end
    
    % display & video
    if rem(t,0.1)==0
        time = t
        graph.updateGraph(1,t);
        if t~=0
            delete(textHandle1);
            textHandle1 = text(0.05,0.96,['Cost: ',num2str(cost/t,5)],'Color','k','FontSize',10);
        end
        delete(textHandle2);
        textHandle2 = text(0.05,0.92,['Time: ',num2str(t,5)],'Color','k','FontSize',10);
        
        frameArray(frameCount) = getframe(gcf);
        frameCount = frameCount + 1;
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
        plot(timeSeries, data(:,i,1),'r')
        hold on 
        plot(timeSeries, data(:,i,2),'b')
        legend('$\phi_i(t)$','$\hat{\phi}_i(t)$','Interpreter','Latex')
        ylabel('$\phi_i(t)$ and $\hat{\phi}_i(t)$','Interpreter','Latex')
        xlabel('Time - $t$','Interpreter','Latex')
        grid on

        subplot(1,3,2)
        plot(timeSeries, data(:,i,1)-data(:,i,2),'b')
        title(['Target ',num2str(i)])
        ylabel('$e_i(t) = (\phi_i(t)-\hat{\phi}_i(t))$','Interpreter','Latex')
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
