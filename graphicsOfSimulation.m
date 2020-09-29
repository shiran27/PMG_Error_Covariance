

% load('results/Case1/ST/Case1_RHC.mat')

%% plotting data: if plotMode
if plotMode
    
    timeSeries = 0:timeResolution:periodT;
    sumOmega = zeros(periodT/timeResolution+1,1);
    sumOmegaActive = zeros(periodT/timeResolution+1,1); % to get ratio of sums
    sumOmegaInactive = zeros(periodT/timeResolution+1,1); % to get ratio of sums
    sumIndividualRatio = zeros(periodT/timeResolution+1,1); % to get sum of ratios
    sumIntIndividualRatio = zeros(periodT/timeResolution+1,1); % to get sum of ratios
    maxCovarianceRecorded = zeros(2,numOfTargets); % maxOmega, time 
    for i = 1:1:numOfTargets
        figure()
        subplot(3,1,1)
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

        subplot(3,1,2)
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

        subplot(3,1,3)
        activePart = data(:,i,5)==1;
        inactivePart = ~activePart;
        plot(timeSeries(activePart), data(activePart,i,3),'.b','DisplayName','Active')
        hold on
        plot(timeSeries(inactivePart), data(inactivePart,i,3),'.r','DisplayName','Inactive')
        plot(timeSeries, data(:,i,3),'k','DisplayName','$\Omega_i(t)$')
        ylabel('$\Omega_i(t) = E(e_i e_i^T)$','Interpreter','Latex')
        xlabel('Time - $t$','Interpreter','Latex')
        legend('Location','NW','Interpreter','Latex')
        grid on
        
        
        Omega_i = data(:,i,3);
        eta_i = data(:,i,5);
        sumOmega = sumOmega + Omega_i;
        sumOmegaActive = sumOmegaActive + Omega_i.*eta_i;
        sumOmegaInactive = sumOmegaInactive + Omega_i.*(1-eta_i);
        
        % integrating over each target 
        Omega_iActive = Omega_i.*eta_i;
        
        Omega_i_old = 0;
        Omega_iActive_old = 0;
        sumVal = 0;
        sumValActive = 0;
        intOmega_i = zeros(periodT/timeResolution+1,1);
        intOmega_iActive = zeros(periodT/timeResolution+1,1);
        count = 1;
        factor = 0.5*timeResolution;
        timePeriod = 0;
        for Omega_i_k = Omega_i'
            Omega_iActive_k = Omega_iActive(count);
            timePeriod = timePeriod + timeResolution; % no need for this here
            sumVal = sumVal + factor*(Omega_i_k + Omega_i_old);
            sumValActive = sumValActive + factor*(Omega_iActive_k + Omega_iActive_old);
            intOmega_i(count) = sumVal;
            intOmega_iActive(count) = sumValActive;
            count = count + 1;
            Omega_i_old = Omega_i_k; 
            Omega_iActive_old = Omega_iActive_k; 
        end
        
        sumIndividualRatio = sumIndividualRatio + Omega_iActive./Omega_i;
        sumIntIndividualRatio = sumIntIndividualRatio + intOmega_iActive./intOmega_i;
        
        
        % max omega tracking
        [val,ind] = max(Omega_i);
        maxCovarianceRecorded(1,i) = val; % max omega for target i
        maxCovarianceRecorded(2,i) = ind*timeResolution; % time which occured the max omega
        
    end
    
    [maxOmega, maxTarget] = max(maxCovarianceRecorded(1,:));
    maxTime = maxCovarianceRecorded(2,maxTarget);
    disp(['Max Omega: ',num2str(maxOmega),'; at Target: ',num2str(maxTarget),'; at Time: ',num2str(maxTime)]);
    
    
    % Integrated costs: Objective function vs RatioOfSums vs SumOfRatios 
    sumOmega_old = 0;
    sumOmegaActive_old = 0;
    sumVal = 0;
    sumValActive = 0;
    intSumOmegaN = zeros(periodT/timeResolution+1,1); % this is divided by t
    intSumOmega = zeros(periodT/timeResolution+1,1);
    intSumOmegaActive = zeros(periodT/timeResolution+1,1);
    count = 1;
    factor = 0.5*timeResolution;
    timePeriod = 0;
    for sumOmega_k = sumOmega'
        sumOmegaActive_k = sumOmegaActive(count);
        timePeriod = timePeriod + timeResolution;
        sumVal = sumVal + factor*(sumOmega_k + sumOmega_old);
        sumValActive = sumValActive + factor*(sumOmegaActive_k + sumOmegaActive_old);
        intSumOmegaN(count) = sumVal/timePeriod;
        intSumOmega(count) = sumVal;
        intSumOmegaActive(count) = sumValActive;
        count = count + 1;
        sumOmega_old = sumOmega_k; 
        sumOmegaActive_old = sumOmegaActive_k; 
    end
     figure()
     plot(timeSeries, intSumOmegaN,'DisplayName','Global objective: $J^A+J^I$')
     ylabel('Cov. cost in $[0,t]$: $J_t$','Interpreter','Latex')
     hold on
     yyaxis right
     plot(timeSeries, -intSumOmegaActive./intSumOmega,'DisplayName','RatOfSum: $-J^A/(J^A+J^I)$')
%      plot(timeSeries, -sumIntIndividualRatio,'DisplayName','SumOfRat: $-\sum_i [J_i^A/(J_i^A+J_i^I)]$')
     ylabel('Ratio cost in $[0,t]$','Interpreter','Latex')
     xlabel('Time: $t$','Interpreter','Latex')
     legend('Location','NW','Interpreter','Latex')
     grid on
     
     
     
     
     % Instantaneous costs: components without normalizing with time
     figure()
     plot(timeSeries, sumOmegaActive,'b','DisplayName','Active: $J^A$')
     hold on
     plot(timeSeries, sumOmegaInactive,'r','DisplayName','Inactive: $J^I$')
     plot(timeSeries, sumOmega,'k','DisplayName','Total: $J^A+J^I$')
     yyaxis right 
     plot(timeSeries, -sumOmegaActive./sumOmegaInactive,'DisplayName','Ratio 1: $-J^A/J^I$')
     ylabel('Ratio cost between $[t,t+\Delta]$','Interpreter','Latex')
     yyaxis left
     ylabel('Cost between $[t,t+\Delta]$','Interpreter','Latex')
     xlabel('Time: $t$','Interpreter','Latex')
     legend('Location','NW','Interpreter','Latex')
     grid on
     
     
     
     % Instantaneous costs : Objective function vs RatioOfSums vs SumOfRatios 
     figure()
     plot(timeSeries, sumOmega,'DisplayName','Global objective: $J^A+J^I$')
     ylabel('Cov. cost between $[t,t+\Delta]$','Interpreter','Latex')
     yyaxis right;
     hold on
     plot(timeSeries, -sumOmegaActive./sumOmega,'DisplayName','RatOfSum: $-J^A/(J^A+J^I)$')
%      plot(timeSeries, -sumIndividualRatio,'DisplayName','SumOfRat: $-\sum_i [J_i^A/(J_i^A+J_i^I)]$')
     ylabel('Ratio cost between $[t,t+\Delta]$','Interpreter','Latex')
     xlabel('Time: $t$','Interpreter','Latex')
     legend('Location','NW','Interpreter','Latex')
     grid on
     
     
     % agent distance
     sumDistance = zeros(periodT/timeResolution+1,1);
     for i = 1:1:numOfAgents
        sumDistance = sumDistance + dataAgent(:,i,1);
     end
     
     sumDistance_old = 0;
     sumVal = 0;
     intSumDistance = zeros(periodT/timeResolution+1,1);
     count = 1;
     factor = 0.5;%*timeResolution;
     timePeriod = 0;
     for sumDistance_k = sumDistance'
        timePeriod = timePeriod + timeResolution;
        sumVal = sumVal + factor*(sumDistance_k + sumDistance_old);
        intSumDistance(count) = sumVal; %sumVal/timePeriod;
        count = count + 1;
        sumDistance_old = sumDistance_k; 
     end
     figure()
     plot(timeSeries, intSumDistance,'k','DisplayName','Total Distance')
     ylabel('Distance','Interpreter','Latex')
     xlabel('Time: $T$','Interpreter','Latex')
     legend('Location','NW','Interpreter','Latex')
     grid on
     disp(['Total Distance = ',num2str(intSumDistance(end))]);
     
     
     % executionTimes
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