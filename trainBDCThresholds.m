%% BDC Trainer 
load('results/Case2/ST/Case2_RHC.mat','data','timeResolution','periodT')

numOfTargets = size(data,2);
timeSeries = 0:timeResolution:periodT;
for i = 1:1:numOfTargets
    
    Omega_i = data(:,i,3);
    eta_i = data(:,i,5);
    
    i
    Omega_iU = Omega_i(diff(eta_i)==1)
    Omega_iL = Omega_i(diff(eta_i)==-1)
    
    
    
    activePart = eta_i==1;
    inactivePart = ~activePart;
    
    figure
    plot(timeSeries(activePart), Omega_i(activePart),'.b','DisplayName','Active')
    hold on
    plot(timeSeries(inactivePart), Omega_i(inactivePart),'.r','DisplayName','Inactive')
    plot(timeSeries, Omega_i,'k','DisplayName','$\Omega_i(t)$')
    ylabel('$\Omega_i(t) = E(e_i e_i^T)$','Interpreter','Latex')
    xlabel('Time - $t$','Interpreter','Latex')
    legend('Location','NW','Interpreter','Latex')
    grid on
    
end
