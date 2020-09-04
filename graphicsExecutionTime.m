clear all 
close all
clc

%% Laptop vs PC: 6 vs 5

%%% Laptop:

load('Results\ET\Case1L_6.mat')

t = [];
et = [];
for agentId = 1:1:2
    t = [t; agents(agentId).executionTimes(:,1)];
    et = [et; agents(agentId).executionTimes(:,3)];
end

figure(1)
plot(t,et,'.','DisplayName','Laptop-Learning')


windowSize = 10; % seconds
stepSize = 0.01;
smoothArray = [];
% countArray = [];
for tVal = 0:stepSize:(t(end)-windowSize)
    window = et(t>tVal & t<tVal+windowSize);
    smoothArray = [smoothArray; tVal,mean(window)];
%     countArray = [countArray; length(window)];
end

figure(2)
plot(smoothArray(:,1),smoothArray(:,2),'DisplayName','Laptop-Learning')



load('Results\ET\Case1_6.mat')

t = [];
et = [];
for agentId = 1:1:2
    t = [t; agents(agentId).executionTimes(:,1)];
    et = [et; agents(agentId).executionTimes(:,3)];
end

figure(1)
hold on
plot(t,et,'.','DisplayName','Laptop')
grid on
legend('Location','NW')
xlabel('Time')
ylabel('Execution time for a RHCP')


smoothArray = [];
for tVal = 0:stepSize:(t(end)-windowSize)
    window = et(t>tVal & t<tVal+windowSize);
    smoothArray = [smoothArray; tVal,mean(window)];
end

figure(2)
hold on
plot(smoothArray(:,1),smoothArray(:,2),'DisplayName','Laptop')
legend('Location','NW')
grid on
xlabel('Time')
ylabel('Average execution time for a RHCP')




%%% PC 

load('Results\ET\Case1L_5.mat')

t = [];
et = [];
for agentId = 1:1:2
    t = [t; agents(agentId).executionTimes(:,1)];
    et = [et; agents(agentId).executionTimes(:,3)];
end

figure(1)
plot(t,et,'.','DisplayName','PC-Learning')


windowSize = 10; % seconds
stepSize = 0.01;
smoothArray = [];
% countArray = [];
for tVal = 0:stepSize:(t(end)-windowSize)
    window = et(t>tVal & t<tVal+windowSize);
    smoothArray = [smoothArray; tVal,mean(window)];
%     countArray = [countArray; length(window)];
end

figure(2)
plot(smoothArray(:,1),smoothArray(:,2),'DisplayName','PC-Learning')



load('Results\ET\Case1_5.mat')

t = [];
et = [];
for agentId = 1:1:2
    t = [t; agents(agentId).executionTimes(:,1)];
    et = [et; agents(agentId).executionTimes(:,3)];
end

figure(1)
hold on
plot(t,et,'.','DisplayName','PC')
grid on
legend('Location','NW')
xlabel('Time')
ylabel('Execution time for a RHCP')


smoothArray = [];
for tVal = 0:stepSize:(t(end)-windowSize)
    window = et(t>tVal & t<tVal+windowSize);
    smoothArray = [smoothArray; tVal,mean(window)];
end

figure(2)
hold on
plot(smoothArray(:,1),smoothArray(:,2),'DisplayName','PC')
legend('Location','NW')
grid on
xlabel('Time')
ylabel('Average execution time for a RHCP')


%% PC vs PC with data accumilating: 5 vs 4

%%% PC

load('Results\ET\Case1L_5.mat')

t = [];
et = [];
for agentId = 1:1:2
    t = [t; agents(agentId).executionTimes(:,1)];
    et = [et; agents(agentId).executionTimes(:,3)];
end

figure(3)
plot(t,et,'.','DisplayName','PC-Learning')


windowSize = 10; % seconds
stepSize = 0.01;
smoothArray = [];
% countArray = [];
for tVal = 0:stepSize:(t(end)-windowSize)
    window = et(t>tVal & t<tVal+windowSize);
    smoothArray = [smoothArray; tVal,mean(window)];
%     countArray = [countArray; length(window)];
end

figure(4)
plot(smoothArray(:,1),smoothArray(:,2),'DisplayName','PC-Learning')



load('Results\ET\Case1_5.mat')

t = [];
et = [];
for agentId = 1:1:2
    t = [t; agents(agentId).executionTimes(:,1)];
    et = [et; agents(agentId).executionTimes(:,3)];
end

figure(3)
hold on
plot(t,et,'.','DisplayName','PC')
grid on
legend('Location','NW')
xlabel('Time')
ylabel('Execution time for a RHCP')


smoothArray = [];
for tVal = 0:stepSize:(t(end)-windowSize)
    window = et(t>tVal & t<tVal+windowSize);
    smoothArray = [smoothArray; tVal,mean(window)];
end

figure(4)
hold on
plot(smoothArray(:,1),smoothArray(:,2),'DisplayName','PC')
legend('Location','NW')
grid on
xlabel('Time')
ylabel('Average execution time for a RHCP')



%%% PC - Data Accumilating

load('Results\ET\Case1L_4.mat')

t = [];
et = [];
for agentId = 1:1:2
    t = [t; agents(agentId).executionTimes(:,1)];
    et = [et; agents(agentId).executionTimes(:,3)];
end

figure(3)
plot(t,et,'.','DisplayName','PC-Mem-Learning')


windowSize = 10; % seconds
stepSize = 0.01;
smoothArray = [];
% countArray = [];
for tVal = 0:stepSize:(t(end)-windowSize)
    window = et(t>tVal & t<tVal+windowSize);
    smoothArray = [smoothArray; tVal,mean(window)];
%     countArray = [countArray; length(window)];
end

figure(4)
plot(smoothArray(:,1),smoothArray(:,2),'DisplayName','PC-Mem-Learning')



load('Results\ET\Case1_4.mat')

t = [];
et = [];
for agentId = 1:1:2
    t = [t; agents(agentId).executionTimes(:,1)];
    et = [et; agents(agentId).executionTimes(:,3)];
end

figure(3)
hold on
plot(t,et,'.','DisplayName','PC-Mem')
grid on
legend('Location','NW')
xlabel('Time')
ylabel('Execution time for a RHCP')


smoothArray = [];
for tVal = 0:stepSize:(t(end)-windowSize)
    window = et(t>tVal & t<tVal+windowSize);
    smoothArray = [smoothArray; tVal,mean(window)];
end

figure(4)
hold on
plot(smoothArray(:,1),smoothArray(:,2),'DisplayName','PC-Mem')
legend('Location','NW')
grid on
xlabel('Time')
ylabel('Average execution time for a RHCP')


%% PC-Mem vs PC-Mem-TicToc : 4 vs 1

load('Results\ET\Case1L_4.mat')

t = [];
et = [];
for agentId = 1:1:2
    t = [t; agents(agentId).executionTimes(:,1)];
    et = [et; agents(agentId).executionTimes(:,3)];
end

figure(5)
plot(t,et,'.','DisplayName','PC-Mem-Learning')


windowSize = 10; % seconds
stepSize = 0.01;
smoothArray = [];
% countArray = [];
for tVal = 0:stepSize:(t(end)-windowSize)
    window = et(t>tVal & t<tVal+windowSize);
    smoothArray = [smoothArray; tVal,mean(window)];
%     countArray = [countArray; length(window)];
end

figure(6)
plot(smoothArray(:,1),smoothArray(:,2),'DisplayName','PC-Mem-Learning')



load('Results\ET\Case1_4.mat')

t = [];
et = [];
for agentId = 1:1:2
    t = [t; agents(agentId).executionTimes(:,1)];
    et = [et; agents(agentId).executionTimes(:,3)];
end

figure(5)
hold on
plot(t,et,'.','DisplayName','PC-Mem')
grid on
legend('Location','NW')
xlabel('Time')
ylabel('Execution time for a RHCP')


smoothArray = [];
for tVal = 0:stepSize:(t(end)-windowSize)
    window = et(t>tVal & t<tVal+windowSize);
    smoothArray = [smoothArray; tVal,mean(window)];
end

figure(6)
hold on
plot(smoothArray(:,1),smoothArray(:,2),'DisplayName','PC-Mem')
legend('Location','NW')
grid on
xlabel('Time')
ylabel('Average execution time for a RHCP')


load('Results\ET\Case1L.mat')

t = [];
et = [];
for agentId = 1:1:2
    t = [t; agents(agentId).executionTimes(:,1)];
    et = [et; agents(agentId).executionTimes(:,3)];
end

figure(5)
plot(t,et,'.','DisplayName','PC-Mem-Tic-Learning')


windowSize = 10; % seconds
stepSize = 0.01;
smoothArray = [];
% countArray = [];
for tVal = 0:stepSize:(t(end)-windowSize)
    window = et(t>tVal & t<tVal+windowSize);
    smoothArray = [smoothArray; tVal,mean(window)];
%     countArray = [countArray; length(window)];
end

figure(6)
plot(smoothArray(:,1),smoothArray(:,2),'DisplayName','PC-Mem-Tic-Learning')



load('Results\ET\Case1.mat')

t = [];
et = [];
for agentId = 1:1:2
    t = [t; agents(agentId).executionTimes(:,1)];
    et = [et; agents(agentId).executionTimes(:,3)];
end

figure(5)
hold on
plot(t,et,'.','DisplayName','PC-Mem-Tic')
grid on
legend('Location','NW')
xlabel('Time')
ylabel('Execution time for a RHCP')


smoothArray = [];
for tVal = 0:stepSize:(t(end)-windowSize)
    window = et(t>tVal & t<tVal+windowSize);
    smoothArray = [smoothArray; tVal,mean(window)];
end

figure(6)
hold on
plot(smoothArray(:,1),smoothArray(:,2),'DisplayName','PC-Mem-Tic')
legend('Location','NW')
grid on
xlabel('Time')
ylabel('Average execution time for a RHCP')


%% PC-Mem-TIcToc vs PC-Mem-TicToc-Matlab solver : 1 vs 2

load('Results\ET\Case1L.mat')

t = [];
et = [];
for agentId = 1:1:2
    t = [t; agents(agentId).executionTimes(:,1)];
    et = [et; agents(agentId).executionTimes(:,3)];
end

figure(7)
plot(t,et,'.','DisplayName','PC-Mem-Tic-Learning')


windowSize = 10; % seconds
stepSize = 0.01;
smoothArray = [];
% countArray = [];
for tVal = 0:stepSize:(t(end)-windowSize)
    window = et(t>tVal & t<tVal+windowSize);
    smoothArray = [smoothArray; tVal,mean(window)];
%     countArray = [countArray; length(window)];
end

figure(8)
plot(smoothArray(:,1),smoothArray(:,2),'DisplayName','PC-Mem-Tic-Learning')



load('Results\ET\Case1.mat')

t = [];
et = [];
for agentId = 1:1:2
    t = [t; agents(agentId).executionTimes(:,1)];
    et = [et; agents(agentId).executionTimes(:,3)];
end

figure(7)
hold on
plot(t,et,'.','DisplayName','PC-Mem-Tic')
grid on
legend('Location','NW')
xlabel('Time')
ylabel('Execution time for a RHCP')


smoothArray = [];
for tVal = 0:stepSize:(t(end)-windowSize)
    window = et(t>tVal & t<tVal+windowSize);
    smoothArray = [smoothArray; tVal,mean(window)];
end

figure(8)
hold on
plot(smoothArray(:,1),smoothArray(:,2),'DisplayName','PC-Mem-Tic')
legend('Location','NW')
grid on
xlabel('Time')
ylabel('Average execution time for a RHCP')


%%% PC-Mem-Tic-Matlab
load('Results\ET\Case1L_2.mat')

t = [];
et = [];
for agentId = 1:1:2
    t = [t; agents(agentId).executionTimes(:,1)];
    et = [et; agents(agentId).executionTimes(:,3)];
end

figure(7)
plot(t,et,'.','DisplayName','PC-Mem-Tic-Matlab-Learning')


windowSize = 10; % seconds
stepSize = 0.01;
smoothArray = [];
% countArray = [];
for tVal = 0:stepSize:(t(end)-windowSize)
    window = et(t>tVal & t<tVal+windowSize);
    smoothArray = [smoothArray; tVal,mean(window)];
%     countArray = [countArray; length(window)];
end

figure(8)
plot(smoothArray(:,1),smoothArray(:,2),'DisplayName','PC-Mem-Tic-Matlab-Learning')



load('Results\ET\Case1_2.mat')

t = [];
et = [];
for agentId = 1:1:2
    t = [t; agents(agentId).executionTimes(:,1)];
    et = [et; agents(agentId).executionTimes(:,3)];
end

figure(7)
hold on
plot(t,et,'.','DisplayName','PC-Mem-Tic-Matlab')
grid on
legend('Location','NW')
xlabel('Time')
ylabel('Execution time for a RHCP')


smoothArray = [];
for tVal = 0:stepSize:(t(end)-windowSize)
    window = et(t>tVal & t<tVal+windowSize);
    smoothArray = [smoothArray; tVal,mean(window)];
end

figure(8)
hold on
plot(smoothArray(:,1),smoothArray(:,2),'DisplayName','PC-Mem-Tic-Matlab')
legend('Location','NW')
grid on
xlabel('Time')
ylabel('Average execution time for a RHCP')




%% Laptop Long run : 7 Vs 3


%%% Laptop

load('Results\ET\Case1_7.mat')

t = [];
et = [];
for agentId = 1:1:2
    t = [t; agents(agentId).executionTimes(:,1)];
    et = [et; agents(agentId).executionTimes(:,3)];
end

figure(9)
plot(t,et,'.','DisplayName','Laptop')


windowSize = 10; % seconds
stepSize = 0.01;
smoothArray = [];
% countArray = [];
for tVal = 0:stepSize:(t(end)-windowSize)
    window = et(t>tVal & t<tVal+windowSize);
    smoothArray = [smoothArray; tVal,mean(window)];
%     countArray = [countArray; length(window)];
end

figure(10)
plot(smoothArray(:,1),smoothArray(:,2),'DisplayName','Laptop')


%%% PC - Mem - Tic

load('Results\ET\Case1_3.mat')

t = [];
et = [];
for agentId = 1:1:2
    t = [t; agents(agentId).executionTimes(:,1)];
    et = [et; agents(agentId).executionTimes(:,3)];
end

figure(9)
hold on
plot(t,et,'.','DisplayName','PC-Mem-Tic')
grid on
legend('Location','NW')
xlabel('Time')
ylabel('Execution time for a RHCP')


smoothArray = [];
for tVal = 0:stepSize:(t(end)-windowSize)
    window = et(t>tVal & t<tVal+windowSize);
    smoothArray = [smoothArray; tVal,mean(window)];
end

figure(10)
hold on
plot(smoothArray(:,1),smoothArray(:,2),'DisplayName','PC-Mem-Tic')
legend('Location','NW')
grid on
xlabel('Time')
ylabel('Average execution time for a RHCP')



%% PC Optimized


load('Results\ET\Case1L_9.mat')

t = [];
et = [];
for agentId = 1:1:2
    count = agents(agentId).executionTimes(end,end);
    t = [t; agents(agentId).executionTimes(1:count,1)];
    et = [et; agents(agentId).executionTimes(1:count,3)];
end

figure(3)
plot(t,et,'.','DisplayName','PC-Optimized-Learn')


windowSize = 10; % seconds
stepSize = 0.01;
smoothArray = [];
% countArray = [];
for tVal = 0:stepSize:(t(end)-windowSize)
    window = et(t>tVal & t<tVal+windowSize);
    smoothArray = [smoothArray; tVal,mean(window)];
%     countArray = [countArray; length(window)];
end

figure(4)
plot(smoothArray(:,1),smoothArray(:,2),'DisplayName','PC-Optimized-Learn')



%%% PC Optimized 
load('Results\ET\Case1_9.mat')

t = [];
et = [];
for agentId = 1:1:2
    count = agents(agentId).executionTimes(end,end);
    t = [t; agents(agentId).executionTimes(1:count,1)];
    et = [et; agents(agentId).executionTimes(1:count,3)];
end

figure(3)
hold on
plot(t,et,'.','DisplayName','PC-Optimized')
grid on
legend('Location','NW')
xlabel('Time')
ylabel('Execution time for a RHCP')


smoothArray = [];
for tVal = 0:stepSize:(t(end)-windowSize)
    window = et(t>tVal & t<tVal+windowSize);
    smoothArray = [smoothArray; tVal,mean(window)];
end

figure(4)
hold on
plot(smoothArray(:,1),smoothArray(:,2),'DisplayName','PC-Optimized')
legend('Location','NW')
grid on
xlabel('Time')
ylabel('Average execution time for a RHCP')


