clear all;
clc;

% % Read value from matlab file
% % Causes on order: Workmen,Materials,Machines,Roller
% [~, causes] = xlsread('Exercise1.xlsx','B2:B5');
% morning = xlsread('Exercise1.xlsx', 'C2:C5')';
% evening = xlsread('Exercise1.xlsx', 'C6:C9')';
% allDay = morning + evening;
% 
% % Histogram
% figure(1);
% subplot(3,1,1);
% histogram('Categories',causes,'BinCounts',morning);
% title('Histogram of Morning shift');
% 
% subplot(3,1,2);
% histogram('Categories',causes,'BinCounts',evening);
% title('Histogram of Evening shift');
% 
% subplot(3,1,3);
% histogram('Categories',causes,'BinCounts',allDay);
% title('Histogram of All day');
% 
% % Create pareto graph -> trova i problemi e gli assegna un peso
% figure(2);
% subplot(3,1,1);
% pareto(morning, causes);
% title('Pareto graph of Morning shift');
% 
% subplot(3,1,2);
% pareto(evening, causes);
% title('Pareto graph of Evening shift');
% 
% subplot(3,1,3);
% pareto(allDay, causes);
% title('Pareto graph of All day shift');

%% 
% "C" Control Chart controlchart

% PER ORA SOLO TEST -> non è la versione finale

% NB: ? is the average number of defective items -> estimate it with Maximum likelihood estimator
lambda = (312 + 245) / 2;  
size = [100 1];
y = poissrnd(lambda, size);
%plot(y);

st = controlchart(y,'charttype',{'c'});

