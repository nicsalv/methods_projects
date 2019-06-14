clear;
clc;

% Read value from matlab file
% Causes on order: Workmen, Materials, Machines, Roller
[~, causes] = xlsread('Exercise1.xlsx','B2:B5');
morning = xlsread('Exercise1.xlsx', 'C2:C5')';
evening = xlsread('Exercise1.xlsx', 'C6:C9')';
allDay = morning + evening;

% Histogram
figure(1);
subplot(3,1,1);
histogram('Categories',causes,'BinCounts',morning);
title('Histogram of Morning shift');

subplot(3,1,2);
histogram('Categories',causes,'BinCounts',evening);
title('Histogram of Evening shift');

subplot(3,1,3);
histogram('Categories',causes,'BinCounts',allDay);
title('Histogram of All day');

% Create pareto graph -> trova i problemi e gli assegna un peso
figure(2);
subplot(3,1,1);
pareto(morning, causes);
title('Pareto graph of Morning shift');

subplot(3,1,2);
pareto(evening, causes);
title('Pareto graph of Evening shift');

subplot(3,1,3);
pareto(allDay, causes);
title('Pareto graph of All day shift');

% "C" Control Chart
% In such a case the incidence of defects might be modeled
% as a Poisson distribution

indexOfSum = [1 2 3 4];

morningTotal = sum(morning(indexOfSum));
eveningTotal = sum(evening(indexOfSum));
allDayTotal = sum(allDay(indexOfSum)) / 2;

lambda = allDayTotal;
size = [50 1];
y = poissrnd(lambda, size);
%plot(y);

% Manipolazione della distribuzione per generare
% situazioni fuori controllo.
minim = min(y);
maxim = max(y);
% std = std(y);

% Limiti di disturbo
m = -80;
M = 80;

numOfCases = 10;
for i = 1 : numOfCases
    % Valori casuali potenzialmente critici da inserire nel dataset
    criticalValue = randi([m M]);
    
    % Inserimento dei valori potenzialmente critici in posizioni casuali
    % del dataset
    criticalIndex = randi([1 length(y)]);
    y(criticalIndex) = y(criticalIndex) + criticalValue;
end

figure(3);
controlchart(y,'charttype',{'c'});
title('Control chart - all Day & all Causes');
