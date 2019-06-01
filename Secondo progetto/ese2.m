clear all;
clc;

% Leggo i tempi di processamento dei job su ogni macchina dal file excel
M1 = xlsread('ese2_data.xlsx', 'B2:K2')';
M2 = xlsread('ese2_data.xlsx', 'B3:K3')';
M3 = xlsread('ese2_data.xlsx', 'B4:K4')';
M4 = xlsread('ese2_data.xlsx', 'B5:K5')';

% In ogni colonna ho i tempi di processamento su quella macchina
M = [M1 M2 M3 M4];

% Inizializzo i job
jobs = 1:length(M1);

% Definisco le macchina virtuali MA e MB per applicare Jonhson
% Per farlo inserisco in MA ed MB gli indici della colonna della rispettiva
% macchina
MA_index = [1 2 3];
MB_index = [4];

% Costruisco le macchine virtuali MA e MB a partire dagli indici
MA = sum(M(:,MA_index), 2);
MB = sum(M(:,MB_index), 2);

% Applico l'algoritmo di Johnson sulle macchine MA e MB
scheduled = JohnsonAlgorithm(jobs, MA', MB');

% Calcola il tempo di computazione totale con questo scheduling Cmax
% Tempi in cui finisce il primo job su MA ed MB
TA = MA(scheduled(1));
TB = MA(scheduled(1)) + MB(scheduled(1));

for i = 2:length(jobs)
    
    TA = TA + MA(scheduled(i));
    
    if (TB < TA)
        % Caso sfortunato in cui la macchina MB finisce prima della macchina 
        % MA e quindi ha un tempo idle
        T_idle = TA - TB;
        TB = TB + T_idle;
    end
    
    TB = TB + MB(scheduled(i));
    
end

Cmax = TB;

