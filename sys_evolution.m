%%%
% Definizione dei parametri del sistema
%%%

clear;

% Costanti numeriche
cc = 2.15;
mc = 2;
Cc = cc * mc;

cf = 1;
mf = 1.3;
Cf = cf * mf;

sc = .095; % squared meters
k_bar_fc = 100;
kfc = k_bar_fc * sc;

sf = 6;
k_bar_af = 1.22;
kaf = k_bar_af * sf;

% Disturbo deterministico
theta_a = xlsread('dati.xlsx', 'F1:F25');
hours = (0 : length(theta_a)-1)';

% Distrubi stocastici
% Deviazioni standard (radici della Varianza)
c_std = 1 / Cc;
f_std = 1 / Cf;

% Parametri temporali
T = 24; % Orizzonte
sampling_time = 1 / 3600; % Un secondo