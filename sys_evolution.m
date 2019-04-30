%%%
% Definizione dei parametri del sistema
%%%

clear;

% Raccolta dei dati dal file Excel
theta_a = xlsread('dati.xlsx', 'F1:F25');
theta_c_star = xlsread('dati.xlsx', 'C1:C25');

% Parametri temporali
T = 24; % Orizzonte
sampling_time = .1; % Un secondo
t = 0 : sampling_time : T; % Asse dei tempi
hours = (0 : length(theta_a)-1)'; % Indice delle ore per le temperature

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

% Matrici in tempo continuo
Ac = [-kfc/Cc, kfc/cc; ...
    kfc/Cf, -(kfc+kaf)/Cf];
Bc = [0 1/Cf]';
C = [0 1];
sys_c = ss(Ac,Bc,C,0);

% Matrici in tempo discreto
sys_d = c2d(sys_c, sampling_time);
A = sys_d.A;
B = sys_d.B;

% Distrubi stocastici
% Deviazioni standard (radici della Varianza)
c_std = 1 / Cc;
f_std = 1 / Cf;

% Definizione delle quantità come nel paper
% Le matrici A e B sono già inserite
z = zeros(2, length(t));
u = zeros(2, length(t)-1);
mi = zeros(2, length(t)-1);
eta = [zeros(1, length(t)); interp1(hours, theta_a, t)];
z_hat = [interp1(hours, theta_c_star, t); zeros(1, length(t))];

plot(hours,theta_c_star,'o', t,z_hat(1,:),':.');
