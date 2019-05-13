%%%
% Definizione dei parametri del sistema
%%%

clear;
clc;

% Raccolta dei dati dal file Excel
theta_a_init = xlsread('dati.xlsx', 'F1:F25');
theta_c_star = xlsread('dati.xlsx', 'C1:C25');

% Parametri temporali
T = 24; % Orizzonte
deltaT = 1 / 3600; % Un secondo
t = 0 : deltaT : T; % Asse dei tempi
hours = (0 : length(theta_a_init)-1)'; % Indice delle ore per le temperature

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
A = [-kfc/Cc, kfc/Cc; ...
    kfc/Cf, -(kfc+kaf)/Cf];
B = [0 1/Cf]';
Ba = [0 kaf/Cf]';
C = [0 1];

% Distrubi stocastici
c_std = deltaT / Cc; % Deviazioni standard (radici della Varianza)
f_std = deltaT / Cf;
csi_std = .5^.5;
omega_c = c_std .* randn(1, length(t)-1);
omega_f = f_std .* randn(1, length(t)-1);
omega = [omega_c; omega_f]; % Rumore sullo stato
csi_f = csi_std .* randn(1, length(t)-1); % Rumore sull'uscita
stats = [mean(omega_f) std(omega_f) var(omega_f)];

% Segnale da inseguire
z_hat = [interp1(hours, theta_c_star, t); zeros(1, length(t))];

% Discretizzazione delle matrici
A = A * deltaT + eye(size(A));
B = deltaT * B;
Ba = deltaT * Ba;

% Calcolo Kalman
Q_var = cov(omega'); % Covarianza del rumore sullo stato
R_var = csi_std^2; % Covarianza del rumore sull'uscita
SIGMA_var = eye(2); % Covarianza della stima iniziale dello stato
[sigma, K_kalm] = kalmanFilter(A, C, Q_var, R_var, SIGMA_var, t);

% Definisco le matrici di costo
M = [1 0; 0 1e-6]; % Si vuole inseguire solo una delle due componenti
N = 1e-3; % Controllo scalare
MT = M;

% Disturbo deterministico (scalato)
theta_a = interp1(hours, theta_a_init, t);
mi = zeros(2, length(t));
for i = 1 : length(t)
   mi(:,i) = theta_a(i) * Ba;
end

% Calcolo di z_d2
z_d2 = zeros(2, length(t));
z_d2(:,1) = [293; 291];
for i = 1 : length(t) - 1
   z_d2(:,i+1) = A * z_d2(:,i) + mi(:,i);
end

% Effettuo il controllo "standard"
[u, z] = standard_control(A, B, C, omega, csi_f, M, MT, N, K_kalm, z_hat, z_d2, t);

% subplot(1,2,1);
% stairs(t, u');
% legend('u');
% 
% subplot(1,2,2);
% stairs(t, [z' z_hat(1,:)' theta_a']);
% legend('theta_c', 'theta_f', 'theta_c*', 'theta_a');

% Performance di controllo
% Errore qudratico medio
mse = sum((z(1,:) - z_hat(1,:)).^2) / length(t);

% Deviazione massima di temperatura
max_abs_dev = max(abs((z(1,:) - z_hat(1,:))));

% Consumo di energia nelle 24 ore
energy = sum(abs(u));
% Consumo celle frigo da internet: da 1.41e7 a 2.47e7 J/24h

% Varianza del controllo
u_var = var(u);

% Risolvo il problema con rolling horizon optimization con time windows
% differenti. 

[u_mpc, z_mpc] = model_predictive_control(A, B, C, omega, csi_f, M, MT, N, K_kalm, z_hat, z_d2, t, 9 * 60);

% Performance di controllo
% Errore qudratico medio
mse_mpc = sum((z_mpc(1,:) - z_hat(1,:)).^2) / length(t);

% Deviazione massima di temperatura
max_abs_dev_mpc = max(abs((z_mpc(1,:) - z_hat(1,:))));

% Consumo di energia nelle 24 ore
energy_mpc = sum(abs(u_mpc));
% Consumo celle frigo da internet: da 1.41e7 a 2.47e7 J/24h

% Varianza del controllo
u_var_mpc = var(u_mpc);

% Differenze tra le prestazioni
diff_mse = mse - mse_mpc;
diff_max_abs_dev = max_abs_dev - max_abs_dev_mpc;
diff_energy = energy - energy_mpc;
diff_u_var = u_var - u_var_mpc;

% Matlab toolbox -> https://it.mathworks.com/help/mpc/gs/control-of-a-multi-input-single-output-plant.html


% figure
% subplot(2,1,1)
% plot(0:Tf-1,y,0:Tf-1,r)
% title('Output')
% grid
% subplot(2,1,2)
% plot(0:Tf-1,u)
% title('Input')
% grid


% stairs(t, [z_cmp' z_hat(1,:)' theta_a']);
% legend('theta_c cmp', 'theta_f cmp', 'theta_c*', 'theta_a');


stairs(t, [z(1,:)' z_mpc(1,:)' z_hat(1,:)' theta_a']);
legend('theta_c standard', 'theta_c cmp', 'theta_c STAR', 'theta_a');
