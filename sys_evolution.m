%%%
% Definizione dei parametri del sistema
%%%

clear;
clc;

% Raccolta dei dati dal file Excel
theta_a = xlsread('dati.xlsx', 'F1:F25');
theta_c_star = xlsread('dati.xlsx', 'C1:C25');

% Parametri temporali
T = 24; % Orizzonte
deltaT = 1 / 3600; % Un secondo
t = 0 : deltaT : T; % Asse dei tempi
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

% Variabili di interesse
q = zeros(1, length(t)-1);
y = zeros(1, length(t));

% Segnale da inseguire
z_hat = [interp1(hours, theta_c_star, t); zeros(1, length(t))];

% Discretizzazione delle matrici
A = A * deltaT + eye(size(A));
B = deltaT * B;
Ba = deltaT * Ba;

% Controllo LQG/LQT come indicato nel paper: parte stocastica del sistema.
% Si risolve mediante un LQG standard,
% trascurando il disturbo deterministico.

Q_var = cov(omega'); % Covarianza del rumore sullo stato
R_var = csi_std^2; % Covarianza del rumore sull'uscita
SIGMA_var = eye(2); % Covarianza della stima iniziale dello stato
[sigma, K_kalm] = kalmanFilter(A, C, Q_var, R_var, SIGMA_var, t);

% Parte stocastica del controllo
u_s = zeros(1, length(t));

% Parte stocastica dello stato
z_s = zeros(2, length(t));
z_estimated = zeros(2, length(t));

% Determinazione della matrice K per il controllo ottimo
M = [1 0; 0 1e-6]; % Si vuole inseguire solo una delle due componenti
N = 1e-5; % Controllo scalare
MT = M;
[P, K, K_infinito] = riccati_P_K(A, B, M, MT, N, t);

% Stima e conti all'istante iniziale
y(1) = C * z_s(:,1) + csi_f(1);
z_estimated(:,1) = z_s(:,1) + K_kalm(:,:,1) * (y(1) - C*z_s(:,1));
u_s(1) = K(:,:,1) * z_estimated(:,1);
z_s(:,2) = A * z_s(:,1) + B * u_s(1) + omega(:,1);

for i = 2 : length(t) - 1
    % Misurazione dell'uscita del sistema
    y(i) = C * z_s(:,i) + csi_f(i);
    
    % Stima dello stato mediante Kalman
    sys_kalm = A * z_estimated(:,i-1) + B * u_s(:,i-1);
    z_estimated(:,i) = sys_kalm + K_kalm(:,:,i) * (y(i) - C * sys_kalm); 
    
    % Calcolo del controllo
    u_s(i) = K(:,:,i) * z_estimated(:,i);
    
    % Evoluzione del sistema affetto dal controllo ottimo
    z_s(:,i+1) = A * z_s(:,i) + B * u_s(i) + omega(:,i);
end

% subplot(2,1,1);
% stairs(t, [z_s(1,:)', z_estimated(1,:)']);
% legend('z_s', 'z_estimated');
% 
% subplot(2,1,2);
% stairs(t, y');
% legend('y');


% Parte deterministica del sistema.
% Si tiene conto del disturbo deterministico theta_a.

% Disturbo deterministico (scalato)
theta_a = interp1(hours, theta_a, t);
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

% Segnale da tracciare. Il riferimento (theta_c_star) è in z_hat.
r = z_hat - z_d2;

[K, Kg, g] = Riccati_nonStandard_LQG(M, N, MT, A, B, r, t);

u_d = zeros(1, length(t));
z_d1 = zeros(2, length(t));

for i = 1 : length(t)-1
    % Calcolo del controllo ottimo
    u_d(:,i) = K(:,:,i) * z_d1(:,i) + Kg(:,:,i) * g(:,i+1);
    
    % Evoluzione del sistema affetto dal controllo
    z_d1(:,i+1) = A * z_d1(:, i) + B * u_d(i);
end

% Ricostruzione del controllo e dello completo
u = u_d + u_s;
z = z_s + z_d1 + z_d2;

subplot(1,2,1);
stairs(t, u');
legend('u');

subplot(1,2,2);
stairs(t, [z' z_hat(1,:)' theta_a']);
legend('theta_c', 'theta_f', 'thetaC*', 'theta_a');