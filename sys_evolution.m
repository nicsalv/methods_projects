%%%
% Definizione dei parametri del sistema
%%%

clear;

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
omega_tilde = [omega_c; omega_f] .* deltaT;
csi_f = csi_std .* randn(1, length(t)-1);
stats = [mean(omega_f) std(omega_f) var(omega_f)];

% Risposta libera del sistema discreto
x = zeros(2, length(t));
q = ones(1, length(t)-1) * 0;
y = zeros(1, length(t));
x(:,1) = [293; 291];

theta_a = interp1(hours, theta_a, t);
eta = zeros(2, length(t)-1);

A_tilde = A .* deltaT + eye(size(A));
B_tilde = deltaT .* B;
Ba_tilde = deltaT .* Ba;

for k = 1 : length(t) - 1
    % Salvataggio della sequenza nota (utile per dopo)
    eta(:,k) = Ba_tilde * theta_a(k);
    % La k-esima colonna di x corrisponde all'istante k*deltaT
    x(:,k+1) = A_tilde * x(:,k) + B_tilde * q(:,k) ...
        + eta(:,k) + omega_tilde(:,k);
    y(k) = C * x(:,k) + csi_f(k);
end
y(end) = y(end-1); % Solamente per plottare correttamente
% stairs(t, [x' theta_a']);
% legend('x1','x2','theta_a');



% Controllo LQG/LQT come indicato nel paper

z = zeros(2, length(t));
z_hat = [278; 0]; % Il tracking, per ora, è solamente un setpoint
z0 = [10; 0];

% Calcolo del controllo
M = [1 0; 0 0]; % Stato di due dimensioni
N = 1; % Controllo scalare
MT = M;
[K,Kg,mu,z_d2,g] = Riccati_nonStandard_LQG...
    (M, N, MT, eta, z_hat, A_tilde, B_tilde, z0, t);

% Evoluzione del sistema controllato
for k = 1 : length(t)-1
    q(k) = K(:,:,k) * (z(:,k) - z_d2(:,k)) + Kg(:,:,k) * g(:,k+1);

    z(:,k+1) = A_tilde * z(:,k) + B_tilde * q(:,k) ...
        + mu(:,k) + omega_tilde(:,k);
    
    y(k) = C * x(:,k) + csi_f(k);

end

stairs(t, [(z(1,:)+z_hat(1))' x(1,:)']);