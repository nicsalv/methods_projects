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
omega_tilde = [omega_c; omega_f];
csi_f = csi_std .* randn(1, length(t)-1);
stats = [mean(omega_f) std(omega_f) var(omega_f)];

% Risposta libera del sistema discreto
x = zeros(2, length(t));
q = ones(1, length(t)-1) * 0;
y = zeros(1, length(t));
x(:,1) = [293; 291];

theta_a = interp1(hours, theta_a, t);

A_tilde = A .* deltaT + eye(size(A));
B_tilde = deltaT .* B;
Ba_tilde = deltaT .* Ba;

for k = 1 : length(t) - 1
    % La k-esima colonna di x corrisponde all'istante k*deltaT
    x(:,k+1) = A_tilde * x(:,k) + B_tilde * q(:,k) ...
        + Ba_tilde * theta_a(k) + omega_tilde(:,k);
    y(k) = C * x(:,k) + csi_f(k);
end
y(end) = y(end-1); % Solamente per plottare correttamente

stairs(t, [x' theta_a' y']);
legend('x1','x2','theta_a','y');

eig(A_tilde)
 
% % Definizione delle quantità come nel paper
% % Le matrici A e B sono già inserite
% z = zeros(2, length(t));
% u = zeros(2, length(t)-1);
% mi = zeros(2, length(t)-1);
% eta = [zeros(1, length(t)); interp1(hours, theta_a, t)];
% z_hat = [interp1(hours, theta_c_star, t); zeros(1, length(t))];
% 
% plot(hours,theta_c_star,'o', t,z_hat(1,:),':.');
