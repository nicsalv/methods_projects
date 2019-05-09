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
omega = [omega_c; omega_f]; % Rumore sullo stato
csi_f = csi_std .* randn(1, length(t)-1); % Rumore sull'uscita
stats = [mean(omega_f) std(omega_f) var(omega_f)];

% Variabili di interesse
z = zeros(2, length(t));
q = zeros(1, length(t)-1);
y = zeros(1, length(t));
z(:,1) = [0; 291]; % Si parte da un punto già ottimo

% Segnale da inseguire
z_hat = [interp1(hours, theta_c_star, t); zeros(1, length(t))];

theta_a = interp1(hours, theta_a, t); % Disturbo deterministico
eta = zeros(2, length(t)-1); % Sequenze note
mi = zeros(2, length(t)-1);

% Discretizzazione delle matrici
A = A .* deltaT + eye(size(A));
B = deltaT .* B;
Ba = deltaT .* Ba;

% Risposta del sistema ad un ingresso nullo
for k = 1 : length(t) - 1
    % Salvataggio delle sequenze note (servono per il controllo)
    eta(:,k) = Ba * theta_a(k);
    
    % La k-esima colonna di x corrisponde all'istante k*deltaT
    z(:,k+1) = A * z(:,k) + mi(:,k) + omega(:,k);
    y(k) = C * z(:,k) + csi_f(k);
end

% Controllo LQG/LQT come indicato nel paper

% Calcolo del controllo
M = [1 0; 0 0.0001]; % Si vuole inseguire solo una delle due componenti
N = 1; % Controllo scalare
MT = M;
%[K,Kg,z_d2,g] = Riccati_nonStandard_LQG(M, N, MT, mi, A, B, z(:,1), t);

z0 = [0 291]';

% Parte stocastica del sistema -> trascuro theta_a -> LQG

% TODO: calcola matrici di covarianza
SIGMA_varianza = eye(2);
[sigma, kalm] = kalmanFilter(A, C, Q_varianza, R_varianza, SIGMA_varianza, t);

us = zeros(1, length(t));

zs = zeros(2, length(t));
z_estimate = zeros(2, length(t));

zs(:,1) = z0;

% TODO: compute k for control -> Riccati_K_G

for i = 1 : length(t) - 1
    
    y(i) = C * zs(:,i) + + csi_f(i);
    
    if ( i == 1 )
        z_estimate(:,i) = zs(:,i) + kalm(:,:,i) * (y(1) - C * zs(:,1));
    else
        z_estimate(:,i) = A * z_estimate(:, i-1) + B * u(:,i-1) + ...
            kalm(:,:,i) * (y(i) - C * ( A * z_estimate(:,i-1) + B * u(:,i-1)));
    end
    
    us(i) = K(i) * z_estimate(:,i);
    
    zs(:,i+1) = A * zs(:,i) + B * us(i) + omega;
    
end

% Parte deterministica del sistema -> considera theta_a -> LQT

% Compute zd2

zd2 = zeros(dim_state, dim);
zd2(:,1) = z0;
for t = 1 : length(t) - 1
    zd2(:,t+1) = A * zd2(:,t) + eta(:,t);
end

% Trovo il riferimento -> il segnale dove sono - dove voglio essere
rt = -zd2 + z_hat;
    % Dove inserisco il riferimento ?!?!?

ud = zeros(1, length(t));
zd1 = zeros(2, length(t));

% TODO: Calcolo le matrici K, Kg e g per il controllo LQT (CAMBIA con quello da slide)
for i = 1 : length(t)-1
    ud(:,i) = K(:,i) * zd1(:,i) + Kg(:,i) * g(i+1);
    
    zd1(:,i+1) = A * zd1(:, i) + B * ud(i);
end

% Ora combino tutte le componenti 

% Evoluzione del sistema controllato -> SBAGLIATO
% zc = zeros(2, length(t));
% zc(:,1) = [0; 291];
% for k = 1 : length(t)-1
%     
%     % Calcolo del controllo
%     q(k) = K(:,:,k) * (zc(:,k) - z_d2(:,k)) + Kg(:,:,k) * g(:,k+1);
%     
%     if (q(k) > 1000)
%         q(k) = 1000;
%     elseif (q(k) < -1000)
%         q(k) = -1000;
%     end
% 
%     % Evoluzione del sistema
%     zc(:,k+1) = A * zc(:,k) + B * q(:,k) + mi(:,k) + omega(:,k);
%     y(k) = C * zc(:,k) + csi_f(k);
% 
% end
% 
% %stairs(t, (zc(1,:) + z_hat(1,:))');
% plot(t, (zc(1,:) + z_hat(1,:))', t, z_hat(1,:)');
% legend('theta_c controlled','z_hat');