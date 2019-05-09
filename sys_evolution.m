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
omega = [omega_c; omega_f] .* deltaT; % Rumore sullo stato
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
    mi(:,k) = eta(:,k) + (A - eye(2)) * z_hat(:,k);
    
    % La k-esima colonna di x corrisponde all'istante k*deltaT
    z(:,k+1) = A * z(:,k) + mi(:,k) + omega(:,k);
    y(k) = C * z(:,k) + csi_f(k);
end

% Controllo LQG/LQT come indicato nel paper

% Calcolo del controllo
M = [.1 0; 0 0]; % Si vuole inseguire solo una delle due componenti
N = 1e-12; % Controllo scalare
MT = M;
[K,Kg,z_d2,g] = Riccati_nonStandard_LQG(M, N, MT, mi, A, B, z(:,1), t);

% Evoluzione del sistema controllato
zc = zeros(2, length(t));
zc(:,1) = [0; 291];
for k = 1 : length(t)-1
    
    % Calcolo del controllo
    q(k) = K(:,:,k) * (zc(:,k) - z_d2(:,k)) + Kg(:,:,k) * g(:,k+1);
    
    if (q(k) > 1000)
        q(k) = 1000;
    elseif (q(k) < -1000)
        q(k) = -1000;
    end

    % Evoluzione del sistema
    zc(:,k+1) = A * zc(:,k) + B * q(:,k) + mi(:,k) + omega(:,k);
    y(k) = C * zc(:,k) + csi_f(k);

end

%stairs(t, (zc(1,:) + z_hat(1,:))');
plot(t, (zc(1,:) + z_hat(1,:))', t, z_hat(1,:)');
legend('theta_c controlled','z_hat');