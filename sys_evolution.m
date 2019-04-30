%%%
% Definizione ed evoluzione del sistema libero
%%%

clear;

%% Costanti numeriche
cc = 2.15;
mc = 2;
Cc = cc * mc;

cf = 1;
mf = 1.3;
Cf = cf * mf;

sc = .095; % squared meteres
k_bar_fc = 100;
kfc = k_bar_fc * sc;

sf = 6;
k_bar_af = 1.22;
kaf = k_bar_af * sf;

% Caratteristiche distrubi
% Deviazioni standard (radici della Varianza)
c_std = 1 / Cc;
f_std = 1 / Cf;

%% Evoluzione del sistema
% Parametri di funzionamento
T = 24;
sampling_time = 1;
t = 0 : sampling_time : T;

% Distrubi
% Inizializzazione della variabile casuale.
% In questo modo l'esperimento è ripetibile
rng(0, 'twister');

omega_c = c_std .* randn(1, length(t));
omega_f = f_std .* randn(1, length(t));
csi_f = .7071 .* randn(1, length(t));

% Raccolta dei disturbi sullo stato in un unica variabile
omega = [omega_c; omega_f];

% Verifica delle caratteristiche
stats = [mean(omega_c) var(omega_c); ...
    mean(omega_f) var(omega_f); ...
    mean(csi_f) var(csi_f)]

theta_a = [291 291 293 295 295 295 297 297 298 298 299 301 303 303 303 ...
    305 307 308 309 306 303 299 297 295 290];

xq = 0 : .1 : T;
vq1 = interp1(t,theta_a,xq);

plot(t,theta_a,'o',xq,vq1,':.');
