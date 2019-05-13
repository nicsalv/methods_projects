function [u, z] = standard_control(A, B, C, omega, csi_f, M, MT, N, K_kalm, z_hat, z_d2, t)
    % Controllo LQG/LQT come indicato nel paper: parte stocastica del sistema.
    % Si risolve mediante un LQG standard,
    % trascurando il disturbo deterministico.
    
    % Uscita del sistema
    y = zeros(1, length(t));

    % Parte stocastica del controllo
    u_s = zeros(1, length(t));

    % Parte stocastica dello stato
    z_s = zeros(2, length(t));
    % z_s(:,1) = [10 -15]'; % Utile per verificare LQG
    z_estimated = zeros(2, length(t));

    % Determinazione della matrice K per il controllo ottimo
    [P, K] = riccati_P_K(A, B, M, MT, N, t);

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
    
    % Utile per apprezzare le differenze tra questo e il receding horizon.
%     subplot(2,1,1);
%     stairs(t, [z_s(1,:)', z_estimated(1,:)']);
%     legend('z_s', 'z_estimated');

%     subplot(2,1,2);
%     stairs(t, y');
%     legend('y');


    % Parte deterministica del sistema.
    % Si tiene conto del disturbo deterministico theta_a.
    % Questo è un problema LQT

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

end