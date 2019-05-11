function [u, z] = model_predictive_control(A, B, C, omega, csi_f, M, MT, N, K_kalm, z_hat, z_d2, t, timeWindow)

    % Indici della finestra temporale
    windowStart = 1;
    windowEnd = timeWindow;
    
    % Inizializzo l'uscita
    y = zeros(1, length(t));

    % Inizializzo i vettori dei contolli
    u_d = zeros(1, length(t));
    u_s = zeros(1, length(t));

    % inizializzo i vettori degli stati
    z_s = zeros(2, length(t));
    z_estimated = zeros(2, length(t));
    z_d1 = zeros(2, length(t));

    time = windowStart:windowEnd;
    target = z_hat(:, time);
    
    % istante iniziale del sistema: parte stocastica
    [P, K_stoc] = riccati_P_K(A, B, M, MT, N, time);
    y(1) = C * z_s(:,1) + csi_f(1);
    z_estimated(:,1) = z_s(:,1) + K_kalm(:,:,1) * (y(1) - C*z_s(:,1));
    u_s(1) = K_stoc(:,:,1) * z_estimated(:,1);
    z_s(:,2) = A * z_s(:,1) + B * u_s(1) + omega(:,1);

    % istante inziale del sistema: parte deterministica
    r = target - z_d2(:, time);
    [K_det, Kg, g] = Riccati_nonStandard_LQG(M, N, MT, A, B, r, time);
    u_d(:,1) = K_det(:,:,1) * z_d1(:,1) + Kg(:,:,1) * g(:,2);
    z_d1(:,2) = A * z_d1(:, 1) + B * u_d(1);

    windowStart = windowStart + 1;
    windowEnd = windowEnd + 1;

    % TODO: check condition
    while windowStart < length(t)

        % Array con valori desirati nella finestra temporale considerata
        time = windowStart:windowEnd;
        target = z_hat(:, time);

        % Parte stocastica 

        % Calcolo il controllo della riccati per l'LQG -> parte stocastica
        [P, K_stoc] = riccati_P_K(A, B, M, MT, N, time);

        % Calcolo la componente stocastica    
        y(windowStart) = C * z_s(:,windowStart) + csi_f(windowStart);

        % Stima dello stato mediante Kalman
        sys_kalm = A * z_estimated(:,windowStart-1) + B * u_s(:,windowStart-1);
        z_estimated(:,windowStart) = sys_kalm + K_kalm(:,:,windowStart) * (y(windowStart) - C * sys_kalm); 

        % Calcolo del controllo
        u_s(windowStart) = K_stoc(:,:,1) * z_estimated(:,windowStart);

        % Evoluzione del sistema affetto dal controllo ottimo
        z_s(:,windowStart+1) = A * z_s(:,windowStart) + B * u_s(windowStart) + omega(:,windowStart);

        % Parte deterministica
        r = target - z_d2(:, time);
        [K_det, Kg, g] = Riccati_nonStandard_LQG(M, N, MT, A, B, r, time);

        u_d(:,windowStart) = K_det(:,:,1) * z_d1(:,windowStart) + Kg(:,:,1) * g(:,2);

        % Evoluzione del sistema affetto dal controllo
        z_d1(:,windowStart+1) = A * z_d1(:, windowStart) + B * u_d(windowStart);

        % Aggiorno il valore degli indici
        windowStart = windowStart + 1;
        if windowEnd < length(t)
            windowEnd = windowEnd + 1;
        end

    end

    % Ricostruzione del controllo e dello completo
    u = u_d + u_s; 
    z = z_s + z_d1 + z_d2;

end