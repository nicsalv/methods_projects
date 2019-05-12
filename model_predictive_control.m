function [u, z] = model_predictive_control(A, B, C, omega, csi_f, M, MT, N, K_kalm, z_hat, z_d2, t, timeWindow)
    
    % Inizializzo l'uscita
    y = zeros(1, length(t));

    % Inizializzo i vettori dei contolli
    u_d = zeros(1, length(t));
    u_s = zeros(1, length(t));

    % inizializzo i vettori degli stati
    z_s = zeros(2, length(t));
    z_estimated = zeros(2, length(t));
    z_d1 = zeros(2, length(t));
    
    % Applichiamo la finestra temporale solo all'LQG -> ovvero la parte
    % stocastica del problema
    
    time = 1:timeWindow;
    
    % Determinazione della matrice K per il controllo ottimo
    [P, K] = riccati_P_K(A, B, M, MT, N, time);
    
    % Stima e conti all'istante iniziale
    y(1) = C * z_s(:,1) + csi_f(1);
    z_estimated(:,1) = z_s(:,1) + K_kalm(:,:,1) * (y(1) - C*z_s(:,1));
    u_s(1) = K(:,:,1) * z_estimated(:,1);
    z_s(:,2) = A * z_s(:,1) + B * u_s(1) + omega(:,1);
    
    % Indici della finestra temporale    
    for i = 2:length(t)-1        
        % Misurazione dell'uscita del sistema
        y(i) = C * z_s(:,i) + csi_f(i);

        % Stima dello stato mediante Kalman
        sys_kalm = A * z_estimated(:,i-1) + B * u_s(:,i-1);
        z_estimated(:,i) = sys_kalm + K_kalm(:,:,i) * (y(i) - C * sys_kalm); 

        [P, K] = riccati_P_K(A, B, M, MT, N, time);
        
        % Calcolo del controllo
        u_s(i) = K(:,:,1) * z_estimated(:,i);

        % Evoluzione del sistema affetto dal controllo ottimo
        z_s(:,i+1) = A * z_s(:,i) + B * u_s(i) + omega(:,i);
        
    end

    
    % Parte deterministica -> rimane inveriata rispetto alla finestra
    % temporale
    r = z_hat - z_d2;

   [K, Kg, g] = Riccati_nonStandard_LQG(M, N, MT, A, B, r, t);

    for i = 1 : length(t)-1
%         if ((length(t) - i) > timeWindow)
%             time = i:(i+timeWindow-1);
%         else
%             time = i:length(t);
%         end
%         rTemp(:,time) = r(:,time);

%         [K, Kg, g] = Riccati_nonStandard_LQG(M, N, MT, A, B, rTemp, time);
        
        % Calcolo del controllo ottimo
        u_d(:,i) = K(:,:,i) * z_d1(:,i) + Kg(:,:,i) * g(:,i+1);
%         u_d(:,i) = K(:,:,1) * z_d1(:,i) + Kg(:,:,1) * g(:,2);

        % Evoluzione del sistema affetto dal controllo
        z_d1(:,i+1) = A * z_d1(:, i) + B * u_d(i);
    end

    % Ricostruzione del controllo e dello completo
    u = u_d + u_s;
    z = z_s + z_d1 + z_d2;

end