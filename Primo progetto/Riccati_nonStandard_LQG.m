function [K, Kg, g] = Riccati_nonStandard_LQG(M, N, MT, A, B, r, t)
    % r: segnale su cui effettuare il tracking
    % t: asse dei tempi
    
    % Define global variable
    dim_time = length(t);
    dim_state = size(A,1);
    dim_control = size(N,1);
    I = eye(dim_state);
    
    % Compute the P matrix (Riccati matrix)
    P = zeros(dim_state, dim_state, dim_time);
    P(:,:,end) = MT;
    for t = dim_time-1 : -1 : 1
        P(:,:,t) = M + A' * P(:,:,t+1) *...
            ((I + (B*(N\B')) * P(:,:,t+1)) \ A);
    end
    
    % Compute g
    g = zeros(dim_state, dim_time);
    g(:,end) = I * MT * r(:,end); % La C non serve: l'uscita è fittizia.
    E = B / N * B';
    for t = dim_time-1 : -1 : 1
        g(:,t) = A' * (I - inv((inv(P(:,:,t+1)) + E)) * E) * g(:,t+1) +...
            I * M * r(:,t);
    end
    
    % Compute Kg matrix
    Kg = zeros(dim_control, dim_state, dim_time - 1);
    for t = 1 : dim_time-1
        Kg(:,:,t) = (N + B' * P(:,:,t+1) * B) \ B';
    end
    
    % Compute K matrix
    K = zeros(dim_control, dim_state, dim_time - 1);
    for t = 1 : dim_time-1
        K(:,:,t) = -(N + B' * P(:,:,t+1) * B) \ (B' * P(:,:,t+1) * A);
    end
end