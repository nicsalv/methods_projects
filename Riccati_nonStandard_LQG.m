function [K, Kg, zd2, g] = Riccati_nonStandard_LQG(M, N, MT, eta, z_hat, A, B, z0, T)
    % Compute K_t and Kg_t matries for optimal control of non standard LQG.
    % The cost function is like : J = sum_t( c(zt, ut) ) + zT' MT zT
    % with c(zt, ut) = zt' M zt + ut' N ut
    % eta represents the deterministic noise 
    % z_hat is the optimal working level 
    % In this formulation : mu = eta + (A - I) z_hat
    % z0 is the initial state
    % T is vector of sampling time
    
    % NB: Controlla bene i calcoli con le matrici
    
    % Define global variable
    dim = length(T);
    dim_state = size(z0,1);
    dim_control = size(N,1);
    I = eye(dim_state);
    
    % Compute the mu vector
    mu = eta + (A - I) * z_hat;
    
    % Compute the zd2 vector
    zd2 = zeros(dim_state, dim);
    zd2(:,1) = z0;
    
    for t = 2 : dim
        zd2(:, t) = A * zd2(:, t-1) + mu(t-1);
    end
    
    % Compute the P matrix (Riccati matrix)
    P = zeros(dim_state, dim_state , dim);
    P(:, :, end) = MT;
    
    for t = (dim-1) : -1 : 1
        P(:, :, t) = M + A' * P(:, :, t+1) * ((I + (B * (N \ B')) * P(:, :, t+1)) \ A);
    end
    
    % Compute the g vector
    g = zeros(dim_state, dim);
    g(:,end) = MT * zd2(:, end);
    
    for t = (dim-1) : -1 : 1
        g(:, t) = (A' - (P(:, :, t+1) * ((I + (B * (N \ B'))) \ P(:, :, t+1))) * (B * (N \ B'))) * g(:, t+1) - M * zd2(:, t);
    end
    
    % Compute Kg matrix
    Kg = zeros(dim_control, dim_state, dim - 1);
    for t = 1 : dim-1
        Kg(:, :, t) = (N + B' * P(:, :, t+1) * B) \ B';
    end
    
    % Compute the K matrix
    K = zeros(dim_control, dim_state, dim - 1);
    for t = 1 : dim-1
        K(:, :, t) = -(N + B' * P(:, :, t+1) * B) \ (B' * P(:, :, t+1) * A);
    end
    
end