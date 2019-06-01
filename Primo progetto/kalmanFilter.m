function [sigma, K] = kalmanFilter(A,C,Q,R,SIGMA,t)
% Q varianza del rumore sullo stato
% R varianza del rumore sull'uscita
% SIGMA varianza della stima iniziale dello stato

    dim_time = length(t);
    dim_state = size(A,1);
    dim_output = size(C,1);
    
    K = zeros(dim_state, dim_output, dim_time);
    sigma = zeros(dim_state, dim_state, dim_time);
    
    sigma(:,:,1) = inv(inv(SIGMA) + (C'/R)*C);
    % C' * inv(R) * C

    for i = 1 : dim_time
        K(:,:,i) = sigma(:,:,i) * C' / R;
        sigma(:,:,i+1) = inv(inv(A*sigma(:,:,i)*A' + Q) + (C' / R * C));
    end
end
