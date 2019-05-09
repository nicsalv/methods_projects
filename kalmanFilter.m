function [sigma, k] = kalmanFilter(A,C,Q_varianza,R_varianza,SIGMA_varianza,T)
%Q_varianza della normale del rumore sullo stato x
%R_varianza della normale del rumore sulla misura y,
%SIGMA_varianza della normale dello stato all'istante 0

    dim = length(T);
    dim_state = size(A,1);
    dim_exit = size(C,1);
    
    k = zeros(dim_state, dim_exit, dim);
    sigma = zeros(dim_state, dim_state, dim);
    
    sigma(:,:,1) = inv(inv(SIGMA_varianza) + (C'/R_varianza)*C);

    for i = 1:dim
        k(:,:,i) = sigma(:,:,i)*(C'/R_varianza);
        sigma(:,:,i+1) = inv(inv(A*sigma(:,:,i)*A'+Q_varianza)+(C'/R_varianza)*C);
    end


end
