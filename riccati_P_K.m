function [P, K, K_infinito] = riccati_P_K(A, B, Q, Qf, R, T)
    N = length(T)-1;

    %calcolo K a tempo finito
    P(:,:,N+1) = Qf;
    for i=N:-1:2
        P(:,:,i) = Q+A'*P(:,:,i+1)*A-A'*P(:,:,i+1)*B*...
             (inv(R+B'*P(:,:,i+1)*B))*B'*P(:,:,i+1)*A;
    end

    for i=1:N
        K(:,:,i) = -inv(R + B'*P(:,:,i+1)*B)*...
                  B'*P(:,:,i+1)*A;
    end

    %calcolo K a tempo infinito
    P_inf(:,:,N+1) = Qf;
    for i=N:-1:1
        P_inf(:,:,i) = Q+A'*P_inf(:,:,i+1)*A-A'*P_inf(:,:,i+1)*B*...
             (inv(R+B'*P_inf(:,:,i+1)*B))*B'*P_inf(:,:,i+1)*A;
    end

    for i=1:N
        K_inf(:,:,i) = -inv(R + B'*P_inf(:,:,i+1)*B)*...
                  B'*P_inf(:,:,i+1)*A;
    end

    K_infinito = K_inf(:,:,1);
end
