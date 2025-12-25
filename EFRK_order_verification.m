%% Interested readers could change the values of 'RK_flag' to verify the EFRK method 
%%  based on a given Butcher tableau. 
RK_flag = 32; % Verify order conditions of RK(4, 4)
switch RK_flag
    case 1 % RK(1, 1)
        A = [0; 1]; order = 1;
    case 2 % RK(2, 2), Heun's second-order scheme
        A = [0 0; 1 0; 1/2 1/2]; order = 2;
    case 3 % RK(3, 3) in Eq. (18)
        A = [0 0 0; 2/3 0 0; 2/9 4/9 0; 1/4 3/16 9/16]; order = 3;
    case 32 % RK(3,3) in Eq. (41)
       A = [0 0 0; 1 0 0; 1/4 1/4 0; 1/6 1/6 2/3]; order = 3;
    case 4 % RK(4, 4)
        A = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0; 1/6 1/3 1/3 1/6]; order = 4;
    case 5 % RK(5, 4) 
        alpha = [0 0 0 0 0;
            1, 0 0 0 0;
            0.444370493651235 0.555629506348765 0 0 0;
            0.620101851488403 0 0.379898148511597 0 0;
            0.178079954393132 0 0 0.821920045606868 0;
            0 0 0.517231671970585 0.096059710526147 0.386708617503268];
        beta = [0 0 0 0 0;
            0.39175222657189  0 0 0 0;
            0 0.368410593050371 0 0 0;
            0 0 0.251891774271694 0 0;
            0 0 0 0.54497475022852  0;
            0 0 0 0.063692468666290 0.226007483236906];
        [A, b] = shuosher2butcher(alpha, beta);
        A = [A; b']; order = 4;
    otherwise
        stage = 4; order = 4;
        A = sym('A', [stage+1, stage]);
        for i = 1:stage
            for j = i:stage
                A(i,j) = 0;
            end
        end
end
syms z; % Denote $z:= \tau L_\kappa$
stage = size(A, 2);
A_hat = sym(zeros(size(A)));
c_hat = sym(zeros(stage+1,1));
psi   = sym(zeros(stage+1,1));

% Calculate coefficients of EFRK method
psi(1) = 1;
for i = 2:stage+1
    psi(i) = 1; c_hat(i) = 0;
    for j = 1:i-1
        psi(i) = psi(i) + z * A(i,j) * psi(j);
    end
    for j = 1:i-1
        A_hat(i,j) = A(i,j) * psi(j) / psi(i);
        for k = j+1:i-1
            A_hat(i,j) = A_hat(i,j) + z * A(i,k) * psi(k) / psi(i) * A_hat(k,j);
        end
        c_hat(i) = c_hat(i) + A_hat(i,j);
    end
end

% Verification of order conditions
b_hat = A_hat(stage+1,:);
A_hat = A_hat(1:stage,:);
c_hat = c_hat(1:stage);
c     = sum(A(1:stage,:), 2);
if order >= 1
    fprintf('Leading Error z^%d:\n', order);
    taylor(sum(b_hat) - 1, z, 0, 'order', order+2)
end
if order >= 2
    fprintf('Leading Error z^%d:\n', order-1);
    for i = 0:1
        taylor(b_hat*(c_hat.^i.*c.^(1-i)) - 1/2, z, 0, 'order', order+1)
    end
end
if order >= 3
    fprintf('Leading Error z^%d:\n', order-2);
    for i = 0:2
        taylor(b_hat*(c_hat.^i.*c.^(2-i)) - 1/3, z, 0, 'order', order)
    end
    for i = 0:1
        taylor(b_hat*A_hat*(c_hat.^i.*c.^(1-i)) - 1/6, z, 0, 'order', order)
    end
end
if order >= 4
    fprintf('Leading Error z^%d:\n', order-3);
    for i = 0:3
        taylor(b_hat*(c_hat.^i.*c.^(3-i)) - 1/4, z, 0, 'order', order-1)
    end
    for i = 0:1
        for j = 0:1
            taylor(b_hat*((c_hat.^i.*c.^(1-i)).*(A_hat*(c_hat.^j.*c.^(1-j)))) - 1/8, z, 0, 'order', order-1)
        end
    end
    for i = 0:2
        taylor(b_hat*A_hat*(c_hat.^i.*c.^(2-i)) - 1/12, z, 0, 'order', order-1)
    end
    for i = 0:1
        taylor(b_hat*A_hat^2*(c_hat.^i.*c.^(1-i)) - 1/24, z, 0, 'order', order-1)
    end
end
