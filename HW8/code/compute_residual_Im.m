function residual = compute_residual_Im(Qh, Q, Fhp, Fhm, area, dx, dt, x, Jp, Jm, res_jmax)
    N = size(Qh, 2)-2;
    Dq_pred = -dt/dx*(Fhp(:,2:end-1)-Fhp(:,1:end-2)+Fhm(:,3:end)-Fhm(:,2:end-1));
    Dq_pred(2,:) = Dq_pred(2,:) + dt/dx*(Q(3,2:end-1) .* ...
                   abs((0.5*(area(3:end)+area(2:end-1))-0.5*(area(2:end-1)+area(1:end-2)))));

    SJ = zeros(N,3,3);
    SJ(:,2,3) = abs((0.5*(area(3:end)+area(2:end-1))-0.5*(area(2:end-1)+area(1:end-2))));
    SJ(:,2,3) = abs((0.5*(area(3:end)+area(2:end-1))-0.5*(area(2:end-1)+area(1:end-2))));
    %gamma = 1.4;
    %SJ(:,2,3) = (gamma-1)*(Qh(3,2:end-1)-0.5*Qh(2,2:end-1).^2./Qh(1,2:end-1)) .* ...
    %            (abs((0.5*(area(3:end)+area(2:end-1))-0.5*(area(2:end-1)+area(1:end-2)))));
    %disp(max(SJ(:,2,3)));
    I = permute(repmat(eye(3),[1,1,N]), [3 1 2]);

    A = -dt/dx * Jp(1:end-2,:,:);
    B = I + dt/dx*(Jp(2:end-1,:,:)-Jm(2:end-1,:,:)) - dt*SJ;
    C = dt/dx * Jm(3:end,:,:);
    D = Dq_pred;
    D(:,1) = D(:,1) - squeeze(A(1,:,:))*zeros(3,1); % Fixed left boundary
    D(:,end) = D(:,end) - squeeze(C(end,:,:))*res_jmax; % Right compability relation

    D = D(:);
    M = block_tridiagonal_matrix_vectorized(A, B, C);

    residual = M \ D;
    residual = reshape(residual, 3, []);
    residual = [zeros(3,1), residual, zeros(3,1)];
end

%%
function M = block_tridiagonal_matrix_vectorized(A, B, C)
    % A: lower diagonal blocks (N×3×3), starts at block row 2
    % B: main diagonal blocks (N×3×3), starts at block row 1
    % C: upper diagonal blocks (N×3×3), starts at block row 1
    
    N = size(B, 1);
    block_size = size(B, 2);  % Assuming 3×3 blocks
    matrix_size = N * block_size;
    
    % Generate local indices for all blocks
    [i_local, j_local] = ndgrid(1:block_size, 1:block_size);
    i_local = i_local(:);
    j_local = j_local(:);
    
    % --- Main diagonal (B) ---
    block_indices = kron((0:N-1)', ones(block_size^2, 1));
    i_main = block_indices*block_size + repmat(i_local, N, 1);
    j_main = block_indices*block_size + repmat(j_local, N, 1);
    v_main = reshape(permute(B, [2, 3, 1]), [], 1);
    
    % --- Upper diagonal (C) ---
    block_indices = kron((0:N-2)', ones(block_size^2, 1));
    i_upper = block_indices*block_size + repmat(i_local, N-1, 1);
    j_upper = (block_indices+1)*block_size + repmat(j_local, N-1, 1);
    v_upper = reshape(permute(C(1:N-1,:,:), [2, 3, 1]), [], 1);
    
    % --- Lower diagonal (A) ---
    block_indices = kron((1:N-1)', ones(block_size^2, 1));
    i_lower = block_indices*block_size + repmat(i_local, N-1, 1);
    j_lower = (block_indices-1)*block_size + repmat(j_local, N-1, 1);
    v_lower = reshape(permute(A(2:N,:,:), [2, 3, 1]), [], 1);
    
    % Combine all indices
    I = [i_main; i_upper; i_lower];
    J = [j_main; j_upper; j_lower];
    V = [v_main; v_upper; v_lower];
    
    % Create sparse matrix
    M = sparse(I, J, V, matrix_size, matrix_size);
end