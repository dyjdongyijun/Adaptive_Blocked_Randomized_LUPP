
function [Q, R, perm, res, time, rank] = RBRP_MGS(A, blk, mode)

    if isgpuarray(A)
      fprintf('RBRP on GPU')
    else
      fprintf('RBRP on CPU with %d threads\n', maxNumCompThreads)
    end

[m, n] = size(A);
perm = 1:n;

diags = vecnorm(A);
res = [];
time = [];
rank = [];
nrmA2 = sum(diags); 

p = min(m,n);

if isgpuarray(A)
  Q = zeros(m, p, 'gpuArray');
  R = zeros(p, n, 'gpuArray');
else
  Q = zeros(m, p);
  R = zeros(p, n);
end

i = 1;
while i<p+1
    t_rbrp=tic;
    bs = min([blk, p-i+1, sum(diags(i:n)>0)]);

    % sampling a subset of indices
    if strcmp(mode, 'greedy') == 1
        % Greedy
        [~, ind] = sort(diags(i:n), 'descend');
        J  = ind(1:bs) + i-1; 
    else
        % Sampling
        w = diags(i:n);

        if isgpuarray(A)
            w = gather(w);
            bs = gather(bs);
        end
        J = datasample(i:n, bs, 'Replace', false, 'Weights', w.*w);
    end

    % select a subset of columns
    X = A(:,J);

    % Gram Schmidt
    I1 = 1:i-1;
    %X = X - Q(:,I1) * R(I1,J);

    % CPQR on the panel
    [Q2, R2, p2] = qr(X, 'econ', 'vector');

    % reorder selected indices
    J = J(p2); 


    % filter
    v = sum(abs(R2).^2, 2);
    v = cumsum(flip(v));
    bs = sum(v > v(end)/blk);

    % working set
    I2 = i:i+bs-1;

    % remaining set
    I3 = i+bs:n;

    % update global QR
    J = J(1:bs);
    Q(:,I2)  = Q2(:,1:bs);
    R(I2,I2) = R2(1:bs,1:bs);


    K1 = setdiff(I2, J, 'stable'); % indices out
    K2 = setdiff(J, I2, 'stable'); % indices in


    % swamp columns of perm
    T = perm(J);
    perm(K2) = perm(K1);
    perm(I2) = T;
    
    % swamp columns of R
    T = R(I1,J);
    R(I1,K2) = R(I1,K1);
    R(I1,I2) = T;

    % swamp columns of A
    A(:,K2) = A(:,K1);

    % swamp columns of diags
    diags(K2) = diags(K1);


    % post processing
    R(I2,I3) = Q2(:,1:bs)' * A(:,I3);
    A(:,I3)  = A(:,I3) - Q2(:,1:bs) * R(I2,I3);

    % update column norms
    diags(I3) = square_difference( diags(I3), sqrt(sum(abs(R(I2,I3)).^2, 1)) );
    tmp = toc(t_rbrp);

    time = [time, tmp];
    res  = [res, sum(diags(I3))];
    rank = [rank, bs];

    if sum(diags(I3)>0) == 0
        fprintf('zero residual\n')
        break
    end


    % next block
    i = i + bs;
end

res = res ./ nrmA2;
time = cumsum(time);
rank = cumsum(rank);
end
