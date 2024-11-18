clear all
close all


%%%%%%%%%%%%%%%%%%%%%%%%
% small matrix sizes
% (good for accuracy tests)
%%%%%%%%%%%%%%%%%%%%%%%%
X = Matrix_Helmholtz(4000,8000);

%X = Matrix_Gaussian_exp(1000);

%X = Matrix_SNN(1000);

%X = Matrix_GMM(2000, 500);

%load mnist.mat
%load cifar10.mat

%X = A(:, randperm(60000, 1000));
%X = X ./ vecnorm(X);

%%%
% Every column of matrix A contains an image from 
%   mnist (length 28x28=782) or cifar10 (length 32x32x3=3072).
% Matrix A has 60,000 columns containing all the training
%   and testing images.
% mnist can be obtained from, e.g., https://github.com/sunsided/mnist-matlab
% cifar10 can be obtained from, e.g., https://www.mathworks.com/matlabcentral/fileexchange/62990-deep-learning-tutorial-series
%%%


%%%%%%%%%%%%%%%%%%%%%%%%
% large matrix sizes
% (good for speed tests)
%%%%%%%%%%%%%%%%%%%%%%%%

%load mnist.mat
%load cifar10.mat
%X = A(:, randperm(60000, 60000));
%X = X ./ vecnorm(X);


%X = Matrix_Gaussian_exp(20000);

%X = Matrix_SNN(20000);

%X = Matrix_GMM(100000, 1000);



% use a gpu if available
if gpuDeviceCount("available")
  X = gpuArray(X);
end


[m,n] = size(X);
fprintf('Matrix size: %d x %d\n', m, n);




%% RBRP

maxNumCompThreads(1);

% block size for RBRP
blk = 64;


tic
[Q2, R2, p2, res_rbrp, time_rbrp, rank_rbrp] = RBRP_MGS(X, blk, 'random');
t_rbrp = toc;


% pick the first few ranks for comparison with other methods
m = 7; 
ranks = rank_rbrp(1:m);


err_rbrp = zeros(1,m);
t_rbrp_id = zeros(1,m);
for i=1:m
    k = rank_rbrp(i);
    tic;
    T = R2(1:k,1:k) \ R2(1:k,k+1:end);
    t_rbrp_id(i) = toc;
    if isgpuarray(X); assert(isgpuarray(T)); end

    sk = p2(1:k);
    rd = p2(k+1:end);
    err_rbrp(i) = norm(X(:,rd)-X(:,sk)*T, 'fro') / norm(X, 'fro');
end
err_rbrp = err_rbrp.^2;

time_rbrp = time_rbrp(1:m) + t_rbrp_id;





%% SVD 


maxNumCompThreads('auto');

if isgpuarray(X)
  display('SVD on GPU')
else
  fprintf('SVD on CPU with %d threads\n', maxNumCompThreads)
end

tic
[~, S, ~] = svd(X,"econ","vector");
t_svd = toc;
fprintf('SVD time: %.2e s\n', t_svd)

err_svd = zeros(1,length(ranks));
for i=1:length(ranks)
    err_svd(i) = norm(S(ranks(i)+1:end)) / norm(S);
end
err_svd = err_svd.^2; % squared




%% CPQR

maxNumCompThreads(1);

if isgpuarray(X)
  display('CPQR on GPU')
else
  fprintf('CPQR on CPU with %d threads\n', maxNumCompThreads)
end


tic; 
[~, R1, p1] = qr(X,'econ','vector'); 
t_cpqr = toc;


% error
res_cpqr = zeros(1,length(ranks));
for i=1:length(ranks)
    res_cpqr(i) = norm(R1(ranks(i)+1:end,ranks(i)+1:end),"fro") / norm(R1,"fro");
end
res_cpqr = res_cpqr.^2;

% time for interpolation matrix
err_cpqr = zeros(1,length(ranks));
t_cpqr_id = zeros(1,length(ranks));
for i=1:length(ranks)
    tic;
    T = R1(1:ranks(i),1:ranks(i)) \ R1(1:ranks(i),ranks(i)+1:end);
    t_cpqr_id(i) = toc;
    if isgpuarray(X); assert(isgpuarray(T)); end

    sk = p1(1:ranks(i));
    rd = p1(ranks(i)+1:end);
    err_cpqr(i) = norm(X(:,rd)-X(:,sk)*T, 'fro') / norm(X, 'fro');    
end
err_cpqr = err_cpqr.^2;

% work per step for Householder QR
[xm, xn] = size(X);
p = min(xm,xn);
work = (xm:-1:xm-p+1) .* (xn:-1:xn-p+1);
work = cumsum(work);
time_cpqr = t_cpqr/work(end)*work(ranks);
time_cpqr = time_cpqr + t_cpqr_id;


%% randCPQR

err_randCPQR = zeros(1,length(ranks));
time_randCPQR = zeros(1,length(ranks));
for i=1:length(ranks)
    tr = tic;
    [sk, rd, T] = randCPQR(X, ranks(i));
    time_randCPQR(i) = toc(tr);

    if isgpuarray(X); assert(isgpuarray(T)); end
    err_randCPQR(i) = norm(X(:,rd)-X(:,sk)*T, 'fro')^2 / norm(X,'fro')^2;
end


err_randCPQR_OS = zeros(1,length(ranks));
time_randCPQR_OS = zeros(1,length(ranks));
for i=1:length(ranks)
    tr = tic;
    [sk, rd, T] = randCPQR_OS(X, ranks(i));
    time_randCPQR_OS(i) = toc(tr);

    if isgpuarray(X); assert(isgpuarray(T)); end
    err_randCPQR_OS(i) = norm(X(:,rd)-X(:,sk)*T, 'fro')^2 / norm(X,'fro')^2;
end



%% randLUPP

err_randLUPP = zeros(1,length(ranks));
time_randLUPP = zeros(1,length(ranks));
for i=1:length(ranks)
    tr = tic;
    [sk, rd, T] = randLUPP(X, ranks(i));
    time_randLUPP(i) = toc(tr);

    if isgpuarray(X); assert(isgpuarray(T)); end
    err_randLUPP(i) = norm(X(:,rd)-X(:,sk)*T, 'fro')^2 / norm(X,'fro')^2;
end


err_randLUPP_OS = zeros(1,length(ranks));
time_randLUPP_OS = zeros(1,length(ranks));
for i=1:length(ranks)
    tr = tic;
    [sk, rd, T] = randLUPP_OS(X, ranks(i));
    time_randLUPP_OS(i) = toc(tr);

    if isgpuarray(X); assert(isgpuarray(T)); end
    err_randLUPP_OS(i) = norm(X(:,rd)-X(:,sk)*T, 'fro')^2 / norm(X,'fro')^2;
end


%% SqNorm


err_sqn = zeros(1,length(ranks));
time_sqn = zeros(1,length(ranks));
for i=1:length(ranks)
    tr = tic;
    [sk, rd, T] = SqNorm(X, ranks(i));
    time_sqn(i) = toc(tr);

    if isgpuarray(X); assert(isgpuarray(T)); end
    err_sqn(i) = norm(X(:,rd)-X(:,sk)*T, 'fro')^2 / norm(X,'fro')^2;
end




%% plot

figure(1)
plt=semilogy(ranks, err_svd, 'ok-',...
    ranks,res_cpqr,'ro-',...
    ranks,err_randCPQR_OS, 'co-', ...
    ranks,err_randLUPP_OS, 'bo-', ...
    ranks,err_sqn, 'yo-', ...
    rank_rbrp(1:m),res_rbrp(1:m),'go-');


grid on
set(plt,'LineWidth',3,'MarkerSize',15);
ax = gca;
set(ax,'fontsize',22);
xlabel('rank');
ylabel('squared error');
legend('SVD', 'CPQR', 'randCPQR-OS', 'randLUPP-OS', 'SqNorm', 'New','Location', 'best')
saveas(gcf,'fig_error','epsc')

%%

figure(2)
plt2=plot(ranks, time_cpqr, 'ro-',...
    ranks,time_randCPQR_OS, 'co-', ...
    ranks,time_randLUPP_OS, 'bo-', ...
    ranks,time_sqn, 'yo-', ...
    rank_rbrp(1:m), time_rbrp, 'go-');

grid on
set(plt2,'LineWidth',3,'MarkerSize',15);
ax = gca;
set(ax,'fontsize',22);
xlabel('rank');
ylabel('Runtime (sec)');
legend('CPQR', 'randCPQR-OS', 'randLUPP-OS', 'SqNorm', 'New','Location', 'Best')
saveas(gcf,'fig_speed','epsc')


%%

fid = fopen('results.m','w');
fprintf(fid, 'SVD\n');
fprintf(fid, '%1.3e\t', err_svd);
fprintf(fid, '\n');

fprintf(fid, 'ranks\n');
fprintf(fid, '%5d\t', ranks);
fprintf(fid, '\n');

fprintf(fid, 'RBRP\n');
fprintf(fid, '%1.3e\t', res_rbrp);
fprintf(fid, '\n');
fprintf(fid, '%1.3e\t', time_rbrp);
fprintf(fid, '\n');

fprintf(fid, 'CPQR\n');
fprintf(fid, '%1.3e\t', res_cpqr);
fprintf(fid, '\n');
fprintf(fid, '%1.3e\t', time_cpqr);
fprintf(fid, '\n');

fprintf(fid, 'randomized CPQR\n');
fprintf(fid, '%1.3e\t', err_randCPQR);
fprintf(fid, '\n');
fprintf(fid, '%1.3e\t', time_randCPQR);
fprintf(fid, '\n');

fprintf(fid, 'randomized CPQR-OS\n');
fprintf(fid, '%1.3e\t', err_randCPQR_OS);
fprintf(fid, '\n');
fprintf(fid, '%1.3e\t', time_randCPQR_OS);
fprintf(fid, '\n');

fprintf(fid, 'randomized LUPP\n');
fprintf(fid, '%1.3e\t', err_randLUPP);
fprintf(fid, '\n');
fprintf(fid, '%1.3e\t', time_randLUPP);
fprintf(fid, '\n');

fprintf(fid, 'randomized LUPP-OS\n');
fprintf(fid, '%1.3e\t', err_randLUPP_OS);
fprintf(fid, '\n');
fprintf(fid, '%1.3e\t', time_randLUPP_OS);
fprintf(fid, '\n');

fprintf(fid, 'randomized SqNorm\n');
fprintf(fid, '%1.3e\t', err_sqn);
fprintf(fid, '\n');
fprintf(fid, '%1.3e\t', time_sqn);
fprintf(fid, '\n');
fclose(fid);
