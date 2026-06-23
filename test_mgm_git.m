%% Clear workspace
 clc
clear 
close all
%% Create example

options = struct();
options.angles = linspace(0,179,180);
% options.phantomImage = image; 
[A,b,x_true,info]=PRtomo(256,options);
% figure, imshow(reshape(x_true,256,256),[])
bsize = info.bSize;
xsize = info.xSize;
rng(5)
e=randn(size(b));
noiseLev = 10;
e=e/norm(e)*norm(b)*noiseLev/100;
b=b+e;
%% Set parameters
maxit=50;
delta=norm(e);
tau=1.01;
err = zeros(maxit,1);

%% Initialize P, R and A
lev = 7;
Ahat = cell(lev+2,1);
Ahat{1} = A;

%% LSQR
tic
x_lsqr_dp=lsqrDP(A,b,tau,delta,300);
toc

it_lsqr = 25;
x_lsqr = lsqr(A,b,it_lsqr);
DP_l = 0;
err_lsqr = zeros(it_lsqr,1);
for i = 1:it_lsqr
    rr=b-A*x_lsqr(:,i); norm(rr);
    if norm(rr)<delta*1.01 && DP_l == 0
    DP_l = i;
    % break
    end
    x_nn = max(x_lsqr,0);
    err_lsqr(i) = norm(x_true-x_nn(:,i))/norm(x_true);
end
tic
x_lsqr = lsqr(A,b,DP_l);
toc

figure
pl = plot(err_lsqr,'LineWidth',1.5,'Color','r');
hold on
[m,I] = min(err_lsqr);
disp(['Err lsqr DP: ',num2str(norm(x_true-x_lsqr(:,DP_l))/norm(x_true))])
disp(['Err lsqr min: ',num2str(m)])
plot(I,m,'-o','Color','r','LineWidth',1.5)
plot(DP_l,err_lsqr(DP_l),'-x','Color','r','LineWidth',1.5)
%%
Min = zeros(4,2);
dp = zeros(4,2);
for k = 1:4
%% Create restrictors and projector
[P,R] = projectors(xsize,bsize,lev,k);
%% Create projected matrix
for i=0:lev
Ahat{i+2}=R{i+1}*Ahat{i+1}*P{i+1};
end

%% Multigrid iterations
tic
x = zeros(xsize(1)*xsize(2),1);
DP = 0;
for i = 1:maxit
    r = b - Ahat{1}*x;
    for j = 1:(lev+1)
        r = R{j}*r;
    end
    h = pinv(full(Ahat{lev+2}))*r;
    for j = (lev+1):-1:1
        h = P{j}*h;
    end
    x = x + h;
    % Post-smoother
     h = lsqr(A,b-A*x,2);
     h = h(:,2);
    x = x + h;
    x = max(x,0); % non-negative cone
    err(i) = norm(x_true-x)/norm(x_true);
    
    % DP 
    rr=b-A*x;
    if norm(rr)<delta*1.01 && DP == 0
    DP = i;
    xDP = x;
    Time_MGM = toc
    % break
    end
    % hh=lsqr(A,rr,1);   
end

%% Show results

p(k) = plot(err);
[Min(k,1),Min(k,2)] = min(err);
plot(Min(k,2),Min(k,1),'-o','Color','k')
plot(DP,err(DP),'-x','Color','k')

dp(k,1) = norm(x_true-xDP)/norm(x_true);
dp(k,2) = DP;
disp(['Err MGM DP: ',num2str(dp(k,1))])
disp(['Err MGM min: ',num2str(Min(k,1))])
end
%%
max_sirt = 80;
x_s = sirt('sart',A,b,1:max_sirt);
err_s = zeros(max_sirt,1);
for i = 1:max_sirt
    err_s(i) = norm(x_true-x_s(:,i))/norm(x_true);
end
ps = plot(err_s,'linewidth',2,'Color','g');
legend('SIRT')
legend([pl,p(1),p(2),p(3),p(4),ps],{'LSQR','[0 1 1]','[1 2 1]','[0 1 3 3 1]', '[1 4 6 4 1]', 'SIRT'})
title(['LSQR vs MGMs vs SIRT with ',num2str(noiseLev),'% noise'])
