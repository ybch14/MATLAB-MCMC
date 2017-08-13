function [Zeta,std1,std2,std3]=RBM_SAMS(W_h10 ,b_h10 ,a_h10 ,...
    W_h20 ,b_h20 ,a_h20 ,...
    W_h100,b_h100,a_h100,...
    W_h500,b_h500,a_h500)
tic;
numvis=[784,784,784,784];numhid=[10,20,100,500];
m=4;N=500000;
% W is the coupling parameter for four RBMs
W=zeros(numvis(1),numhid(4),m);
W(:,1:numhid(1),1)=W_h10 ;W(:,1:numhid(2),2)=W_h20;
W(:,1:numhid(3),3)=W_h100;W(:,1:numhid(4),4)=W_h500;
% b is the biases of visible units for four RBMs
b=zeros(4,numvis(1));
b(1,:)=b_h10;b(2,:)=b_h20;b(3,:)=b_h100;b(4,:)=b_h500;
% a is the biases of hidden units for four RBMs
a=zeros(4,numhid(4));
a(1,1:numhid(1))=a_h10 ;a(2,1:numhid(2))=a_h20;
a(3,1:numhid(3))=a_h100;a(4,1:numhid(4))=a_h500;
% pai is fixed positive weights
pai=1/m*ones(1,m);Zeta=zeros(1,4);t0=randi([8,400]);
% Gamma is the transtion matrix of the label L
Gamma=[0,1,0,0;0.5,0,0.5,0;0,0.5,0,0.5;0,0,1,0];
L=1;v=zeros(1,numvis(1));

Zeta_all=zeros(m,N);
% start the interate procedure
for t=1:N
    % generate j from Gamma(L,:)
    temp_Gamma=cumsum(Gamma(L,:));
    u=rand;buffer=find(temp_Gamma>=u);
    j=buffer(1);
    % calculate the accept probability of j
    W_j=W(:,1:numhid(j),j);b_j=b(j,:);a_j=a(j,1:numhid(j));
    log_qj=b_j*v'+sum(log(1+exp(v*W_j+a_j)));
    W_L=W(:,1:numhid(L),L);b_L=b(L,:);a_L=a(L,1:numhid(L));
    log_qL=b_L*v'+sum(log(1+exp(v*W_L+a_L)));
    chi=min([1,Gamma(j,L)/Gamma(L,j)*exp(-Zeta(j)+log_qj+Zeta(L)-log_qL)]);
    u=rand;
    L=(u<chi)*j+(u>=chi)*L;
    % generate v_(t+1) using Gibbs sample
    W_L=W(:,1:numhid(L),L);b_L=b(L,:);a_L=a(L,1:numhid(L));
    prob_h_v=sigm(v*W_L+a_L);
    h=prob_h_v>rand(size(prob_h_v));
    prob_v_h=sigm(h*W_L'+b_L);
    v=prob_v_h>rand(size(prob_v_h));
    % update Zeta
    delta=((1:m)==L);gamma_t=t0/max([t0,t]);
    zeta=Zeta+gamma_t*(delta-pai);
    Zeta=zeta-zeta(1);
    Zeta_all(:,t)=Zeta';
end
toc;
std1=sqrt(mean((Zeta_all(2,N-N/100+1:N)-mean(Zeta_all(2,N-N/100+1:N))).^2));
std2=sqrt(mean((Zeta_all(3,N-N/100+1:N)-mean(Zeta_all(3,N-N/100+1:N))).^2));
std3=sqrt(mean((Zeta_all(4,N-N/100+1:N)-mean(Zeta_all(4,N-N/100+1:N))).^2));
figure;
hold on;
plot(1:N,Zeta_all(1,:),'blue');
plot(1:N,Zeta_all(2,:),'red');
plot(1:N,Zeta_all(3,:),'green');
plot(1:N,Zeta_all(4,:),'yellow');
title('SAMC算法结果变化趋势');
xlabel('迭代次数');
ylabel('相对值');
legend('h10模型','h20模型','h100模型','h500模型');
end