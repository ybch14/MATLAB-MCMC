function [logZZ,logZZ_up,logZZ_down]=RBM_RTS(parameter_W,parameter_a,parameter_b,beta,batchdata)
tic;
a_B=parameter_a;
a_A=0*parameter_a;
W_B=parameter_W;
% W_A=0*parameter_W
b_B=parameter_b;
[numcases,numdims,numbatches]=size(batchdata);
count_int=zeros(numdims,1);
for batch=1:numbatches
    xx=sum(batchdata(:,:,batch));
    count_int=count_int+xx';
end
lp=5;
p_int=(count_int+lp*numbatches)/(numcases*numbatches+lp*numbatches);
log_base_rate = log(p_int) - log(1-p_int);
b_A=log_base_rate';

[numvis,numhid]=size(parameter_W);
K=length(beta);
log_ZA=sum(log(1+exp(b_A)))+sum(log(1+exp(a_A)));
Z=ones(1,K);

runtime=200;numsweep=50;N=100;
c_all=zeros(runtime,K);
error=zeros(runtime,1);
for run=1:runtime
    c=zeros(1,K);
    v=rand(1,numvis)>0.5;
    beta_next=beta(randi([1,K]));
    for n=1:N
        for i=1:numsweep
            prob_hB_v=1./(1+exp(-beta_next*(v*W_B+a_B)));
            hB=prob_hB_v>rand(size(prob_hB_v));
            prob_v_h=1./(1+exp(-(1-beta_next)*b_A-beta_next*(hB*W_B'+b_B)));
            v=prob_v_h>rand(size(prob_v_h));
        end
        log_fk=zeros(1,K);
        v_bA=v*b_A';v_bB=v*b_B';Wh=v*W_B+a_B;
        for k=1:K
            log_fk(k)=(1-beta(k))*v_bA+sum(log(1+exp((1-beta(k))*a_A)))+beta(k)*v_bB+sum(log(1+exp(beta(k)*Wh)));
            log_fk(k)=log_fk(k)-log(Z(k))-(1-beta(k))*log_ZA;
        end
        prob_beta_v=exp(log_fk-logsum(log_fk));
        beta_next=beta(randi([1,K]));
        c=c+prob_beta_v/N;
    end
    c_all(run,:)=c;
    Z=[Z(1),Z(2:K).*c(2:K)/c(1)];
    error(run)=max(abs(1/K-c));
end
figure;
plot(1:runtime,error);
title(['h',num2str(numhid),'.mat 模型收敛情况']);
xlabel('迭代次数');
ylabel('max(abs(1/K-c))');
logZZ=log(Z(end));
c_stable=c_all(runtime/2+1:runtime,:);
sigma_1=mean((c_stable(:,1)-mean(c_stable(:,1))).^2);
sigma_k=mean((c_stable(:,K)-mean(c_stable(:,K))).^2);
sigma_1k=mean((c_stable(:,1)-mean(c_stable(:,1))).*(c_stable(:,K)-mean(c_stable(:,K))));
var_log=sigma_1/(c_stable(runtime/2,1)^2)+sigma_k/(c_stable(runtime/2,K)^2)-2*sigma_1k/(c_stable(runtime/2,1))*(c_stable(runtime/2,K));
logZZ_up=logZZ+sqrt(var_log);
logZZ_down=logZZ-sqrt(var_log);
toc;
end