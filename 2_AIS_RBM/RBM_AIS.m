function [logZZ,logZZ_up,logZZ_down]=RBM_AIS(parameter_W,parameter_a,parameter_b,numruns,beta,batchdata)
% parameter_W is the matrix W_ij [numvis, numhid]
% parameter_a is the row vector of hidbiases a_j [1, numhid]
% parameter_b is the row vector of visbiases b_j [1, numvis]
% numruns is number of AIS run times
% beta is a row vector contain beta's [1, num_beta]
tic;
[numcases,numdims,numbatches]=size(batchdata);
 count_int = zeros(numdims,1);
 for batch=1:numbatches
    xx = sum(batchdata(:,:,batch));
    count_int = count_int + xx';
 end
 lp=5;
 p_int = (count_int+lp*numbatches)/(numcases*numbatches+lp*numbatches);
 log_base_rate = log( p_int) - log(1-p_int);

close all;
[numvis,numhid]=size(parameter_W);
bes_A=log_base_rate';
b_A=repmat(bes_A,numruns,1);
a_B=repmat(parameter_a,numruns,1);
b_B=repmat(parameter_b,numruns,1);

logw=zeros(numruns,1);
v=repmat(1./(1+exp(-bes_A)),numruns,1);
v=v>rand(numruns,numvis);
logw=logw-(v*bes_A'+numhid*log(2));

Wh=v*parameter_W+a_B;
Bv_base=v*bes_A';
Bv=v*parameter_b';

for bb=beta(2:end-1)
    %     fprintf('beta=%f\n',bb);
    expWh=exp(bb*Wh);
    logw=logw+(1-bb)*Bv_base+bb*Bv+sum(log(1+expWh),2);
    prob_hB_v=expWh./(1+expWh);
    h_B=prob_hB_v>rand(numruns,numhid);
    v=1./(1+exp(-(1-bb)*b_A-bb*(h_B*parameter_W'+b_B)));
    v=v>rand(numruns,numvis);
    Wh=v*parameter_W+a_B;
    Bv_base=v*bes_A';
    Bv=v*parameter_b';
    expWh=exp(bb*Wh);
    logw=logw-((1-bb)*Bv_base+bb*Bv+sum(log(1+expWh),2));
end

expWh=exp(Wh);
logw=logw+v*parameter_b'+sum(log(1+expWh),2);

r_AIS=log(sum(exp(logw)))-log(numruns);
aaa=mean(logw(:));
logstd=log(std(exp(logw-aaa)))+aaa-log(numruns)/2;

logZZ_base=sum(log(1+exp(bes_A)))+(numhid)*log(2);
logZZ=r_AIS+logZZ_base;

logZZ_up=logsum([logstd;r_AIS])+logZZ_base;
logZZ_down=logdiff([logstd;r_AIS])+logZZ_base;
toc;
end