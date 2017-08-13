close all;clear all;clc;
% this is the main procedure for SAMS algorithm
load data/test.mat;
load data/h10.mat;
W_h10=parameter_W;a_h10=parameter_a;b_h10=parameter_b;
load data/h20.mat;
W_h20=parameter_W;a_h20=parameter_a;b_h20=parameter_b;
load data/h100.mat;
W_h100=parameter_W;a_h100=parameter_a;b_h100=parameter_b;
load data/h500.mat;
W_h500=parameter_W;a_h500=parameter_a;b_h500=parameter_b;
[Zeta,std_h20,std_h100,std_h500]=RBM_SAMS(W_h10 ,b_h10 ,a_h10 ,...  
    W_h20 ,b_h20 ,a_h20 ,...
    W_h100,b_h100,a_h100,...
    W_h500,b_h500,a_h500);
logZZ_h10=calculate_true_partition(W_h10,a_h10,b_h10);
logZZ=logZZ_h10+Zeta;

loglik_h10=calculate_logprob(W_h10,a_h10,b_h10,logZZ(1),testbatchdata);
fprintf(1,'Estimated log-partition function of h10.mat: %f(%f,%f)\n',logZZ(1),logZZ(1),logZZ(1));
fprintf(1,'Average estimated log_prob on the test data: %f\n', loglik_h10);

loglik_h20=calculate_logprob(W_h20,a_h20,b_h20,logZZ(2),testbatchdata);
fprintf(1,'Estimated log-partition function of h10.mat: %f(%f,%f)\n',logZZ(2),logZZ(2)-std_h20,logZZ(2)+std_h20);
fprintf(1,'Average estimated log_prob on the test data: %f\n', loglik_h20);

loglik_h100=calculate_logprob(W_h100,a_h100,b_h100,logZZ(3),testbatchdata);
fprintf(1,'Estimated log-partition function of h10.mat: %f(%f,%f)\n',logZZ(3),logZZ(3)-std_h100,logZZ(3)+std_h100);
fprintf(1,'Average estimated log_prob on the test data: %f\n', loglik_h100);

loglik_h500=calculate_logprob(W_h500,a_h500,b_h500,logZZ(4),testbatchdata);
fprintf(1,'Estimated log-partition function of h10.mat: %f(%f,%f)\n',logZZ(4),logZZ(4)-std_h500,logZZ(4)+std_h500);
fprintf(1,'Average estimated log_prob on the test data: %f\n', loglik_h500);

z=[logZZ';...
    loglik_h10;loglik_h20;loglik_h100;loglik_h500];
save z.mat z;