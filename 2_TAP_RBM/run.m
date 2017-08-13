close all;clear all;clc;
% this is the main procedure for TAP algorithm
load data/test.mat

load data/h10.mat;
[logZZ_h10,logZZ_h10_up,logZZ_h10_down]=RBM_TAP(parameter_W,parameter_b,parameter_a);
loglik_h10=calculate_logprob(parameter_W,parameter_a,parameter_b,logZZ_h10,testbatchdata);
fprintf(1,'Estimated log-partition function of h10.mat: %f(%f,%f)\n',logZZ_h10,logZZ_h10_down,logZZ_h10_up);
fprintf(1,'Average estimated log_prob on the test data: %f\n', loglik_h10);

load data/h20.mat;
[logZZ_h20,logZZ_h20_up,logZZ_h20_down]=RBM_TAP(parameter_W,parameter_b,parameter_a);
loglik_h20=calculate_logprob(parameter_W,parameter_a,parameter_b,logZZ_h20,testbatchdata);
fprintf(1,'Estimated log-partition function of h20.mat: %f(%f,%f)\n',logZZ_h20,logZZ_h20_down,logZZ_h20_up);
fprintf(1,'Average estimated log_prob on the test data: %f\n', loglik_h20);

load data/h100.mat;
[logZZ_h100,logZZ_h100_up,logZZ_h100_down]=RBM_TAP(parameter_W,parameter_b,parameter_a);
loglik_h100=calculate_logprob(parameter_W,parameter_a,parameter_b,logZZ_h100,testbatchdata);
fprintf(1,'Estimated log-partition function of h100.mat: %f(%f,%f)\n',logZZ_h100,logZZ_h100_down,logZZ_h100_up);
fprintf(1,'Average estimated log_prob on the test data: %f\n', loglik_h100);

load data/h500.mat;
[logZZ_h500,logZZ_h500_up,logZZ_h500_down]=RBM_TAP(parameter_W,parameter_b,parameter_a);
loglik_h500=calculate_logprob(parameter_W,parameter_a,parameter_b,logZZ_h500,testbatchdata);
fprintf(1,'Estimated log-partition function of h500.mat: %f(%f,%f)\n',logZZ_h500,logZZ_h500_down,logZZ_h500_up);
fprintf(1,'Average estimated log_prob on the test data: %f\n', loglik_h500);

z=[logZZ_h10;logZZ_h20;logZZ_h100;logZZ_h500;...
    loglik_h10;loglik_h20;loglik_h100;loglik_h500];
save z.mat z;