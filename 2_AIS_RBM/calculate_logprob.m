function [logprob] = calculate_logprob(vishid,hidbiases,visbiases,logZ,batchdata)
[~,~,numbatches]=size(batchdata);

 data = [];
 for ii=1:numbatches
   data = [data; batchdata(:,:,ii)];
 end
 numcases = size(data,1);

 pd = data*visbiases' + sum(log(1+exp(ones(numcases,1)*hidbiases + data*vishid)),2);
 logprob = sum(pd)/numcases  - logZ;



