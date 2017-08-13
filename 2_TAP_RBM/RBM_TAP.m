function [logZZ,logZZ_up,logZZ_down]=RBM_TAP(parameter_coupling,parameter_vis,parameter_hid)
tic;
W=parameter_coupling;[numvis,numhid]=size(W);
a=parameter_vis;b=parameter_hid;
mh=0.1*ones(1,numhid);mv=0.1*ones(1,numvis);
N=10000;
mh_all=zeros(N,numhid);mv_all=zeros(N,numvis);
for i=1:N
    mh=sigm(b+mv*W-(mh-1/2).*((mv-mv.^2)*(W.^2)));
    mv=sigm(a+mh*W'-(mv-1/2).*((mh-mh.^2)*(W.^2)'));
    mh_all(i,:)=mh;
    mv_all(i,:)=mv;
end
mh_stable=mh_all(N-N/10+1:N,:);
mv_stable=mv_all(N-N/10+1:N,:);
logZZ_stable=zeros(N/10,1);
for i=1:(N/10)
    mht=mh_stable(i,:);mvt=mv_stable(i,:);
    st=mht*log(mht)'+(1-mht)*log(1-mht)'+(mvt*log(mvt)'+(1-mvt)*log(1-mvt)');
    logZZ_stable(i)=-(st-a*mvt'-b*mht'-mvt*W*mht'-(mvt-mvt.^2)*((W.^2)/2)*(mht-mht.^2)');
end
logZZ=mean(logZZ_stable);
logZZ_std=std(logZZ_stable);
logZZ_up=logZZ+logZZ_std;
logZZ_down=logZZ-logZZ_std;
toc;
end