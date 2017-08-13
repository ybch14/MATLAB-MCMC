close all;clear all;clc;
step=0.1;
n_x=0:step:10;
n_y=0:step:20;
% each column is a state also a point on plane
S=zeros(2,length(n_x)*length(n_y));
for i=1:length(n_x)
    S(1,(i-1)*length(n_y)+1:i*length(n_y))=n_x(i)*ones(1,length(n_y));
end
S(2,:)=repmat(n_y,1,length(n_x));
% the matrix is used to fix numerical error
s_find=round(10*S);
% Pai is stationary distribution,also target distribution
mu=[5;10];sigma=[1,1;1,4];
Pai=zeros(1,length(S));
for i=1:length(S)
    x=S(:,i)-mu;
    Pai(i)=1/(2*pi*sqrt(det(sigma)))*exp(-1/2*(x'/sigma*x))*step^2;
end
% T is referring matrix
% each row of T is uniform distribution U([0,0],[10,20])
% the probability is all the same, no need to calculate
% the number of sample times is 10000
total_number=1e5;
state=zeros(2,total_number);state(:,1)=[5;10];
for i=2:total_number
    Xn=state(:,i-1)';
    % generate Y uniformly
    Y(1)=round(rand()*10*10);
    Y(2)=round(rand()*20*10);
    Y_Loc=find(s_find(1,:)==Y(1)&s_find(2,:)==Y(2));
    Xn_Loc=find(s_find(1,:)==Xn(1)*10&s_find(2,:)==Xn(2)*10);
    % all the T(.) are the same so no need to multiple T
    r=min([1,Pai(Y_Loc)/(Pai(Xn_Loc))]);
    u=rand;
    if u<r
        state(:,i)=Y'/10;
    else
        state(:,i)=Xn';
    end
end
% plot the sample result compared with theorical distribution
data=mvnrnd(mu,sigma,total_number);
plot(data(:,1),data(:,2),'b.');
axis image;
axis([0,10,0,20]);
hold on;
plot(state(1,:),state(2,:),'r.');
axis image;
axis([0,10,0,20]);
figure;subplot(1,2,1);
hist3(state',[100,200]);
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
title('采样情况');
subplot(1,2,2);
hist3(data,[100,200]);
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
title('理论取值');

% calculate the correlation coefficient
EX=mean(state(1,:));
VarX=mean((state(1,:)-EX).^2);
EY=mean(state(2,:));
VarY=mean((state(2,:)-EY).^2);
Cov=mean((state(1,:)-EX).*(state(2,:)-EY));
p=Cov/sqrt(VarX*VarY);