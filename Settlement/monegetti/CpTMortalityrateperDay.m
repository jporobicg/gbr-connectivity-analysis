%% alive and competent curve
clear all

nolarvae=1e5;
% [palive_comp]= Mortality_Competency_COTS0123(Npart);
% % Second, acquistion, loss, and survival parameters for the piecewise model
% Second, acquistion, loss, and survival parameters for the piecewise model
awe = 1.290078 ;    % acquisition rate
b1we = 0.00188621;  % loss rate parameter pre-change point
b2we = 0.3969605 ;  % loss rate parameter post-change point
v1we = 0.3652014 ;  % loss shape parameter for pre-change point
tcwe = 3.331252;    % length of obligate pre-competent period
Tcpwe = 69.91391;   % change point for loss of competence
u1ww = 0.4002655 ;  % mortality rate parameter pre-change point
u2ww = 0.01923417 ; % mortality rate parameter post-change point
v1ww = 2.892410 ;   % mortality shape parameter pre-change point
v2ww = 1.715599 ;   % mortality shape parameter post-change point
Tcpww = 2.583284;   % change point for mortality

%% lifetime 
time=0:1/24:45;
for t=1:length(time)
if time(t)<=Tcpww
    
  alive(t) = exp(-(u1ww*time(t))^v1ww);
else
    alive(t)= exp(-(u1ww*Tcpww)^v1ww)*exp(-(u2ww*(time(t)-Tcpww))^v2ww);
end
end

plifetime1 = wblrnd(1/u1ww, v1ww, nolarvae,1);
plifetime2 = wblrnd(1/u2ww, v2ww, nolarvae,1) + Tcpww;
totlifetime=zeros(nolarvae,1)*nan;
ll1=find(plifetime1 < Tcpww);ll2=find(plifetime1 >= Tcpww);
totlifetime(ll1)=plifetime1(ll1);totlifetime(ll2)=plifetime2(ll2);
clear ll1 ll2

time=1/24:1/24:45;

for i=1:length(time)
    ll=find(palive_comp(:,1)<=time(i) & palive_comp(:,2)>=time(i));
    AliveComp(i)=length(ll)/Npart;
end

plot(time,AliveComp*100)

% https://stackoverflow.com/questions/12358860/exponential-random-numbers-with-a-bound-in-matlab
% sizeOut = [1, 1000]; % sample size
% mu = 100; % parameter of exponential 
% r1 = 50;  % lower bound
% r2 = 150; % upper bound
% 
% r = exprndBounded(mu, sizeOut, r1, r2); % bounded output    
% 
% function r = exprndBounded(mu, sizeOut, r1, r2);
% 
% minE = exp(-r1/mu); 
% maxE = exp(-r2/mu);
% 
% randBounded = minE + (maxE-minE).*rand(sizeOut);
% r = -mu .* log(randBounded);