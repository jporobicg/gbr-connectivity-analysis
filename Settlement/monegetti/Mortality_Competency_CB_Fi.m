%% Specify parameters %%
function [palive_comp,npalive_comp]= Mortality_Competency_CB_Fi(Nlarvae)
% First competence acquistion, loss, and survival parameters for the non-piecwise model
% parameters for the *non*-piecewise model (a la Connolly & Baird)
ae = 0.6587968;     % acquisition rate
be = 0.01431578 ;   % loss rate
tce = 2.907979 ;    % length of obligate pre-competent period
uw = 0.09861929 ;   % mortality rate parameter
vw = 0.5906351 ;    % mortality shape parameter

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

%% Specify the number of larvae released %%

nolarvae = Nlarvae;

%% Calculate when larvae are alive and competent to settle %%

%% %%%%%%%% Non-piecewise model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This model is the best fitting model given the approach of
% Connolly & Baird 2010 and the Figueiredo et al papers

npacquisition = exprnd( 1/ae,nolarvae,1) + tce;
% Following the pre-competent period, larvae acquire competence at a
% rate a.e. The time that larvae acquire competence can be drawn from
% an exponential distribution. Therefore, the time that a larva is
% competent to settle is this randomly drawn time plus the pre-competent
% period, tc.e.
% np.acquisition is a vector of times of acquisition of settlement
% competence for the non-piecewise model.The ith element is the time
% that larva i acquires competence.

nploss = exprnd(1/be,nolarvae,1) + npacquisition;
% Once competence to settle has been acquired, larvae lose settlement
% competence at a rate b.e. The time that larvae lose competence can be
% drawn from an exponential distribution. Therefore, the time that a
% larva loses its settlement competence is this randomly drawn time plus
% the time of acquisition of settlement competence.
% np.loss is a vector of the times of loss of settlement competence for
% the non-piecewise model. The ith element is the time that larva i
% loses competence.

nplifetime = wblrnd(1/uw, vw,nolarvae,1 );
% Larvae die at temporally varying rate which increases with time
% according to the rate parameter u.e (note rate = 1/scale parameter) and
% shape parameter v.e. The time that larvae
% die can be drawn from a weibull distribution.
% np.lifetime is a vector of lifetimes for the non-piecewise model. The
% ith element is the lifetime for larva i.

% When are the larvae alive and competent to settle?
% npalive_comp = matrix(length(nolarvae),  2)
npalive_comp=zeros(nolarvae,2)*nan;
ll1=find(npacquisition<nplifetime);ll2=find(npacquisition>=nplifetime);
npalive_comp(ll1,1)=npacquisition(ll1);
npalive_comp(ll2,1)=0;
clear ll1 ll2

ll1=find(npacquisition>nplifetime);
npalive_comp(ll1,2)=0;
ll2=find(npacquisition<=nplifetime);% && nploss < nplifetime);
% npalive_comp2(ll2,2)=nploss(ll2)
% ll3=find(npacquisition>=nplifetime && nploss > nplifetime);
% npalive_comp2(ll3,2)=nploss(ll3)


for i=1:length(ll2)
    
%     if npacquisition(i) < nplifetime(i)
%         npalive_comp(i,1)= npacquisition(i);
%     else
%         npalive_comp(i,1)=0;
%     end
    % npalive_comp(i,1) = ifelse(np.acquisition[i] < np.lifetime[i], np.acquisition[i], 0)
    % If competence acquisition occurs before the larva dies, record the competence
    % acquisition time as np.acquisition[i]. Otherwise, we will put zero here
    % AND in the competence loss time (i.e., competence duration is zero)
    
%     if npacquisition(ll2(i)) > nplifetime(i)
%         npalive_comp(i,2) = 0;
%     else
        if nploss(ll2(i)) < nplifetime(ll2(i))
            npalive_comp(ll2(i),2) = nploss(ll2(i));
        else
            npalive_comp(ll2(i),2) = nplifetime(ll2(i));
        end
%     end
    %   npalive_comp[i,2] = ifelse(np.acquisition[i] > np.lifetime[i], 0,
    % If acquisition time is after larva has died, record competence loss time as zero;
    % as noted above, this ensures that there is no competence period
    %                   ifelse(np.loss[i] < np.lifetime[i], np.loss[i],
    %                          np.lifetime[i]))
    %   % Otherwise, if the larva loses competence before it dies, record competence
    % loss time as np.loss; whereas if it dies before competence is lost, record
    % the end of its competence time as its time of death
    
    
end
% np.alive_comp is a matrix defining the time period for which larvae
% are alive AND competent to settle for the non-piecewise model. The
% first column represents when larvae are alive and first become
% competent to settle, and the second column represents when larvae die
% or lose competence (whichever occurs first).
% Therefore, time np.alive_comp[i,1] to time np.alive_comp[i,2]
% represents the time period for which larva i is alive and competent
% to settle.

%% Piecewise model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Piecewise model

pacquisition = exprnd(1/awe,nolarvae,1) + tcwe;
% Following the pre-competent period, larvae acquire competence at a
% rate a.we. The time that larvae acquire competence can be drawn from
% an exponential distribution. Therefore, the time that a larva is
% competent to settle is this randomly drawn time plus the pre-competent
% period, tc.we.
% p.acquisition is a vector of times of acquisition of settlement
% competence for the piecewise model.The ith element is the time
% that larva i acquires competence.
% NOTE: This is the same as in the other model, just with a different
% acquisition rate parameter


% To simulate the competence loss times, we first generate the time competence
% is lost if the pre-switch time competence duration was the true competence
% duration (that's p.loss1). For any larvae for which p.loss1 < Tcp, we
% will record that larva's competence loss time as being whatever is given
% in p.loss1. For any larvae for which p.loss1 > Tcp, their competence loss
% time becomes determined by the second vector, p.loss2. We implement
% this by first simulating *both* competence loss times for each larvae,
% and THEN checking to see which competence loss time applies to which larvae.
% This check is implemented in the for-loop below.

ploss1 = wblrnd( 1/b1we, v1we,nolarvae,1) + pacquisition;
% Once competence to settle has been acquired, larvae lose settlement
% competence at a temporally varying rate. Prior to Tcp.we, larvae
% loss competence at a rate according to the shape parameters v1.we
% and RATE parameter b1.we. The time that larvae lose competence can be
% drawn from a Weibull distribution. Therefore, the time that a
% larva loses its settlement competence is this randomly drawn time plus
% the time of acquisition of settlement competence.
% p.loss1 is a vector of the times of loss of settlement competence for
% the piecewise model according to the rate of loss of competence prior
% to Tcp.we. The ith element is the time that larva i
% loses competence.

ploss2 = exprnd(1/b2we,nolarvae,1) + Tcpwe;
% Following time Tcp.we, larvae lose competence at a constant rate b2.we
% The time that larvae lose competence can be drawn from an exponential
% distribution. Therefore, assuming that a larva remains competent
% past time Tcp.we, the time that a larva loses its settlement
% competence is this randomly drawn time plus the time of the
% change-point, Tcp.we.
% p.loss2 is a vector of the times of loss of settlement competence for
% the piecewise model, assuming that larvae are remain competent past
% time Tcp.we. The ith element is the time that larva i loses competence.

% To simulate lifetimes with the piecewise model, we proceed analogously
% to the piecewise competence loss specified above

plifetime1 = wblrnd(1/u1ww, v1ww, nolarvae,1);
% Larvae die at temporally varying rate. Prior to Tcp.ww larvae die at
% a rate according to the parameters u.ww and v.ww. The time that
% larvae die can be drawn from a weibull distribution.
% p.lifetime1 is a vector of lifetimes for the piecewise model
% according to the rate of mortality prior to Tcp.ww. The ith element
% is the lifetime for larvae i.

plifetime2 = wblrnd(1/u2ww, v2ww, nolarvae,1) + Tcpww;
% Following time Tcp.ww, larvae die at a different rate according to
% the parameters u2.ww and v2.ww. The time that larvae die can be
% drawn from a Weibull distribution. The lifetime for larvae that
% survive at least to time Tcp.ww is this randomly drawn lifetime plus
% time Tcp.ww.
% p.lifetime2 is a vector of lifetimes for the piecewise model,
% assuming larvae die post Tcp.ww. The ith element is the lifetime for
% larvae i.

% When are the larvae alive and competent to settle?

% palive_comp = matrix((nolarvae),  2);

% totlifetime = c()
% totcomp = c()
totlifetime=zeros(nolarvae,1)*nan;
ll1=find(plifetime1 < Tcpww);ll2=find(plifetime1 >= Tcpww);
totlifetime(ll1)=plifetime1(ll1);totlifetime(ll2)=plifetime2(ll2);
clear ll1 ll2
% totlifetime[i] = ifelse(plifetime1[i] < Tcpww, plifetime1[i], plifetime2[i])
    % tot.lifetime[i] is the complete lifetime for larva i.
    % If a larva dies before the changepoint Tcp.ww its lifetime is
    % p.liftime1. If a larva survives to the time of the changepoint,
    % its lifetime is p.lifetime2.
totcomp=zeros(nolarvae,1)*nan;
ll1=find(ploss1<Tcpwe);ll2=find(ploss1>=Tcpwe);
totcomp(ll1)=ploss1(ll1);totcomp(ll2)=ploss2(ll2);
clear ll1 ll2
    %totcomp[i] = ifelse(p.loss1[i] < Tcpwe, ploss1[i], ploss2[i])
    % tot.comp[i] is the the time that larva i loses competence to settle.
    % If a larva loses competence befor the changepoint Tcp.we the time
    % that it lost competence is from p.loss1. If a larva lost competence
    % after the changepoint, the time it lost competence is from p.loss2.
 palive_comp=zeros(nolarvae,2)*nan;   
ll1=find(pacquisition<totlifetime);ll2=find(pacquisition>=totlifetime);
palive_comp(ll1,1)=pacquisition(ll1);palive_comp(ll2,1)=0;
clear ll1 ll2
% p.alive_comp[i,1] is the time that larva i is alive and becomes
    % competent to settle. If the larva dies before it would have become
    % competent, time 0 is noted.
ll1=find(pacquisition>totlifetime); ll2=find(pacquisition<=totlifetime);  
palive_comp(ll1,2)=0;
for i=1:length(ll2)
%     if plifetime1(i) < Tcpww
%         totlifetime(i)=plifetime1(i);
%     else
%         totlifetime(i)=plifetime2(i);
%     end
%     if ploss1(i) < Tcpwe
%         totcomp(i)=ploss1(i);
%     else
%         totcomp(i)=ploss2(i);
%     end
%     if pacquisition(i) < totlifetime(i)
%         palive_comp(i,1)=pacquisition(i);
%     else
%         palive_comp(i,1)=0;
%     end
    %   palive_comp[i,1] = ifelse(pacquisition[i] < totlifetime[i],
    %                              pacquisition[i], 0)
    
    
%     if pacquisition(i) > totlifetime(i)
%         palive_comp(i,2)=0;
%     else
        if totcomp(ll2(i)) < totlifetime(ll2(i))
            palive_comp(ll2(i),2)=totcomp(ll2(i));
        else
            palive_comp(ll2(i),2)=totlifetime(ll2(i));
        end
%     end
    
    %   palive_comp[i,2] = ifelse(pacquisition[i] > totlifetime[i], 0,
    %                   ifelse(totcomp[i] < totlifetime[i], totcomp[i],
    %                          totlifetime[i]))
    % p.alive_comp[i,2] is the time that larva i dies or loses competence.
    % If the larva dies before it would have become competent a 0 is
    % noted. If the larva loses competence before it dies then tot.comp
    % is noted. If the larva dies before it would have lost competence
    % then tot.lifetime is noted.
end  % p.alive_comp is a matrix defining the time period for which larvae
% are alive AND competent to settle for the piecewise model. The first
% column represents when larvae are alive and first become competent to
% settle, and the second column represents when larvae die or lose
% competence (whichever occurs first).
% Therefore, time np.alive_comp[i,1] to time np.alive_comp[i,2]
% represents the time period for which larva i is alive and competent
% to settle.

