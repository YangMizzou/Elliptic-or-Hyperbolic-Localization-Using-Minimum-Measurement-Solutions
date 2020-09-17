function [u_hat]=BLUEest(s,r,Q,LocType,SleSign)
% function [u_hat]=BLUEest(s,r,Q,LocType,SleSign)
% 
% This function realizes the BLUE for Elliptic and Hyperbolic Positioning 
% Using Minimum Measurement Solutions.
%
% Input parameter list:
% s: (Dim x M), sensor positions
%     s(:,1) is transmitter position for elliptic positioning
%     s(:,1) is reference sensor position for hyperbolic positioning
% r: ((M-1) x 1), input measurements
% Q: ((M-1) x (M-1)), measurement noise covariance
% LocType : either (+1) or (-1)
%           (+1) elliptic positioning
%           (-1) hyperbolic positioning
% SleSign: either 1 or 0
%            1  with individual selection method
%            0  without individual selection method
%     
% Output parameter list:
% u_hat:   object position estimate.
%
% The program can be used for 2-D(Dim=2) or 3-D(Dim=3) localization.
%
% Reference:
% Sanaa S. A. Al-Samahi, Yang Zhang, and K. C. Ho, "Elliptic and hyperbolic 
% localizations using minimum measurement solutions", Elsevier Signal Process., 
% vol. 167, Feb. 2020.
% 
% Yang Zhang, K. C. Ho and Sanaa S.A. Al-Samahi     02-28-2020
% 
%       Copyright (C) 2020
%       Computational Intelligence Signal Processing Laboratory
%       University of Missouri
%       Columbia, MO 65211, USA.
%       hod@missouri.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Dim,M]=size(s);                       % M = number of receivers
% Dim=dimension

s_ref=s(:,1);                          % refernece sensor
s_rem=s(:,2:end);                      % remaining sensors

SoR=1:M-1;
Th=6.635;                              % residual square error threshold
Csol=nchoosek(SoR,Dim);                % possible combinations in obtaining the minimum solution
[NoIS,~]=size(Csol);                   % number of individual solutions
sol=zeros(Dim,2*NoIS);
CorSolInd=ones(1,NoIS);
RptCnt = 1;                             % number of repetitions to recompute matrices B and G for the BLUE
for i=1:NoIS
    r_ind=r(Csol(i,:));
    s_rem_ind=s_rem(:,Csol(i,:));
    [sol(:,2*i-1),sol(:,2*i),CorSolInd(i)]=IndvLocSol(s_rem_ind,s_ref,r_ind,LocType);
end
[Sol_ind,Sol_ind_sub,~]=SolDetect(s_rem,s_ref,r,sol,Q,Csol,Th,LocType);

%=================================== BLUE ====================================
NoM=ceil((M-1)/Dim);
Combidx=Config_Comb(Dim,M-1);              % all possible combination choice for the BLUE
if SleSign==0
    CombChoice=Combidx(1,:);
elseif SleSign==1
    [CombChoice]=IndvSolSel(s_rem,s_ref,Sol_ind,CorSolInd,Q,Csol,Combidx,LocType);
end
NoCorSol=length(find(CorSolInd==1));     % number of correct individual solutions
if NoCorSol==0
    u_hat=1e10*ones(Dim,1);
else
    LocComb=Csol(CombChoice,:);            % locate combination of measurements
    Sol_CombAll=[Sol_ind(:,CombChoice);Sol_ind_sub(:,CombChoice)];
    u_hatVec=Inf(Dim,2^NoM);
    Res=Inf(1,2^NoM);
    
    for i=0:(2^NoM)-1
        bi=de2bi(i,NoM);
        bi=bi+ones(1,NoM);
        SolComb=[];
        for j=1:NoM
            SolComb=[SolComb,Sol_CombAll((bi(j)-1)*Dim+1:bi(j)*Dim,j)];
        end
        if max(SolComb,[],2)==Inf(Dim,1)
            continue;
        else
            u_hat=mean(SolComb,2);
            for k=1:RptCnt
                B=[];
                G=[];
                for j=1:NoM
                    s_ci=s_rem(:,LocComb(j,:));
                    G_cius_ci=repmat(u_hat,1,Dim)-s_ci; G_cius_ci=G_cius_ci./(ones(Dim,1)*sqrt(sum(G_cius_ci.^2)));
                    G_cius_ref=(u_hat-s_ref)/norm(u_hat-s_ref);
                    G_ci=[G_cius_ci+LocType*G_cius_ref*ones(1,Dim)]';
                    iG_ci=inv(G_ci);
                    B=blkdiag(B,iG_ci);
                    G=[G;eye(Dim)];
                end
                C=[];
                LocVec=reshape(LocComb',1,[]);
                for j=1:M-1,
                    C=[C;LocVec==j];
                end
                W=pinv(B*C'*Q*C*B');
                h=reshape(SolComb,[],1);
                u_hat=inv(G'*W*G)*G'*W*h;
            end
            u_hatVec(:,i+1)=u_hat;
            m_hat=sqrt(sum((repmat(u_hat,1,M-1)-s_rem).^2))'+LocType*norm(u_hat-s_ref);
            Res(i+1)=(r-m_hat)'*inv(Q)*(r-m_hat);
        end
    end
    [~,BLUEuIdx]=min(Res);
    u_hat=u_hatVec(:,BLUEuIdx);
end

end