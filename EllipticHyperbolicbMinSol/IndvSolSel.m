function [CombChoice]=IndvSolSel(s_rem,s_ref,SolInd,CorSolInd,Q,Csol,CombIdx,LocType)
% function [CombChoice]=IndvSolSel(s_rem,s_ref,SolInd,CorSolInd,Q,Csol,CombIdx,LocType)
%
% This function realizes the minimum solutions selection.
%
% Input parameter list:
% s_rem:  (Dim x N), receiver position matrix, N is the number of receivers.
% s_ref:  (Dim x 1), transmitter position for elliptic positioning or
%                    reference sensor position for hyperbolic positioning or
% SolInd:     individual solutions
% CorSolInd:  correct solution indicator
% Q:          (N x N), Covariance matrix of the range measurements
% Csol:       all possible combinations in obtaining the minimum measurement solution
% CombIdx:    possible combination choices of individual solutions for BLUE
% LocType :   either (+1) or (-1)
%               (+1) elliptic positioning; 
%               (-1) hyperbolic positioning
%     
% Output parameter list:
% CombChoice:  returned combination choice
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
%       Columbia, MO 65211, USA
%       hod@missouri.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Dim,~]=size(s_rem);
[NoIS,~]=size(Csol);                 % NoIS = number of individual solutions
for i=1:NoIS                         % compute determinant of individual covariance
    if CorSolInd(i)==0
        RIndr(i)=Inf;
    else
        s_Ind=s_rem(:,Csol(i,:));
        GInd_us_rem=repmat(SolInd(:,i),1,Dim)-s_Ind; GInd_us_rem=GInd_us_rem./(ones(Dim,1)*sqrt(sum(GInd_us_rem.^2)));
        GInd_us_ref=(SolInd(:,i)-s_ref)/norm(SolInd(:,i)-s_ref);
        GInd=(GInd_us_rem+LocType*GInd_us_ref*ones(1,Dim))';
        QInd=Q(Csol(i,:),Csol(i,:));
        CovInd=inv(GInd)*QInd*inv(GInd)';
        [~,DInd]=eig(CovInd);
        DInd=diag(DInd);
        if DInd(1)<0||DInd(2)<0
            RIndr(i)=Inf;
        else
            RIndr(i)=prod(DInd);
        end
    end
end

TotalComb=size(CombIdx,1);
Qc=Inf(1,TotalComb);
for i=1:TotalComb                      % evaluate the quality factor
    Qc(i)=prod(RIndr(CombIdx(i,:)));
end
[~,MinQc]=min(Qc);
CombChoice=CombIdx(MinQc,:);

end

