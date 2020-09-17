function [sol,sol_sub,ResErr]=SolDetect(s_rem,s_ref,r,solAll,Q,Csol,Th,LocType)
% function [sol,sol_sub,ResErr]=SolDetect(s_rem,s_ref,r,solAll,Q,Csol,Th,LocType)
%
% This function realizes the elimination of the minimum solution ambiguity. 
%
% Input parameter list:
% s_rem:  (Dim x N), receiver position matrix, N is the number of receivers.
% s_ref:  (Dim x 1), transmitter position for elliptic positioning or
%                    reference sensor position for hyperbolic positioning or
% r:      (N x 1), noisy measurements.
% solAll: minimum measurement solutions
% Q:      (N x N), Covairance matrix of the range measurements.
% Csol:   all possible combinations in obtaining minimum measurement solution
% Th:     threshold
% LocType :   either (+1) or (-1)
%               (+1) elliptic positioning; 
%               (-1) hyperbolic positioning
%     
% Output parameter list:
% sol:     proper solutions.
% sol_sub: the other set of solutions whose residual square errors are 
%          relatively larger but smaller than the threshold
% ResErr:  residual square errors
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

[Dim,N]=size(s_rem);            % N = number of receivers
                                % Dim = dimension
[~,n]=size(solAll);             % n = number of individual solutions
sol=zeros(Dim,n/2);
sol_sub=Inf(Dim,n/2);
ResErr=zeros(2,n/2);
for i=1:(n/2),
    %=========== compute corresponding residual square error ============
    for j=1:2
        r_i=r;
        r_i(Csol(i,:))=[];
        Q_i=Q;
        Q_i(:,Csol(i,:))=[];
        Q_i(Csol(i,:),:)=[];
        r_hat=sqrt(sum((repmat(solAll(:,2*(i-1)+j),1,N)-s_rem).^2))'+LocType*norm(solAll(:,2*(i-1)+j)-s_ref);
        r_hati=r_hat;
        r_hati(Csol(i,:))=[];
        ResErr(j,i)=(r_i-r_hati)'*inv(Q_i)*(r_i-r_hati);
        if isnan(ResErr(j,i))
            ResErr(j,i)=Inf;
        end
    end
    
    %=========== compare residual square error, choose proper solution ============
    if (ResErr(1,i)<ResErr(2,i))
        sol(:,i)=solAll(:,2*i-1);
        if ResErr(2,i)<=Th
            sol_sub(:,i)=solAll(:,2*i);
        end
    else
        sol(:,i)=solAll(:,2*i);
        if ResErr(1,i)<=Th
            sol_sub(:,i)=solAll(:,2*i-1);
        end
    end
end

end

