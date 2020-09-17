function [u_hat1,u_hat2,CorSolInd]=IndvLocSol(s_rem,s_ref,r,Loctype)
% function [u_hat1,u_hat2,CorSolInd]=IndvLocSol(s_rem,s_ref,r,Loctype)
%
% This function realizes the individual solution using minimum number of 
% measurements. 
%
% Input parameter list:
% s_rem:  (Dim x M), receiver position matrix, N is the number of receivers.
% s_ref:  (Dim x 1), transmitter position for elliptic positioning or
%                    reference sensor position for hyperbolic positioning or
%                    
% r:      (M x 1), noisy range measurements.
% Loctype:      either (+1) or (-1)
%               (+1) elliptic positioning; 
%               (-1) hyperbolic positioning
%     
% Output parameter list:
% u_hat1:       quadratic solution 1
% u_hat2:       quadratic solution 2
% CorSolInd:    CorSolInd indicates if the two solutions give positive 
%               right side of eqn.(4).
%               1 at least one solution is correct; 
%               0 neither solution is correct
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Dim,N]=size(s_ref);           % N = number of receivers
                               % Dim = dimension
alpha=zeros(Dim);
k=zeros(1,Dim);
CorSolInd=1;
for i=1:Dim
    alpha(:,i)=Loctype*(s_rem(:,i)-s_ref)/r(i);
    k(i)=Loctype*(r(i)^2+norm(s_ref)^2-norm(s_rem(:,i))^2)/(2*r(i));
end

A=[alpha(:,1)*alpha(:,1)'-eye(Dim),k(1)*alpha(:,1)+s_ref;
    k(1)*alpha(:,1)'+s_ref', k(1)^2-norm(s_ref)^2];
if Dim==2
    B1=[alpha(2,2)-alpha(2,1),k(2)-k(1)]/(alpha(1,1)-alpha(1,2));
elseif Dim==3
    B1=inv([alpha(1:2,1)'-alpha(1:2,2)';alpha(1:2,1)'-alpha(1:2,3)'])...
        *([alpha(3,2)-alpha(3,1), k(2)-k(1);alpha(3,3)-alpha(3,1), k(3)-k(1)]);
end
B=[B1;eye(2)];
H=B'*A*B;

Delta=H(1,2)^2-H(1,1)*H(2,2);   
if Delta<0
    Delta=0;                  % approximate complex solutins with real part 
end

uLE=(-H(1,2)+sqrt(Delta))/H(1,1);
ur=B*([uLE;1]);
u_hat1=ur(1:Dim);

uLE=(-H(1,2)-sqrt(Delta))/H(1,1);
ur=B*([uLE;1]);
u_hat2=ur(1:Dim);

uPost1=k(1)+alpha(:,1)'*u_hat1;
uPost2=k(1)+alpha(:,1)'*u_hat2;

if uPost1<0 && uPost2<0
    CorSolInd=0;
elseif uPost1<0 && uPost2>=0
    u_hat1=Inf;
elseif uPost1>=0 && uPost2<0
    u_hat2=Inf;
end

end

