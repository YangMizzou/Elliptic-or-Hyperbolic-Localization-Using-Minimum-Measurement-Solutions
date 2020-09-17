% This program gives an example of how to call the routine BLUEest.m 
% using hyperbolic measurements to reproduce Fig. 15 or 16 in 
% Sanaa S. A. Al-Samahi, Yang Zhang, K. C. Ho, "Elliptic and hyperbolic 
% localizations using minimum measurement solutions", Elsevier Signal 
% Process., vol. 167, Feb. 2020.
%
% Sanaa S.A. Al-Samahi, Yang Zhang and K. C. Ho   02-28-2020
%
%       Copyright (C) 2020
%       Computational Intelligence Signal Processing Laboratory
%       University of Missouri
%       Columbia, MO 65211, USA
%       hod@missouri.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all; warning('off');

% ==================== settings ====================
% -- Choose between 2-D (Fig. 15) or 3-D (Fig. 16) --
% -- 2-D --
M=4;                               % number of receivers
dim=2;                             % Dimension
% % -- 3-D --
% M=5;                               % number of receivers
% dim=3;                             % Dimension

IndGeoResult=0+1;                    % show individual geometry result (0 or 1)

% NumEn=2000;                      % number of ensemble runs
NumEn=200;                         % number of ensemble runs

NoG=20;                            % number of RANDOM geometrIES
% ==================================================


% ===== generate localization geometries randomly =====
rand('state',0);

dd=1;
while (dd<=NoG)                   % generate sensor and source positions
    nn=0;
    ssg=100*(rand(dim,M)-0.5);
    for i=1:M-1,
        for j=i+1:M,
            if norm(ssg(:,i)-ssg(:,j))<20
                nn=nn+1;
            end
        end
    end
    if nn>0
        continue;
    else
        ss(:,:,dd)=ssg;
        uu(:,dd)=100*rand(dim,1)-50;
        dd=dd+1;
    end
end
% =====================================================

nsePwrAlldB=[-40:5:10];                       % measurement noise power in log-scale
CRLB=zeros(length(nsePwrAlldB),NoG);
mse_BLUEd=zeros(length(nsePwrAlldB),NoG);

for D=1:NoG,                                  % loop through geometries
    s=ss(:,:,D);
    u=uu(:,D);
    ro=sqrt(sum((repmat(u,1,M-1)-s(:,2:M)).^2))'-norm(u-s(:,1));   % true range measurement
    mse_BLUE=zeros(1,length(nsePwrAlldB));
    randn('state',319);                              % noise initialization
    noise=randn(M-1,NumEn);                          % random measurement noise
    
    for nsePwrIdx=1:length(nsePwrAlldB)
        fprintf('Geometry / 10log(NsePwr):  %d / %d\n',D,nsePwrAlldB(nsePwrIdx));
        nsePwr=10^(nsePwrAlldB(nsePwrIdx)/10);
        Q=(eye(M-1)+ones(M-1))/2;
        Q=Q*nsePwr;                                  % measurement noise covariance
        
        % ------------------------------ CRLB evaluation ----------------------------
        rho_us=repmat(u,1,M)-s; rho_us=rho_us./(ones(dim,1)*sqrt(sum(rho_us.^2)));
        G=[rho_us(:,2:end)-rho_us(:,1)*ones(1,M-1)]';
        crlb=inv(G'*inv(Q)*G);
        CRLB(nsePwrIdx,D)=trace(crlb);
        for ii=1:NumEn
            r=ro+sqrtm(Q)*noise(:,ii);
            
            % ---------------------------- Proposed Solution ----------------------------
            [u_hat]=BLUEest(s,r,Q,-1,1);
            mse_BLUE(nsePwrIdx)=mse_BLUE(nsePwrIdx)+sum(abs(u_hat-u).^2);
        end
    end
    mse_BLUEd(:,D)=mse_BLUE/NumEn;
    
    % -- show performance in individual geometry --
    if (IndGeoResult)
        figure;
        plot(nsePwrAlldB,10*log10(CRLB(:,D)),'-r');
        hold  on;
        plot(nsePwrAlldB,10*log10(mse_BLUEd(:,D)),'-bo');
        hold off;
        xlabel('10 log(\sigma^2)'); ylabel('10 log(MSE)');
        grid on;
        legend('CRLB','Proposed Solution','location','northwest');
        if (dim==2)
            axis([-40 10 -30 60]);
        else
            axis([-40 10 -30 100]);
        end;
    end;
    
end

mCRLB=mean(CRLB,2);
mBLUE=mean(mse_BLUEd,2);
figure;
plot(nsePwrAlldB,10*log10(mCRLB),'-r');
hold  on;
plot(nsePwrAlldB,10*log10(mBLUE),'-bo');
hold off;
xlabel('10 log(\sigma^2)'); ylabel('10 log(MSE)');
grid on;
legend('CRLB','Proposed Solution','location','northwest');
if (dim==2) 
    axis([-40 10 -30 60]);
else
    axis([-40 10 -30 100]);
end;


