% This program gives an example of how to call the routine BLUEest.m
% using elliptic measurements to reproduce Fig. 12 or 13 in
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
randn('state',0);                              % noise initialization

% ==================== settings ====================
% -- Choose between 2-D (Fig. 12) or 3-D (Fig. 13) --
% -- 2-D --
M=4;                               % number of (receivers + transmitter)
dim=2;                             % Dimension
% % -- 3-D --
% M=5;                               % number of (receivers + transmitter)
% dim=3;                             % Dimension

IndGeoResult=0;                    % show individual geometry result (0 or 1)

% NumEn=2000;                      % number of ensemble runs
NumEn=200;                         % number of ensemble runs

NoG=20;                            % number of RANDOM geometrIES
% ==================================================


N=M-1;                             % number of measurements

%%%%%[uu,ts]=RdmGeoGenEllip(dim,N,NoG);             % generate random geometries

% ===== generate localization geometries randomly =====
rand('state',0);
ts=zeros(dim,M,NoG); uu=zeros(dim,NoG);

% -- 2-D --
if dim==2
    dd=1;
    while (dd<NoG+1)
        nn=0;
        ssg=200*(rand(dim,M)-0.5);
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
            ts(:,:,dd)=ssg;
            dd=dd+1;
        end
    end
    
    nn=1;
    for dd=1:NoG
        while (nn~=0)
            nn=0;
            uu(:,dd)=200*(rand(dim,1)-0.5);
            if prod(uu(:,dd)>-50 & uu(:,dd)<50)
                nn=1;
                continue;
            end
        end
        nn=1;
    end
    
    % -- 3-D --
elseif dim==3
    nn=1;
    for dd=1:NoG
        while (nn~=0)
            nn=0;
            uu(:,dd)=200*(rand(dim,1)-0.5);
            if prod(uu(:,dd)>-50 & uu(:,dd)<50)
                nn=1;
                continue;
            end
        end
        nn=1;
    end
    
    dd=1;
    while (dd<NoG+1)
        nn=0;
        ssg=200*(rand(dim,M)-0.5);
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
            ts(:,:,dd)=ssg;
            dd=dd+1;
        end
    end
    
end
% =====================================================

nsePwrAlldB=[-40:5:20];                        % measurement noise power in log-scale
CRLB=zeros(length(nsePwrAlldB),NoG);
mse_BLUEd=zeros(length(nsePwrAlldB),NoG);

for D=1:NoG
    t=ts(:,1,D);
    s=ts(:,2:end,D);
    u=uu(:,D);
    ro=sqrt(sum((repmat(u,1,N)-s).^2))'+norm(u-t);     % true range measurement
    mse_BLUE=zeros(1,length(nsePwrAlldB));
    
    for nsePwrIdx=1:length(nsePwrAlldB)
        fprintf('Geometry / 10log(NsePwr):  %d / %d\n',D,nsePwrAlldB(nsePwrIdx));
        nsePwr=10^(nsePwrAlldB(nsePwrIdx)/10);
        Q=eye(N)*nsePwr;                              % measurement noise covariance
        
        % ------------------------------ CRLB evaluation ----------------------------
        rho_us=repmat(u,1,N)-s; rho_us=rho_us./(ones(dim,1)*sqrt(sum(rho_us.^2)));
        rho_uso=(u-t)/norm(u-t);
        G=[rho_us+rho_uso*ones(1,N)]';
        crlb=inv(G'*inv(Q)*G);
        CRLB(nsePwrIdx,D)=trace(crlb);
        for ii=1:NumEn
            noise=randn(N,1);
            r=ro+sqrt(nsePwr)*noise;
            % ---------------------------- Proposed Solution ----------------------------
            [u_hat]=BLUEest([t,s],r,Q,1,1);
            mse_BLUE(nsePwrIdx)=mse_BLUE(nsePwrIdx)+sum(abs(u_hat-u).^2);
        end
    end
    mse_BLUEd(:,D)=mse_BLUE/NumEn;
    
    % -- show performance in individual geometry --
    if (IndGeoResult)
        figure;
        plot(nsePwrAlldB,10*log10(CRLB(:,D)),'-r');
        hold on;
        plot(nsePwrAlldB,10*log10(mse_BLUEd(:,D)),'-bo');
        hold off;
        xlabel('10 log(\sigma^2)'); ylabel('10 log(MSE)');
        grid on;
        legend('CRLB','Proposed Solution','location','northwest');
        if (dim==2)
            axis([-40 20 -30 45]);
        else
            axis([-40 20 -30 60]);
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
    axis([-40 20 -30 45]);
else
    axis([-40 20 -30 60]);
end;



