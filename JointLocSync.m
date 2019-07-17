close all
clear all;
% Near-Field Joint Localization and Synchronization
% (c) Henk Wymeersch, 2019
% Usage: this code generates position estimates



% System parameters
% -----------------
K = 1024;                   % number of subcarriers - 1
N = 128;                    % number of antennas -1
fc = 28;                    % carrier [GHz]
c = 0.3;                    % speed of light [m/ns]
lambda = c/fc;              % carrier wavelength [m]                
Delta = lambda/2;           % antenna spacing in [m] -> setting too large creates sidelobes
W = 0.1;                    % bandwidth [GHz]
Pt=1;                       % mW
N0 = 290*1.381e-23*1e3*1e9; % noise PSD in mW/GHz (290 Kelvin * Boltzmann constant in W/Hz)
bias=20;                     % UE clock bias in m
useSubArray=1;
P = Pt/W*ones(1,K+1);          % energy per subcarrier   
Deltaf=W/(K+1);                % subcarrier spacing
rf=Deltaf/fc;                  % ratio 



steps=15;                   % number of distances we will consider
sims=20;                     % number of realization per distance, set to larger value for smoother curve
xvec=logspace(-1,1.5,steps);
yvec=logspace(-1,1.5,steps);
xaxisvec=sqrt(xvec.^2+yvec.^2);
disp('start simulation')
for l=1:steps
    disp(['step ' num2str(l) ' of '  num2str(steps) ' completed.']);
    xUE=[xvec(l) yvec(l)];
    % Generate the received signal
    % ----------------------------    
    d=norm(xUE);                   % distance to UE
    x=xUE(1);                      % x coordinate
    y=xUE(2);                      % y coordinate
    theta=acos(x/d);               % AOA    
    iin=-N/2:1:N/2;                % array used for indexing
    iik=-K/2:1:K/2;                % array used for indexing
    D=sqrt((x-iin*Delta).^2+y^2);  % distance to each antenna element
    alphad=lambda./((4*pi*D));     % channel amplitudes        
    for s=1:sims   
        psi=rand(1)*2*pi;              % channel phase        
        Y=zeros(K+1,N+1);
        % generate QPSK data
        si=2*(round(rand(1,K+1)))-1;
        sq=2*(round(rand(1,K+1)))-1;        
        signals=sqrt(P/2).*(si+1i*sq);        
        S=diag(signals);
        for n=-N/2:N/2
            ni=n+N/2+1;
            % general model:
            Y(:,ni)=alphad(ni)*exp(1i*psi)*signals.'.*exp(-1i*2*pi/lambda*(D(ni)-D(N/2+1))-1i*2*pi/lambda*iik'*rf*(D(ni)-bias))+(randn(K+1,1)+1i*randn(K+1,1))*sqrt(N0/2);            
            % nearfield, narrowband model:
            %Y(:,ni)=alphad(ni)*exp(j*psi)*sqrt(P').*exp(+j*2*pi/lambda*(ni*Delta*cos(theta))-j*2*pi/lambda*iik'*rf*(D(N/2+1)-bias))+(randn(K+1,1)+j*randn(K+1,1))*sqrt(N0/2);                                      
        end        
        [x_hat,B_hat]=estimateLocation((inv(S'*S))*S'*Y,lambda,Delta,W,c,useSubArray,d);        
        LocError(l,s)=norm(xUE-x_hat);
        if (isempty(B_hat))
            BiasError(l,s)=+Inf;
        else
            BiasError(l,s)=abs(bias-B_hat);
        end                
    end            
end

disp('plot results')
figure(1)
h=loglog(xaxisvec,mean(LocError'),'r-',xaxisvec,mean(BiasError'),'b--');
grid
set(gca,'FontSize',12);
set(h,'Linewidth',2,'MarkerSize',8);
xl=xlabel('distance [m]');
yl=ylabel('RMSE [m]');
l=legend('UE location','UE bias');
set(xl,'Interpreter','latex','FontSize',12);
set(yl,'Interpreter','latex','FontSize',12);
set(l,'Interpreter','latex','FontSize',12,'Location','SouthEast');
pbaspect([2 1 1])
set(gcf, 'Color', 'w');
