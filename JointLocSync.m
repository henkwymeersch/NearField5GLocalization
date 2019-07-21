close all
clear all;
% Near-Field Joint Localization and Synchronization
% (c) Henk Wymeersch, 2019
% Usage: this code generates position estimates



% System parameters
% -----------------
K =256;                   % number of subcarriers - 1
N = 128;                    % number of antennas -1
fc = 28;                    % carrier [GHz]
c = 0.3;                    % speed of light [m/ns]
lambda = c/fc;              % carrier wavelength [m]                
Delta = lambda/2;           % antenna spacing in [m] -> setting too large creates sidelobes
W = 0.1;                    % bandwidth [GHz]
Pt=1;                       % mW
N0 = 290*1.381e-23*1e3*1e9; % noise PSD in mW/GHz (290 Kelvin * Boltzmann constant in W/Hz)
bias=100;                     % UE clock bias in m
P = Pt/W*ones(1,K+1);          % energy per subcarrier   
Deltaf=W/(K+1);                % subcarrier spacing
rf=Deltaf/fc;                  % ratio 
visualize=1;                % visualize the FFTs and the position estimates


steps=50;                       % number of distances we will consider
sims=100;                       % number of realization per distance, set to larger value for smoother curve
xaxisvec=logspace(-1,2,steps);  % vector of distances
disp('start simulation')
iin=-N/2:1:N/2;                % array used for indexing
iik=-K/2:1:K/2;                % array used for indexing
for l=1:steps
    disp(['step ' num2str(l) ' of '  num2str(steps) ' completed.']);
    
    % Generate the received signal
    % ----------------------------           
    for s=1:sims   
        d=xaxisvec(l);                  % distance to UE
        theta=pi/2+(2*rand(1)-1)*pi/4;        % AOA            
        x=d*cos(theta);                 % x coordinate
        y=d*sin(theta);                 % y coordinate    
        xUE=[x y];
        D=sqrt((x-iin*Delta).^2+y^2);  % distance to each antenna element
        alphad=lambda./((4*pi*D));     % channel amplitudes        
        bias=rand(1)*100;                     % UE clock bias in m
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
        
        % using sub-array approach
        useSubArray=1;
        x_hat=[];
        B_hat=[];
        [x_hat,B_hat]=estimateLocation((inv(S'*S))*S'*Y,lambda,Delta,W,c,useSubArray,d,0,xUE,visualize);                        
        
        if (isempty(x_hat))
            LocError(l,s)=+Inf;
        else
            LocError(l,s)=norm(xUE-x_hat);
        end                  
        if (isempty(B_hat))
            BiasError(l,s)=+Inf;
        else
            BiasError(l,s)=abs(bias-B_hat);
        end                
        
        % using conventional approach, needs to know the bias
        [x_hat2,~]=estimateLocation((inv(S'*S))*S'*Y,lambda,Delta,W,c,0,d,bias,xUE,visualize);  
        LocError2(l,s)=norm(xUE-x_hat2);
        
    end            
end

disp('plot results')
figure(1)
h=loglog(xaxisvec,mean(LocError'),'b-',xaxisvec,mean(BiasError'),'g--',xaxisvec,mean(LocError2'),'r-');
grid
set(gca,'FontSize',12);
set(h,'Linewidth',2,'MarkerSize',8);
xl=xlabel('distance [m]');
yl=ylabel('RMSE [m]');
l=legend('UE location - sub-array','UE bias - sub-array','UE location - standard');
set(xl,'Interpreter','latex','FontSize',12);
set(yl,'Interpreter','latex','FontSize',12);
set(l,'Interpreter','latex','FontSize',12,'Location','SouthEast');
pbaspect([2 1 1])
set(gcf, 'Color', 'w');
