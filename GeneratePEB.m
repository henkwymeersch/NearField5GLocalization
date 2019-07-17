close all
clear all;
clc
% Near-Field Joint Localization and Synchronization
% (c) Henk Wymeersch, 2019
% Usage: this code generates position error bounds


freeparam='Delta';              % can be W, Delta, distance, N, K, 
K = 1024;                       % number of subcarriers - 1
N = 128;                        % number of antennas -1
fc = 28;                        % carrier [GHz]
c = 0.3;                        % speed of light [m/ns]
lambda = c/fc;                  % carrier wavelength [m]                
Delta = lambda/2;               % antenna spacing in [m]
W = 0.1;                        % bandwidth [GHz]
Pt=1;                           % mW
N0 = 290*1e3*1.381e-23*1e9;     % noise PSD in mW/GHz


% determine the variable to sweep over: 
steps=25;
xvec=1;
yvec=8;
switch (freeparam)
    case 'W'
        Wvec=logspace(-4,0,steps);
        xaxisvec=Wvec;
    case 'Delta'
        Deltavec=lambda*(logspace(-1,2,steps));        
        xaxisvec=Deltavec/lambda;        
    case 'N'
        Nvec=2.^(floor(linspace(1,13,steps)));       
        xaxisvec=Nvec;
    case 'distance'
        xvec=logspace(-1,2,steps);
        yvec=logspace(0,2,steps);
        xaxisvec=sqrt(xvec(1).^2+yvec.^2);
    case 'K'
        Kvec=2.^(floor(linspace(1,13,steps)));       
        xaxisvec=Kvec;
end

xUE = [xvec(1),yvec(1)];           % user location [m,m]
disp('start simulation')
for l=1:steps    
    disp(['step ' num2str(l) ' of '  num2str(steps) ' completed.']);
    switch (freeparam)
        case 'W'
            W=Wvec(l);
        case 'Delta'
            Delta=Deltavec(l);
        case 'N'
            N=Nvec(l);
        case 'distance'            
            xUE(2)=yvec(l);
        case 'K'
            K=Kvec(l);
            W=Deltaf*(K+1);
    end        
    d=norm(xUE);                    % distance between array center and user        
    x=xUE(1);                       % x coordinate 
    y=xUE(2);                       % y coordinate
    theta=acos(x/d);                % AOA
    P = Pt/W*ones(1,K+1);          % energy per subcarrier 
    iin=-N/2:1:N/2;                % array used for indexing
    iik=-K/2:1:K/2;                % array used for indexing   
    for m=1:6    
        EK(m)=sum(P.*iik.^(m-1));     
    end
    iin=-N/2:1:N/2;
    for m=1:6
        EN(m)=sum(iin.^(m-1));     
    end
    Deltaf=W/(K+1);         % subcarrier spacing
    rf=Deltaf/fc;           % ratio 
    alpha=lambda/((4*pi*d));
    gamma=abs(alpha)^2/N0*(2*pi/lambda)^2;
    
    T=[1 0 0 0; 0 x/d y/d 0 ; 0 -y/d^2 x/d^2 0; 0 0 0 1];   % Jacobian    
    
    % Case 1: FIM for standard case
    % -----------------------------   
    PEBstandard(:,l)=getPEBStandard(lambda,EK,EN,Delta,theta,rf,gamma,T);
    
    % Case 2: Near field
    % ------------------
    PEBNearField(:,l)=getPEBNearField(lambda,EK,EN,Delta,theta,rf,gamma,T,x,y,N,d,N0);
                        
    % Remaining Cases
    % ------------------
    [PEBGeneral(:,l), PEBWideband(:,l)]=getPEBGeneral(lambda,EK,EN,Delta,theta,rf,gamma,T,x,y,N,d,K,N0,P);
            
end


disp('plot results')
figure(1)
h=loglog(xaxisvec,PEBstandard(1,:),'r+-',xaxisvec,PEBNearField(1,:),'gs-',xaxisvec,PEBWideband(1,:),'b.-',xaxisvec,PEBGeneral(1,:),'k-',xaxisvec,PEBNearField(2,:),'gs--',xaxisvec,PEBWideband(2,:),'b.--',xaxisvec,PEBGeneral(2,:),'k--');
grid
set(gca,'FontSize',12);
set(h,'Linewidth',2,'MarkerSize',8);
switch (freeparam)
        case 'W'                               
            xl=xlabel('W [GHz]');
        case 'Delta'                        
            xl=xlabel('$\Delta$/$\lambda$');
        case 'N'                        
            xl=xlabel('# antennas []');
        case 'K'                        
            xl=xlabel('number of subcarriers []');
        case 'distance'
            xl=xlabel('distance [m]');
end
yl=ylabel('PEB [m]');
set(xl,'Interpreter','latex','FontSize',12);
set(yl,'Interpreter','latex','FontSize',12);
l=legend('standard model, $B$ known','near-field narrowband, $B$ known', 'far-field wideband, $B$ known','general model, $B$ known','near-field narrowband, $B$ unknown', 'far-field wideband, $B$ unknown','general model, $B$ unknown');
set(l,'Interpreter','latex','FontSize',12,'Location','NorthWest');
pbaspect([2 1 1])
set(gcf, 'Color', 'w');
