function [PEBGeneral PEBWideband]=getPEBGeneral(lambda,EK,EN,Delta,theta,rf,gamma,T,x,y,N,d,K,N0,P)
    Jm=zeros(4,4);        
    Jm2=Jm;    
    ii=-N/2:1:N/2;
    D=sqrt((x-ii*Delta).^2+y^2);
    alphad=lambda./((4*pi*D));      % distance-dependent channel gain
    for k=-K/2:K/2
        for n=-N/2:N/2
            dn=D(n+N/2+1);
            % general
            nabla1 = lambda/(2*pi);                           
            nabla2 = (1+k*rf)*(d-n*Delta*cos(theta))/dn-1;                
            nabla3=  (1+k*rf)*n*Delta*sin(theta)*d/dn;
            nabla4= -k*rf;
            nabla=[nabla1 nabla2 nabla3 nabla4];
            Jtemp=abs(alphad(n+N/2+1))^2/N0*(2*pi/lambda)^2*P(k+K/2+1)*nabla'*nabla;               
            Jm=Jm+Jtemp;

            % only wideband                 
            nabla1 = lambda/(2*pi);               
            nabla2 = (k*rf)*(d-n*Delta*cos(theta))/dn;                
            nabla3=  n*Delta*sin(theta);
            nabla4= -k*rf;
            nabla=[nabla1 nabla2 nabla3 nabla4];
            Jtemp=abs(alphad(n+N/2+1))^2/N0*(2*pi/lambda)^2*P(k+K/2+1)*nabla'*nabla;  
            Jm2=Jm2+Jtemp;            
        end
    end
    % general PEB
    J=Jm;
    JP=T'*J*T;           
    tmp=inv(JP(1:3,1:3));    
    PEBGeneral(1)=sqrt(trace(tmp(2:3,2:3)));    
    tmp=inv(JP);    
    PEBGeneral(2)=sqrt(trace(tmp(2:3,2:3)));
    
    
    % wideband PEB
    J=Jm2;
    JP=T'*J*T;           
    tmp=inv(JP(1:3,1:3));    
    PEBWideband(1)=sqrt(trace(tmp(2:3,2:3)));    
    tmp=inv(JP);    
    PEBWideband(2)=sqrt(trace(tmp(2:3,2:3))); 