function PEB=getPEBStandard(lambda,EK,EN,Delta,theta,rf,gamma,T)
    % computes the PEB with and without synchronization for the standard
    % case: far-field and narrowband
    % 
    J1=(2*pi/lambda)^(-2)*EK(1)*EN(1)*[1 0 0 0]'*[1 0 0 0];
    J2=EK(1)*EN(3)*Delta^2*(sin(theta))^2*[0 0 1 0]'*[0 0 1 0];
    J3=EK(3)*EN(1)*rf^2*[0 1 0 -1]'*[0 1 0 -1];
    J=gamma*(J1 + J2 + J3);        % FIM from the paper        
    JP=T'*J*T;                    % transformed to position domain
    
    % now compute PEB for case with and without synchronization
    
    tmp=inv(JP(1:3,1:3));    
    PEB(1)=sqrt(trace(tmp(2:3,2:3)));
    
    tmp=inv(JP);    
    PEB(2)=sqrt(trace(tmp(2:3,2:3)));