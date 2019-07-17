function PEB=getPEBNearField(lambda,EK,EN,Delta,theta,rf,gamma,T,x,y,N,d,N0)
    % computes the PEB with and without synchronization for the near-field
    % case
    % 
          
    A=zeros(4,4);
    ii=-N/2:1:N/2;
    D=sqrt((x-ii*Delta).^2+y^2);
    alphad=lambda./((4*pi*D));      % distance-dependent channel gain
    for i=1:4
        for j=1:4
            A(i,j)=sum(ii.^(i-1).*(d./D).^(j-1).*abs(alphad).^2/abs(alphad(N/2+1))^2);
        end
    end        
    
    J1=(2*pi/lambda)^(-2)*EK(1)*A(1,1)*[1 0 0 0]'*[1 0 0 0];
    J2=EK(1)*A(3,1)*Delta^2*(sin(theta))^2*[0 0 1 0]'*[0 0 1 0];
    J3=EK(3)*A(1,1)*rf^2*[0 1 0 -1]'*[0 1 0 -1];  
    
    t1= -Delta/d*cos(theta)*A(2,2)+A(1,2)-A(1,1);    
    t2=A(2,2)*Delta*sin(theta);
    J4=(lambda/(2*pi))*EK(1)*[0 t1 t2 0;t1 0 0 0;t2 0 0 0; 0 0 0 0];                                    
    t3=EK(1)*(A(1,1)+A(1,3)-2*(Delta/d*cos(theta)*A(2,3)+A(1,2)-Delta/d*cos(theta)*A(2,2)));        
    t3=EK(1)*(A(1,1)+A(1,3)-2*(Delta/d*cos(theta)*A(2,3)+A(1,2)-Delta/d*cos(theta)*A(2,2))+Delta^2*A(3,3)*(cos(theta)/d)^2);                
    t4=EK(1)*Delta*sin(theta)*((A(2,3)-A(2,2))-Delta*A(3,3)*cos(theta)/d);     
    J8=[0 0 0 0; 0 t3 t4 0; 0 t4 0 0; 0 0 0 0];
    J=abs(alphad(N/2+1))^2/N0*(2*pi/lambda)^2*(J4+J1+J2+J3+J8);
    JP=T'*J*T;       
    % now compute PEB for case with and without synchronization    
    tmp=inv(JP(1:3,1:3));    
    PEB(1)=sqrt(trace(tmp(2:3,2:3)));    
    tmp=inv(JP);    
    PEB(2)=sqrt(trace(tmp(2:3,2:3)));