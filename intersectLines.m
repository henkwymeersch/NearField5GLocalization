function x_hat = intersectLines(loc1,loc2,theta1,theta2)
% computes the intersection of a line starting at [loc1(1) 0] with AOA theta1 and a
% line starting at [loc2(1) 0] with AOA theta2
    T1=tan(theta1);
    T2=tan(theta2);
    x=(loc1(1)*T1-loc2(1)*T2)/(T1-T2);
    y=x*T1-loc1(1)*T1;
    x_hat=[x y];