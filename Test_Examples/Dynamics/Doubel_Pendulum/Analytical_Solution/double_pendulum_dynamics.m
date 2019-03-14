function dy = double_pendulum_dynamics(t,y,m1,l1,m2,l2)


% System of ODEs for a double pendulum (mass m and link length L)
g = 9.8;

th1 = y(1);       % angle 1
th2 = y(2);       % angle 2
thd1 = y(3);           % angle 1 dot
thd2 = y(4);           % angle 2 dot

% The derivatives
dy(1) = thd1;
dy(2) = thd2;
dy(3) = -(m2*l1*l2*thd2^2*sin(th1-th2)/2+3*m2*l1^2*thd1^2*sin(th1-th2)*cos(th1-th2)/4-3*m2*l1*g*sin(th2)*cos(th1-th2)/4+m1*g*l1*sin(th1)/2+m2*g*l1*sin(th1))/(m1*l1^2/3+m2*l1^2-3*m2*l1^2*cos(th1-th2)*cos(th1-th2)/4);
dy(4) = ((3*l1)/(2*l2))*thd1^2*sin(th1-th2)-((3*g)/(2*l2))*sin(th2)-((3*l1)/(2*l2))*dy(3)*cos(th1-th2);
dy = dy(:);