fprintf('We have taken the angle of attack as zero for ensuring simplified calculations \n');
Vinf = input(' Enter the freestream velocity: \n');
R = input(' Enter the radius: \n');
n= input('Enter the number of panels \n');
%angle substended by each panel at center
dtheta=2*pi/n;
%central angle and values of theta to give values of x1 y1 and so on
theta=pi+pi/n:-dtheta:-pi+pi/n;
%Boundary Points of Panel
X=R*cos(theta);
Y=R*sin(theta);
%we will compute the potential due to 
%contx and conty are the control points of the panel
%S is the length of each panel
%Beta is the angle between the freestream velocity direction & normal of the panel
for i=1:n
 phi(i)= atan2((Y(i+1)-Y(i)),(X(i+1)-X(i)));
 beta(i)=phi(i)+pi/2;
 contx(i)=(X(i+1)+X(i))/2;
 conty(i)=(Y(i+1)+Y(i))/2;
 S(i)=sqrt((Y(i+1)-Y(i))^2+(X(i+1)-X(i))^2);
end
%We will compute the required integrals for computing the potential at ith panel 
for i=1:n
    for j=1:n
        if i == j
            I(i,j)=0;
        else
        A=-(contx(i)-X(j))*cos(phi(j))-(conty(i)-Y(j))*sin(phi(j));
        B=(contx(i)-X(j))^2+(conty(i)-Y(j))^2;
        C=sin(phi(i)-phi(j));
        D=(conty(i)-Y(j))*cos(phi(i))-(contx(i)-X(j))*sin(phi(i));
        E=sqrt(B-A^2);
        I(i,j)=C/2*log((S(j)^2+2*A*S(j)+B)/B)+(D-A*C)/E*(atan2((S(j)+A),E)-atan2(A,E));
        J(i,j)=(D-A*C)/2/E*log((S(j)^2+2*A*S(j)+B)/B) -C*(atan2((S(j)+A),E)-atan2(A,E));
        end
    Vn(i,1)=Vinf*cos(beta(i)); %Component of Vinf normal to the ith panel is
    end
end
M=(I/(2*pi))+(eye(n)/2);   %eye(n) stand for identity matrix of order n
lambda=-inv(M)*Vn;

% check whether the sum of source strenghs is 0
fprintf('The sum of all sources is %f',dot(lambda,S));  
Vs = Vinf*sin(beta) + J*lambda/2/pi; % tangent velocity
Cp = 1-(Vs/Vinf).^2;                  % pressure coefficient
%Plot
close all
angle = 0:0.01:2*pi;              % angle for plotting  
Cp_exact = 1 - 4*sin(angle).^2;  % pressure coefficient, analytical solution
Vs_exact = 2*Vinf*sin(angle);   % tangent velocity, analytical solution

figure

subplot(1,2,1);
plot(R*cos(0:0.01:2*pi),R*sin(0:0.01:2*pi),'b', X,Y,'k',contx,conty,'om'); 
axis equal; legend('Shape of the Body','Panel distribution','Control points')
xlabel('x, m'); ylabel('y, m'); title('Fig. 1 Source Panel Distribution');

subplot(1,2,2);
plot(angle,Cp_exact,'g',beta,Cp,'*c');
grid;
legend('C_p (exact', 'C_p (Source Panel Method )');