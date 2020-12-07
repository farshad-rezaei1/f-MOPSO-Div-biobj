%________________________________________________________________________%
% Diversity-enhanced fuzzy Multi-Objective Particle Swarm Optimizer      %
% (f-MOPSO/Div)                                                          %
%                                                                        %
% Developed in MATLAB R2018b                                             %
%                                                                        %
% Author and programmer: Farshad Rezaei, PhD                             %
%                                                                        %
% e-Mail: farshad.rezaei@gmail.com                                       %
%         f.rezaei@alumni.iut.ac.ir                                      %
%                                                                        %
% Homepage: https://www.linkedin.com/in/farshad-rezaei-5a92559a/         %
%                                                                        %
% Main paper: Rezaei, F., Safavi, H.R.,(2020) "f-MOPSO/Div: an improved  %
% extreme-point-based multi-objective PSO algorithm applied to a         %
% socio-economic-environmental conjunctive water use problem",           %
% Environmental Monitoring and Assessment, 192(12): 1-27                 %
%________________________________________________________________________%
function [lb,ub,nx,fobj] = Objective_Function(F)

switch F
    case 'F1'
        fobj=@Constr_Ex;
        nx=2;
        lb=[0.1 0];
        ub=[1 5];
        
    case 'F2'
        fobj=@Fonseca;
        nx=10;
        lb=(-4)*ones(1,nx);
        ub=4*ones(1,nx);
        
    case 'F3'
        fobj=@Kita;
        nx=2;
        lb=0*ones(1,nx);
        ub=7*ones(1,nx);
        
    case 'F4'
        fobj=@Binh;
        nx=2;
        lb=[0 0];
        ub=[5 3];
              
    case 'F5'
        fobj=@DTLZ2;
        nx=11;
        lb=0*ones(1,nx);
        ub=1*ones(1,nx);
        
    case 'F6'
        fobj=@Kursawe;
        nx=3;
        lb=-5*ones(1,nx);
        ub=5*ones(1,nx);
                
    case 'F7'
        fobj=@Deb;
        nx=2;
        lb=0.1*ones(1,nx);
        ub=1*ones(1,nx);  
end

% F1

function [z1,z2]=Constr_Ex(x)
m1=100;
m2=100;    
v1=x(1);
v2=x(2);
z1=v1+m1*abs(min(v2+9*v1-6,0))+m2*abs(min(-v2+9*v1-1,0));
z2=(1+v2)/v1+m1*abs(min(v2+9*v1-6,0))+m2*abs(min(-v2+9*v1-1,0));    
end

% F2

function [z1,z2]=Fonseca(x)
sum1=0;sum2=0;
for i=1:nx
    sum1=sum1+(x(i)-1/sqrt(nx))^2;
    sum2=sum2+(x(i)+1/sqrt(nx))^2;
end
z1=1-exp(-sum1);
z2=1-exp(-sum2);
end

% F3

function [z1,z2]=Kita(x)
m=100;  
pen1=m*(1+sign((1/6)*x(1)+x(2)-13/2));
pen2=m*(1+sign((1/2)*x(1)+x(2)-15/2));
pen3=m*(1+sign(5/x(1)+x(2)-30));
pen4=m*(2-2*sign(x(1)));
pen5=m*(2-2*sign(x(2)));
z1=-(-x(1)^2+x(2))+pen1+pen3+pen4;
z2=-(0.5*x(1)+x(2)+1)+pen2+pen3+pen5; 
end

% F4

function [z1,z2]=Binh(x)
m1=1000;
m2=1000;
v1=x(1);
v2=x(2);
z1=4*v1^2+4*v2^2+m1*abs(max((v1-5)^2+v2^2-25,0)) + ...
    m2*abs(min((v1-8)^2+(v2+3)^2-7.7,0));
z2=(v1-5)^2+(v2-5)^2+m1*abs(max((v1-5)^2+v2^2-25,0)) + ...
    m2*abs(min((v1-8)^2+(v2+3)^2-7.7,0));
end

% F5

function [z1,z2]=DTLZ2(x)
sum3=0;
for i=2:11
    sum3=sum3+(x(i)-0.5)^2;
end
z1=(1+sum3)*cos(x(1)*(pi/2));
z2=(1+sum3)*sin(x(1)*(pi/2));
end

% F6

function [z1,z2]=Kursawe(x)
sum1=0;sum2=0;
for i=1:2
    sum1=sum1+(-10*exp(-0.2*sqrt(x(i)^2+x(i+1)^2)));
end
for i=1:3
    sum2=sum2+((abs(x(i)))^0.8+5*sin(x(i)^3));
end
z1=sum1;
z2=sum2;
end

% F7

function [z1,z2]=Deb(x)
z1=x(1);
z2=(2-exp(-((x(2)-0.2)/0.004)^2)-0.8*exp(-((x(2)-0.6)/0.4)^2))/x(1);
end
end