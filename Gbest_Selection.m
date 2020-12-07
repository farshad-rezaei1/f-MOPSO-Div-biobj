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
function [g_elem]=Gbest_Selection(z1_ppbest,z2_ppbest)
np=size(z1_ppbest,1);
part=zeros(np);
zz11=zeros(np);
zz22=zeros(np);
flag=ones(1,np);
div=zeros(np);
ss=zeros(np);

% Preserving the Pbest Particles Once Detected to be Non-dominated
for j=1:np
    if j~=1
        for y=1:j-1
            if flag(1,y)==1
                if and(z1_ppbest(j)<=z1_ppbest(y),z2_ppbest(j)<=z2_ppbest(y))
                    flag(1,y)=0;
                elseif and(z1_ppbest(j)>=z1_ppbest(y),z2_ppbest(j)>=z2_ppbest(y)) && ...
                        (z1_ppbest(j)~=z1_ppbest(y) || z2_ppbest(j)~= z2_ppbest(y))
                    flag(1,j)=0;
                    break
                end
            end
        end
    end
end

% Descriminating the Non-dominated Pbests
d=0;
for j=1:np
    if flag(1,j)==1
        d=d+1;
        part(d)=j;
        zz11(d)=z1_ppbest(j);
        zz22(d)=z2_ppbest(j);
    end
end

% Evaluating the Diversity of the Pbests
for i=1:d
    min_ss=inf;
    for j=1:d
        ss(j)=abs(zz11(i)-zz11(j))+abs(zz22(i)-zz22(j));
        if i~=j && ss(j)<min_ss
            min_ss=ss(j);
        end
    end
    div(i)=(d~=1)*min_ss+d==1*inf;
end

% Selecting the Most-diversified Pbest as the Gbest Particle
max_div=-inf;
for i=1:d
    if div(i)>=max_div
        max_div=div(i);
        t=part(i);
    end
end
g_elem=t; % g_elem is the position of the Gbest particle in the current population 
end