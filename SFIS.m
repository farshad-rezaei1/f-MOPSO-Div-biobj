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
function [dom,dom_ppbest]=SFIS(np_archive,zz1,zz2,z1_ppbest,z2_ppbest,z1_archive,z2_archive,z1_min,z1_max,z2_min,z2_max,st,it)
np=size(zz1,1);
z11=zeros(np);
z22=zeros(np);
z=zeros(4/st+1);
mu=zeros(4/st+1);
dom=zeros(np);
dom_ppbest=zeros(np);

% Deriving the Statistical Parameters
sigma_z1_high=std(z1_archive(1,1:floor(np_archive/3)))*(1/(z1_max-z1_min));
sigma_z1_mid=std(z1_archive(1,floor(np_archive/3)+1:2*floor(np_archive/3)))*(1/(z1_max-z1_min));
sigma_z1_low=std(z1_archive(1,2*floor(np_archive/3)+1:np_archive))*(1/(z1_max-z1_min));
mu_z1_high=(mean(z1_archive(1,1:floor(np_archive/3)))-z1_min)/(z1_max-z1_min);
mu_z1_mid=(mean(z1_archive(1,floor(np_archive/3)+1:2*floor(np_archive/3)))-z1_min)/(z1_max-z1_min);
mu_z1_low=(mean(z1_archive(1,2*floor(np_archive/3)+1:np_archive))-z1_min)/(z1_max-z1_min);
sigma_z2_high=std(z2_archive(1,1:floor(np_archive/3)))*(1/(z2_max-z2_min));
sigma_z2_mid=std(z2_archive(1,floor(np_archive/3)+1:2*floor(np_archive/3)))*(1/(z2_max-z2_min));
sigma_z2_low=std(z2_archive(1,2*floor(np_archive/3)+1:np_archive))*(1/(z2_max-z2_min));
mu_z2_high=(mean(z2_archive(1,1:floor(np_archive/3)))-z2_min)/(z2_max-z2_min);
mu_z2_mid=(mean(z2_archive(1,floor(np_archive/3)+1:2*floor(np_archive/3)))-z2_min)/(z2_max-z2_min);
mu_z2_low=(mean(z2_archive(1,2*floor(np_archive/3)+1:np_archive))-z2_min)/(z2_max-z2_min);    

% Applying Fuzzy Inference System to Obtain the Comprehensive Dominance Index
if it==1
    iii=1;
else
    iii=2;
end
j = 1;
while j<=np
    for ii=1:iii
        if ii==1
            z11(j)=(zz1(j)-z1_min)/(z1_max-z1_min);
            z22(j)=(zz2(j)-z2_min)/(z2_max-z2_min);
        elseif ii==2
            z11(j)=(z1_ppbest(j)-z1_min)/(z1_max-z1_min);
            z22(j)=(z2_ppbest(j)-z2_min)/(z2_max-z2_min);
        end
    w1=0; w2=1;    
    sum2=0; sum1=0; i=1;
    while i<=(4/st+1-4)
        if w1<0.5
            mf1=1/(1+exp(4/(sigma_z1_high*sqrt(exp(1)))*(z11(j)-mu_z1_high)));
            mf2=1/(1+exp(4/(sigma_z2_mid*sqrt(exp(1)))*(z22(j)-mu_z2_mid)));
            mu(i)=mf1*mf2;
            z(i)=w1*z11(j)+w2*z22(j);
            sum2=sum2+mu(i)*z(i);sum1=sum1+mu(i);
            i=i+1;
            mf1=1/(1+exp(4/(sigma_z1_high*sqrt(exp(1)))*(z11(j)-mu_z1_high)));
            mf2=1/(1+exp(4/(sigma_z2_low*sqrt(exp(1)))*(z22(j)-mu_z2_low)));
            mu(i)=mf1*mf2;
            z(i)=w1*z11(j)+w2*z22(j);
            sum2=sum2+mu(i)*z(i);sum1 = sum1+mu(i);
            i=i+1;
            mf1=1/(1+exp(4/(sigma_z1_mid*sqrt(exp(1)))*(z11(j)-mu_z1_mid)));
            mf2=1/(1+exp(4/(sigma_z2_low*sqrt(exp(1)))*(z22(j)-mu_z2_low)));
            mu(i)=mf1*mf2;
            z(i)=w1*z11(j)+w2*z22(j);
            sum2=sum2+mu(i)*z(i);sum1=sum1+mu(i);
            i=i+1;
            mf1=1/(1+exp(4/(sigma_z1_low*sqrt(exp(1)))*(z11(j)-mu_z1_low)));
            mf2=1/(1+exp(4/(sigma_z2_low*sqrt(exp(1)))*(z22(j)-mu_z2_low)));
            mu(i)=mf1*mf2;
            z(i)=w1*z11(j)+w2*z22(j);
            sum2=sum2+mu(i)*z(i);sum1=sum1+mu(i);
            i=i+1;
            w1=w1+st;w2=w2-st;
        elseif w1==0.5
            mf1=1/(1+exp(4/(sigma_z1_low*sqrt(exp(1)))*(z11(j)-mu_z1_low)));
            mf2=1/(1+exp(4/(sigma_z2_low*sqrt(exp(1)))*(z22(j)-mu_z2_low)));
            mu(i)=mf1*mf2;
            z(i)=w1*z11(j)+w2*z22(j);
            sum2=sum2+mu(i)*z(i);sum1=sum1+mu(i);
            i=i+1;
            w1=w1+st;w2=w2-st;
        elseif w1>0.5
            mf1=1/(1+exp(4/(sigma_z1_low*sqrt(exp(1)))*(z11(j)-mu_z1_low)));
            mf2=1/(1+exp(4/(sigma_z2_mid*sqrt(exp(1)))*(z22(j)-mu_z2_mid)));
            mu(i)=mf1*mf2;
            z(i)=w1*z11(j)+w2*z22(j);
            sum2=sum2+mu(i)*z(i);sum1=sum1+mu(i);
            i=i+1;
            mf1=1/(1+exp(4/(sigma_z1_low*sqrt(exp(1)))*(z11(j)-mu_z1_low)));
            mf2=1/(1+exp(4/(sigma_z2_high*sqrt(exp(1)))*(z22(j)-mu_z2_high)));
            mu(i)=mf1*mf2;
            z(i)=w1*z11(j)+w2*z22(j);
            sum2=sum2+mu(i)*z(i);sum1=sum1+mu(i);
            i=i+1;
            mf1=1/(1+exp(4/(sigma_z1_mid*sqrt(exp(1)))*(z11(j)-mu_z1_mid)));
            mf2=1/(1+exp(4/(sigma_z2_high*sqrt(exp(1)))*(z22(j)-mu_z2_high)));
            mu(i)=mf1*mf2;
            z(i)=w1*z11(j)+w2*z22(j);
            sum2=sum2+mu(i)*z(i);sum1 = sum1+mu(i);
            i=i+1;
            mf1=1/(1+exp(4/(sigma_z1_low*sqrt(exp(1)))*(z11(j)-mu_z1_low)));
            mf2=1/(1+exp(4/(sigma_z2_low*sqrt(exp(1)))*(z22(j)-mu_z2_low)));
            mu(i)=mf1*mf2;
            z(i)=w1*z11(j)+w2*z22(j);
            sum2=sum2+mu(i)*z(i);sum1 = sum1+mu(i);
            i=i+1;
            w1=w1+st;w2=w2-st;
        end
    end
    if ii==1
        dom(j)=sum2/(sum1+eps);
    elseif ii==2
        dom_ppbest(j)=sum2/(sum1+eps);
    end
    end
    j=j+1;
end
end