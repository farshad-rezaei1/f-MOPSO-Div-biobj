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
function [pp,pv]=Initialization(np,nx,varmax,varmin,velmax,velmin)
pp=zeros(np,nx); % Particles' positions
pv=zeros(np,nx); % Particles' velocities

for j=1:np
    pp(j,1:nx)=(varmax-varmin).*rand(1,nx)+varmin; % Initializing the Particles' position
    pv(j,1:nx)=(velmax-velmin).*rand(1,nx)+velmin; % Initializing the Particles' velocity
end
end