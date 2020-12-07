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
function [index]=Archive(index,z1_archive,z2_archive,np_archive)
% Preserving a certain solution in the external archive once detected to be
% a non-dominated solution and removing the ones previously recorded in the
% archive once detected to be dominated by the current one

for y=1:np_archive-1
    if index(1,y)==1
        if and(z1_archive(np_archive)<=z1_archive(y),z2_archive(np_archive)<=z2_archive(y))
            index(1,y)=0;
        elseif and(z1_archive(np_archive)>=z1_archive(y),z2_archive(np_archive)>=z2_archive(y)) && ...
                (z1_archive(np_archive)~=z1_archive(y) || z2_archive(np_archive)~=z2_archive(y))
            index(1,np_archive)=0;
            break
        end
    end
end
end