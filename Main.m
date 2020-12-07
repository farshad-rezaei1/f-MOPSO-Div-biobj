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
% f-MOPSO/Div algorithm 
clc
clear
close all
tic
run=10; % Maximum number of the algorithm runnings conducted
maxit=200; % Maximum number of the iterations
stallit=50; % Number of the iterations over which the global guides are not getting better
np=40; % Population size
Function_name='F5'; % Name of the test problem suite which can be replaced with any desired problem
[lb,ub,nx,fobj]=Objective_Function(Function_name); % Load details of the selected benchmark function
varmin=lb; % Upper bound defined for the positions which can generally be a desired vector
varmax=ub; % Lower bound defined for the positions which can generally be a desired vector
limvel=0.1; % A ratio of the maximum distance in the search space to yield the maximum velocity
velmax=limvel*(varmax(1,1:nx)-varmin(1,1:nx)); % Upper bound defined for the velocities
velmin=-velmax; % Lower bound defined for the velocities
st=0.1; % Weight step which can be a round ratio of 1 and that the smallest this ratio, the better
stallpower=6; % Power of 10 denoting the maximum allowable error between the objective values over the stall iterations
index_main=zeros(run,maxit*np); % An index calculated to delineate if a solution is non-dominated
nond_archive_main=zeros(run,maxit*np,nx); % Elements of a solution stored in the non-dominated archive
z1_archive_main=zeros(run,maxit*np); % First objective of a solution in the non-dominated archive
z2_archive_main=zeros(run,maxit*np); % Second objective of a solution in the non-dominated archive
fitness_main=zeros(run,maxit); % Average Dominance Indices of the population at each iteration
k_max=0.9; % Maximum value of k used to tune the constriction coefficient at each iteration
k_min=0.4; % Minimum value of k used to tune the constriction coefficient at each iteration
pareto=zeros(maxit*np,nx); % Elements of a non-dominated solution stored in the archive
obj1=zeros(1,maxit*np); % First objective of a non-dominated solution stored in the archive
obj2=zeros(1,maxit*np); % Second objective of a non-dominated solution stored in the archive
x1=zeros(maxit); % values of the X axis of one of the charts
y1=zeros(maxit); % Values of the Y axis of one of the charts
ss=zeros(1,maxit*np); % Pre-allocating ss values
d=zeros(1,maxit*np); % Pre-allocating d values
mean_ideal_dis=zeros(run); % Pre-allocating for speed
spacing=zeros(run); % Pre-allocating for speed
no_solutions=zeros(run); % Pre-allocating for speed
ideal_point=zeros(run,2); % Pre-allocating for speed
for nrun=1:run
    [index,nond_archive,z1_archive,z2_archive,fitness,np_archive,maxit_final]=f_MOPSO_Div(maxit,np,nx,st,k_max,k_min,varmin,varmax,velmin,velmax,stallit,stallpower,fobj);
    index_main(nrun,:)=index(1,:);
    nond_archive_main(nrun,:,1:nx)=nond_archive(:,1:nx);
    z1_archive_main(nrun,:)=z1_archive(1,:);
    z2_archive_main(nrun,:)=z2_archive(1,:);
    fitness_main(nrun,:)=fitness(1,:);
    
    % Displaying Final Results and the Chart
    d1=0;
    for yy = 1:np_archive
        if index(1,yy)==1
            d1=d1+1;
            pareto(d1,:)=nond_archive_main(nrun,yy,:);
            obj1(1,d1)=z1_archive_main(nrun,yy);
            obj2(1,d1)=z2_archive_main(nrun,yy);
            if d1==1
                obj1_min=obj1(1,d1);
                obj2_min=obj2(1,d1);
            elseif obj1(1,d1)<obj1_min
                obj1_min=obj1(1,d1);
            elseif obj2(1,d1)<obj2_min
                obj2_min=obj2(1,d1);
            end
        end
    end
    ideal1=min(obj1(1,1:d1));
    ideal2=min(obj2(1,1:d1));
    
    % The Chart of the Pareto-front
    x2=obj1(1,1:d1);
    y2=obj2(1,1:d1);
    figure(nrun);
    plot(x2,y2,'ro')
    xlabel('Z1');
    ylabel('Z2');
    hold on
    
    % The Chart of the Dominance Index (DI) vs. Iterations
%     for r=1:maxit_final
%         x1(r)=r;
%         y1(r)=max(fitness_main(nrun,r),0);
%     end
%     figure(nrun+run);
%     semilogy(x1,y1,'-b')
%     xlabel('Iteration');
%     ylabel('Dominance Index');
%     hold on
%     disp('The Pareto = ');
%     for i = 1:d1
%         for j = 1:nx
%             disp(num2str(pareto(i,j)));
%         end
%     end
    sum1=0;
    for i=1:d1
        sum1=sum1+sqrt((obj1(1,i)-ideal1)^2+(obj2(1,i)-ideal2)^2);
    end
    mid=sum1/d1;
    for i=1:d1
        min_ss=inf;
        for j=1:d1
            ss(1,j)=abs(obj1(1,i)-obj1(1,j))+abs(obj2(1,i)-obj2(1,j));
            if i~=j && ss(1,j)<min_ss
                min_ss=ss(1,j);
            end
        end
        d(1,i)=min_ss;
    end
    d_average=mean(d(1,1:d1));
    sum2=0;
    for i=1:d1
        sum2=sum2+(d_average-d(1,i))^2;
    end
    s=sqrt(sum2/(d1-1));
    ns=d1;
    mean_ideal_dis(nrun)=mid;
    spacing(nrun)=s;
    no_solutions(nrun)=ns;
    ideal_point(nrun,1)=ideal1;
    ideal_point(nrun,2)=ideal2;
end
disp('Performance Metrics in total runs =');

% Computing the Performance Criteria

% MID = Mean Ideal Distance
disp(['Average MID = ',num2str(mean(mean_ideal_dis(1:run)))]); 
disp(['Best MID = ',num2str(min(mean_ideal_dis(1:run)))]); 
disp(['Worst MID = ',num2str(max(mean_ideal_dis(1:run)))]); 
disp(['Std of MID = ',num2str(std(mean_ideal_dis(1:run)))]); 

% S = Spacing
disp(['Average S = ',num2str(mean(spacing(1:run)))]); 
disp(['Best S = ',num2str(min(spacing(1:run)))]); 
disp(['Worst S = ',num2str(max(spacing(1:run)))]); 
disp(['Std of S = ',num2str(std(spacing(1:run)))]); 

% NS = Number of Solutions
disp(['Average NS = ',num2str(mean(no_solutions(1:run)))]); 
disp(['Best NS = ',num2str(max(no_solutions(1:run)))]); 
disp(['Worst NS = ',num2str(min(no_solutions(1:run)))]); 
disp(['Std of NS = ',num2str(std(no_solutions(1:run)))]); 

% Ideal Point
for i=1:2
    if i==1
        disp(['Average Minimum of the 1st Objective = ',num2str(mean(ideal_point(1:run,i)))]);
    elseif i==2
        disp(['Average Minimum of the 2nd Objective = ',num2str(mean(ideal_point(1:run,i)))]);
    end
end
toc