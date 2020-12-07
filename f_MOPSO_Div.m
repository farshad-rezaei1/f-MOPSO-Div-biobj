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
function [index,nond_archive,z1_archive,z2_archive,fitness,np_archive,maxit_final]=f_MOPSO_Div(maxit,np,nx,st,k_max,k_min,varmin,varmax,velmin,velmax,stallit,stallpower,fobj) 
% This code is specifically developed for solving the bi-objective optimization problems 

% MOPSO Parameters
maxit_final=maxit;
pp_pbest=zeros(np,nx);
zz1=zeros(np);
zz2=zeros(np);
nond_archive=zeros(maxit*np,nx);
z1_archive=zeros(1,maxit*np);
z2_archive=zeros(1,maxit*np);
dom_archive=zeros(maxit);
z1_ppbest=zeros(np);
z2_ppbest=zeros(np);
fitness=zeros(maxit);
z1_min=inf;
z1_max=-inf;
z2_min=inf;
z2_max=-inf;
% p_max=1;
% p_min=1/nx;
% nm=20;
index=ones(1,maxit*np);
it=1;
np_archive=0;

% Initialization process of the algorithm
[pp,pv]=Initialization(np,nx,varmax,varmin,velmax,velmin);
for j=1:np
    np_archive=np_archive+1;
    x(1:nx)=pp(j,1:nx);
    
    % Function Evaluation   
    [z1,z2]=fobj(x);
    
    zz1(j)=z1;
    zz2(j)=z2;
    z1_ppbest(j)=z1;
    z2_ppbest(j)=z2;
    pp_pbest(j,1:nx)=pp(j,1:nx);
    if z1<z1_min
        z1_min=z1;
    elseif z1>z1_max
        z1_max=z1;
    elseif z2<z2_min
        z2_min=z2;
    elseif z2>z2_max
        z2_max=z2;
    end
    z1_archive(1,np_archive)=z1;
    z2_archive(1,np_archive)=z2;
    nond_archive(np_archive,1:nx)=x(1:nx);
    if np_archive~=1
        [index]=Archive(index,z1_archive,z2_archive,np_archive);
    end
end

% SFIS Implementation
[dom,dom_ppbest]=SFIS(np_archive,zz1,zz2,z1_ppbest,z2_ppbest,z1_archive,z2_archive,z1_min,z1_max,z2_min,z2_max,st,it);

% Determining the Global Best (Gbest) Particle
[g_elem]=Gbest_Selection(z1_ppbest,z2_ppbest);
pp_gbest(1:nx)=pp_pbest(g_elem,1:nx);
dom_gbest=dom_ppbest(g_elem);

% Sending the Gbest Particle to an Archive
z1_archive(it)=z1_ppbest(g_elem);
z2_archive(it)=z2_ppbest(g_elem);
dom_archive(it)=dom_gbest;

% Calculating the Mean Dominance Index of the Population
sumf=0;
for j=1:np
    sumf=sumf+dom(j);
end
fitness(1,it)=sumf/np;

% Main Loop

while it<maxit
    it = it+1;
    k=k_max-(k_max-k_min)*(it-2)/(maxit-2); 
%     disp(['Number of Iterations= ',num2str(it)]);
    for j = 1:np
        phi1=2.05*rand(1,nx);
        phi2=2.05*rand(1,nx);
        phi=phi1+phi2;
        khi=(2*k)./abs(2-phi-sqrt(phi.*(phi-4))); 
        
        % Update the velocity of the particles
        pv(j,1:nx)=khi.*(pv(j,1:nx)+phi1.*(pp_pbest(j,1:nx)-pp(j,1:nx))+phi2.*(pp_gbest(1:nx)-pp(j,1:nx)));
        
        % Return the velocity of the particles if going beyond the velocity boundaries
        flag4lbv=pv(j,:)<velmin(1,:);
        flag4ubv=pv(j,:)>velmax(1,:);
        pv(j,:)=(pv(j,:)).*(~(flag4lbv+flag4ubv))+velmin.*flag4lbv+velmax.*flag4ubv;
        
        % Update the position of the particles
        pp(j,:)=pp(j,:)+pv(j,:);
        
%         % Imposing Mutation
%         for i=1:nx
%             rm=rand(1,1);
%             pm=p_max-((p_max-p_min)/maxit)*it;
%             if rm<pm
%                 r=rand(1,1);
%                 if r<0.5
%                     delta=(2*r)^(1/(nm + 1))-1;
%                     pp(j,i)=pp(j,i)+delta*(varmax(1,i)-varmin(1,i));
%                 else
%                     delta=1-((2*(1 - r))^(1/(nm + 1)));
%                     pp(j,i)=pp(j,i)+delta*(varmax(1,i)-varmin(1,i));
%                 end    
%             end
%         end
        
        % Return the position and velocity of the particles if going beyond the position boundaries
        flag4lbp=pp(j,:)<varmin(1,:);
        flag4ubp=pp(j,:)>varmax(1,:);
        pp(j,:)=(pp(j,:)).*(~(flag4lbp+flag4ubp))+varmin.*flag4lbp+varmax.*flag4ubp; 
        pv(j,:)=(pv(j,:)).*(ones(1,nx)-2*(flag4lbp+flag4ubp));        
        np_archive=np_archive+1;
        x(1:nx)=pp(j,1:nx); 
        
        % Function Evaluation   
        [z1,z2]=fobj(x);
        
        zz1(j)=z1;
        zz2(j)=z2;
        if z1<z1_min
            z1_min=z1;
        elseif z1>z1_max
            z1_max=z1;
        elseif z2<z2_min
            z2_min=z2;
        elseif z2>z2_max
            z2_max=z2;
        end
        z1_archive(1,np_archive)=z1;
        z2_archive(1,np_archive)=z2;
        nond_archive(np_archive,1:nx)=x(1:nx);
        if np_archive~=1
            [index]=Archive(index,z1_archive,z2_archive,np_archive);
        end
    end
    
    % SFIS Implementation
    [dom,dom_ppbest]=SFIS(np_archive,zz1,zz2,z1_ppbest,z2_ppbest,z1_archive,z2_archive,z1_min,z1_max,z2_min,z2_max,st,it);
    
    % Determining the Personal Best Particles
    for j=1:np
        if dom(j)<dom_ppbest(j)
            dom_ppbest(j)=dom(j);
            z1_ppbest(j)=zz1(j);
            z2_ppbest(j)=zz2(j);
            pp_pbest(j,1:nx)=pp(j,1:nx);
        end
    end
    
    % Determining the Global Best (Gbest) Particle
    [g_elem]=Gbest_Selection(z1_ppbest,z2_ppbest);
    pp_gbest(1:nx)=pp_pbest(g_elem,1:nx);
    dom_gbest=dom_ppbest(g_elem);
     
    % Sending the Gbest Particle to an Archive    
    z1_archive(it)=z1_ppbest(g_elem);
    z2_archive(it)=z2_ppbest(g_elem);
    dom_archive(it)=dom_gbest;
    
    % Calculating the Mean Dominance Index of the Population
    sumf=0;
    for j=1:np
        sumf=sumf+dom(j);
    end
    fitness(1,it)=sumf/np;
    
    % Termination Criterion
    if it>=stallit
        if and(abs(z1_archive(it)-z1_archive(it-stallit+1))<10^(-stallpower),abs(z2_archive(it)-z2_archive(it-stallit+1))<10^(-stallpower))
            maxit_final=it;
            break
        end
    end
end