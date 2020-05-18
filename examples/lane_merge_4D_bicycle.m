clear
close all

% example for a vehicle lane change model. Dynamics are 4D bicycle model
% found in vehicle_lange_merge_dyn.m
% the trajectories are parameterizes by a commanded velcoity and lateral
% position (k(1) and k(2))

%Author: Sean Vasov
%Date: May 17, 2020

%-------------------------------------------------------------------------%

m = [15,10,10,10];

T = 4;

grid_lower_bounds = [-10, -6, -0.3, 8];
grid_upper_bounds = [110, 6, 0.3 , 22];

h = (grid_upper_bounds-grid_lower_bounds)./m;

N = 1e3;

cfl_num = 0.8;

k_test = [19,4];


%% compute average velocities
disp('computing average velocities for cell')

fbar = compute_average_velocity_functions(@(z,k) vehicle_lane_merge_dyn(0,z,k),2,h);

tic
ubar = compute_average_velocities(k_test,fbar,m,grid_lower_bounds,grid_upper_bounds);
toc

dt = min([cfl_num*h(1)./max(abs(ubar{1})),...
        cfl_num*h(2)./max(abs(ubar{2})),...
        cfl_num*h(3)./max(abs(ubar{3})),...
        cfl_num*h(4)./max(abs(ubar{4}))]);

Ntstep = floor(T/dt);


%% compute initial condition for ego vehicle
mu = [0,0,0,15];
s =  [0.5,0.5,1*pi/180,0.5];



p = cell(1,floor(1/dt));

p{1} = [randRange(mu(1)-3*s(1),mu(1)+3*s(1),mu(1),s(1),1,N);...
        randRange(mu(2)-3*s(2),mu(2)+3*s(2),mu(2),s(2),1,N);...
        randRange(mu(3)-3*s(3),mu(3)+3*s(3),mu(3),s(3),1,N);...
        randRange(mu(4)-3*s(4),mu(4)+3*s(4),mu(4),s(4),1,N)];
    
disp('computing initial condition')
tic

ph = x_val_to_sub(p{1}',h,grid_lower_bounds);
idx0 = sub2indnd(m,ph);
Q0 = accumarray(idx0,1,[prod(m),1]);
Q0 = sparse(Q0/N/prod(h));

toc



%% matrix with max fluxes for upwind
disp(' precompute and store sparse matrices for some operations')
tic
[UP,DQI,DQi,DF] = build_matrices_for_fvm_advection(ubar,m);
toc

disp('advection beginning')
%% forward advect
tic
Q = fvm_advection(Q0,dt,h,ubar,UP,DQI,DQi,DF,T);
toc

disp(['initial volume: ',num2str(sum(sum(Q{1}))*prod(h))])
disp(['final volume: '  ,num2str(sum(sum(Q{end}))*prod(h))])

[meanQ,stdQ] = compute_mean_and_stddev(Q,m,grid_lower_bounds,grid_upper_bounds);


%% get obstacle as particles
mu_obs = [0,4,0,18];
s_obs = [0.1,0.1,0.01,0.1];

Qobs = cell(1,floor(T/dt));

pobs = cell(1,floor(T/dt));

pobs{1} = [randRange(mu_obs(1)-3*s_obs(1),mu_obs(1)+3*s_obs(1),mu_obs(1),s_obs(1),1,N); randRange(mu_obs(2)-3*s_obs(2),mu_obs(2)+3*s_obs(2),mu_obs(2),s_obs(2),1,N);...
           randRange(mu_obs(3)-3*s_obs(3),mu_obs(3)+3*s_obs(3),mu_obs(3),s_obs(3),1,N); randRange(mu_obs(4)-3*s_obs(4),mu_obs(4)+3*s_obs(4),mu_obs(4),s_obs(4),1,N)];
disp('computing obstacle condition')

tic

for i = 1:Ntstep
    
    subobs = x_val_to_sub(pobs{i}',h,grid_lower_bounds);
    
    idx0obs = sub2indnd(m,subobs);
    
    Qobs{i} = accumarray(idx0obs,1,[prod(m),1]);
    
    Qobs{i} = sparse(Qobs{i}/N/prod(h));
    
    if i < Ntstep
        pcur = pobs{i};
        pnext = NaN(size(pcur));
        for j = 1:N
            
            vobs   =   pcur(4,j);
            psiobs =   pcur(3,j);
            
            fobs = [vobs*cos(psiobs)+psiobs*1.5*sin(psiobs);vobs*sin(psiobs)-psiobs*1.5*cos(psiobs);-psiobs;atan(18-vobs)];
            
            pnext(:,j) = pcur(:,j)+dt*fobs;
            
        end
        pobs{i+1} = pnext;
    end

end

toc


%% monte carlo
disp('running monte carlo')
tic
mean_mc = NaN(length(Q),4);
std_mc = NaN(length(Q),4);

for j = 1:4
mean_mc(1,j) = mean(p{1}(j,:));
std_mc(1,j) =  std(p{1}(j,:));
end

for i = 1:Ntstep-1
    pttt = NaN(4,N);
    pit = p{i};
    for j = 1:N
%         [~,ptemp] = ode45(@(t,z) f(z),[0 dt],pit(:,j));
        pttt(:,j) = pit(:,j)+dt*vehicle_lane_merge_dyn(0,pit(:,j),k_test);
    end
    p{i+1} = pttt;
    
    for j = 1:4
        mean_mc(i+1,j) = mean(p{i+1}(j,:));
        std_mc(i+1,j) = std(p{i+1}(j,:));
    end

end

toc
%% caclulate probability of collision at each timestep

Pcol = zeros(1,length(Q));
Ecol = zeros(1,length(Q));
Pcolsample = zeros(1,length(Q));
Ecolsample = zeros(1,length(Q));
Pcolmc = zeros(1,length(Q));
Ecolmc = zeros(1,length(Q));

ftprint = [4,2];
ftprint_scale = ceil(ftprint./h(1:2));


mo = 1500;
mc = 1500;

ni = NaN(1,(floor(1/dt)));
disp('calculating probability and energy of collision with fvm')
tic
for n = 1:Ntstep
    
    ind_o = find(Qobs{n});
    sub_o = ind2subnd([m m m m],ind_o);
%     
    
    is_to_check = find(Q{n});
    
    sub_i = ind2subnd([m m m m],is_to_check);
    
    ni(n) = length(is_to_check);
    
    for i = 1:size(sub_i,1)
        

      %footprint
        L =  (abs(sub_o(:,1) - sub_i(i,1)) <= ftprint_scale(1)) & (abs(sub_o(:,2) - sub_i(i,2)) <= ftprint_scale(2));
        
         Pcol(n) = Pcol(n)+sum(Q{n}(is_to_check(i))*prod(h)^2*Qobs{n}(ind_o(L)));

                
        pl =  (sub_i(i,3)-1/2)*h(3)+grid_lower_bounds(3)+[-h(3)/2 h(3)/2];
        plo = (sub_o(L,3)-1/2)*h(3)+grid_lower_bounds(3)+[-h(3)/2 h(3)/2];
        
        vl =  (sub_i(i,4)-1/2)*h(4)+grid_lower_bounds(4)+[-h(4)/2 h(4)/2];
        vlo = (sub_o(L,4)-1/2)*h(4)+grid_lower_bounds(4)+[-h(4)/2 h(4)/2];
        
        
        %en approximately v^2-2*cos(psi-psio)*v*vo+vo^2
        pdiff = cos(pl(1)-plo(:,1))-cos(pl(1)-plo(:,2))-cos(pl(2)-plo(:,1))+cos(pl(2)-plo(:,2));
        
        
        Ecol(n) = Ecol(n)+sum(Q{n}(is_to_check(i))*Qobs{n}(ind_o(L)).*...
                    ((vlo(:,2).^3/3-vlo(:,1).^3/3)*prod(h(1:3))*prod(h)+(vl(2)^3/3-vl(1)^3/3)*prod(h(1:3))*prod(h)+...
                     -2*(vlo(:,2).^2/2-vlo(:,1).^2/2)*(vl(2)^2/2-vl(1)^2/2).*pdiff*prod(h(1:2))^2));
      
    end
    
end
toc

disp('calculating probability and energy of collision for monte carlo and sampling ')
tic
for n = 1:Ntstep
    for i = 1:N
           
        xn = p{n}(1,i);
        yn = p{n}(2,i);
        psin = p{n}(3,i);
        vn = p{n}(4,i);

        Lcol = abs(xn-pobs{n}(1,:)) < ftprint(1) & abs(yn-pobs{n}(2,:))<ftprint(2);
           
        Pcolmc(n) = Pcolmc(n)+nnz(Lcol)/N^2;
        
        ern = (vn*cos(psin)-pobs{n}(4,Lcol).*cos(pobs{n}(3,Lcol))).^2+(vn*sin(psin)-pobs{n}(4,Lcol).*sin(pobs{n}(3,Lcol))).^2;
           
        Ecolmc(n) = Ecolmc(n)+sum(ern)/N^2;
        
    end
end
toc

energy_mult = mc^2*mo/2/(mc+mo)^2;
Ecol =  energy_mult*Ecol;
Ecolmc = energy_mult*Ecolmc;


%% calculating and plotting statistics

tvec = (0:(Ntstep-1))*dt;

figure(1)
subplot(2,2,1)
plot(tvec,Pcol)
hold on
% plot((0:length(Q)-1)*dt*T,Pcolsample)
plot(tvec,Pcolmc,'o')

ylabel('P_{col}')
xlabel('time')
grid on

subplot(2,2,2)
plot(tvec,Pcol-Pcolmc)
hold on
% plot((0:length(Q)-1)*dt*T,Pcolsample-Pcolmc)
plot(0,0)

xlabel('time')
ylabel('Error in P_{col}')
grid on

subplot(2,2,3)
plot(tvec,Ecol)
hold on
% plot((0:length(Q)-1)*dt*T,Ecolsample)
plot(tvec,Ecolmc,'o')

xlabel('time')
ylabel('E_{col}')
grid on

subplot(2,2,4)
plot(tvec,Ecol-Ecolmc)
hold on
% plot((0:length(Q)-1)*dt*T,Ecolsample-Ecolmc)
plot(0,0)

grid on
xlabel('time')
ylabel('Error in E_{col}')



for i = 1:4
    
    figure(2)
    subplot(2,2,i)
    plot(tvec,meanQ(:,i))
    hold on
    plot(tvec,mean_mc(:,i),'o')
    
    figure(3)
    subplot(2,2,i)
    plot(tvec,stdQ(:,i))
    hold on
    plot(tvec,std_mc(:,i),'o')
    
end

figure(4)

xgrid2D = get_grid_points(m(1:2),grid_lower_bounds(1:2),grid_upper_bounds(1:2));


 Qprojxy = project_onto_dimension(Q,[1 2],m,grid_lower_bounds,grid_upper_bounds);
 Qobsprojxy = project_onto_dimension(Qobs,[1 2],m,grid_lower_bounds,grid_upper_bounds);
 
 XX = reshape(xgrid2D(:,1),m(1:2));
 VV = reshape(xgrid2D(:,2),m(1:2));
 
for n = 1:1:Ntstep
   
    QQ = reshape(Qprojxy{n},m(1:2));
    
    QQO = reshape(Qobsprojxy{n},m(1:2));

    hold on
    mesh(XX,VV,QQ)
    mesh(XX,VV,QQO)
end





