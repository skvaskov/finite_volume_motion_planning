clear
close all

% example for a longitudinal vehicle model. Dynamics are 2D
% d/dt [x;v] = [v;-2*atan(v-k)] where k is a commanded velocity
% the average initial condition for the vehicle is [0,15]
% the average initial condition for the obstacle is [6,10], the obstacle
% accelerates to 13 m/s

%Author: Sean Vasov
%Date: May 17, 2020

%-------------------------------------------------------------------------%

%system dimension
d = 2;

%number of cells on each dimension
m = [100, 25];

%grid lower and upper bounds
grid_lower_bounds = [-10 8];
grid_upper_bounds = [90 22];

%cell width in each dimension
h = (grid_upper_bounds-grid_lower_bounds)./m;

%dynamics
f =@(z,k) [z(2);-2*atan(z(2)-k)];

%number of particles to use for initial condition and monte carlo
N = 1e3;

%time horizon
T = 4;

%CFL number (between 0 and 1) lower means smaller timestep
cfl_num = 0.8;

%velocity command to test
k_test = 19;

%% compute average velocities
disp('computing average velocities for cell')

fbar = compute_average_velocity_functions(f,1,h);

tic
ubar = compute_average_velocities(k_test,fbar,m,grid_lower_bounds,grid_upper_bounds);
toc

dt = min(cfl_num*h(1)./max(abs(ubar{1})),cfl_num*h(2)./max(abs(ubar{2})));

Ntstep = floor(T/dt);


%% compute initial condition

mean_mc = NaN(2,Ntstep);
std_mc =  NaN(2,Ntstep);

mu = [0,15];
s =  [0.1,0.2];

p = cell(1,Ntstep);

%uniform
diffx = sqrt(s(1)*12);
diffv = sqrt(s(2)*12);

p{1} = [randRange(mu(1)-diffx/2,mu(1)+diffx/2,[],[],1,N);...
        randRange(mu(2)-diffv/2,mu(2)+diffv/2,[],[],1,N)];


mean_mc(:,1) = [mean(p{1}(1,:)); mean(p{1}(2,:))];
std_mc(:,1) = [std(p{1}(1,:));std(p{1}(2,:))];

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
tic


%% forward advect
Q = fvm_advection(Q0,dt,h,ubar,UP,DQI,DQi,DF,T);

disp(['initial volume: ',num2str(sum(sum(Q{1}))*prod(h))])
disp(['final volume: '  ,num2str(sum(sum(Q{end}))*prod(h))])


%% calculate probability of obstacle and also try our dynamics with gaussian
mu_obs = NaN(floor(1/dt),2);
mu_ego = NaN(floor(1/dt),2);

mu_obs(1,:) = [6,10];
mu_ego(1,:) = mu;

S_ego =cell(floor(1/dt),1);
S_obs = cell(floor(1/dt),1);
std_ego = NaN(floor(1/dt),2);

S_ego{1} = diag(s);
S_obs{1} = diag([0.5,0.5]);

std_obs(1,:) = sqrt([0.5,0.5]);
std_ego(1,:) = sqrt([S_ego{1}(1,1),S_ego{1}(2,2)]);

%% generate obstacle densities

zsym = sym('z',[2 1],'real');
g = matlabFunction(eye(2)*zsym+dt*f(zsym,k_test),'vars',{zsym});

G = matlabFunction(jacobian(g(zsym),zsym),'vars',{zsym});

slow_flag = false;

disp('forward propogating obstacle and ego vehicle gaussians')
tic
for n = 1:Ntstep-1
    
    if mu_obs(n,2)<13
       acc = 3.5;
    else
        acc = -3;
        slow_flag = true;
    end
    
    if slow_flag && mu_obs(n,2) > 12
        acc = -2;
    end
%     
    
    %obstacle
    mu_obs(n+1,:) = mu_obs(n,:) + dt*[mu_obs(n,2),acc];
    Go = eye(2)+dt*[0 1;0 0];
    
    S_obs{n+1} = Go*S_obs{n}*Go';
    std_obs(n+1,:) = sqrt([S_obs{n+1}(1,1),S_obs{n+1}(2,2)]);
    
    %us
    
    mu_ego(n+1,:) = g(mu_ego(n,:)')';
    
    Gn = G(mu_ego(n,:)');
    
    S_ego{n+1} = Gn*S_ego{n}*Gn';
    
    std_ego(n+1,:) = sqrt([S_ego{n+1}(1,1),S_ego{n+1}(2,2)]);
end
toc


%% pplot densities 
xgrid = get_grid_points(m,grid_lower_bounds,grid_upper_bounds);
XX = reshape(xgrid(:,1),m);
VV = reshape(xgrid(:,2),m);

for n = 1:10:Ntstep
    
        QQ = reshape(Q{n},m);

        hold on
        mesh(XX,VV,QQ)
        
         mesh(XX,VV,reshape(mvnpdf([XX(:) ,VV(:)],mu_obs(n,:),S_obs{n}),m));

end
% 
 xlabel('x')
 ylabel('v')

 title('Densities (obstacle on bottom)')


%% calculating and plotting statistics
figure(2)

tvec = (0:(Ntstep-1))*dt;

[meanQ,stdQ] = compute_mean_and_stddev(Q,m,grid_lower_bounds,grid_upper_bounds);

subplot(2,2,1)
plot(tvec,meanQ(:,1))
hold on

title('mean')
xlabel('time')
ylabel('x')

subplot(2,2,2)
plot(tvec,meanQ(:,2))
hold on

title('mean')
xlabel('time')
ylabel('v')

subplot(2,2,3)
plot(tvec,stdQ(:,1))
hold on

title('std dev')
xlabel('time')
ylabel('x')


subplot(2,2,4)
plot(tvec,stdQ(:,2))
hold on

title('std dev')
xlabel('time')
ylabel('v')



%% particle
disp('running monte carlo')

for i = 1:Ntstep-1
    pttt = NaN(2,N);
    pit = p{i};
    parfor j = 1:N
        [~,ptemp] = ode45(@(t,z) f(z,k_test),[0 dt],pit(:,j));
        pttt(:,j) = ptemp(end,:)'
    end
    p{i+1} = pttt;
    
mean_mc(:,i+1) = [mean(p{i+1}(1,:)); mean(p{i+1}(2,:))];
std_mc(:,i+1) =  [std(p{i+1}(1,:)) ;  std(p{i+1}(2,:))];
    
end



figure(2)
subplot(2,2,1)
plot(tvec,mean_mc(1,:),'o')
plot(tvec,mu_ego(:,1))
legend({'finite volume','monte carlo','gaussian'})

subplot(2,2,2)
plot(tvec, mean_mc(2,:),'o')
plot(tvec, mu_ego(:,2))


subplot(2,2,3)
plot(tvec,std_mc(1,:),'o')
plot(tvec,std_ego(:,1))

subplot(2,2,4)
plot(tvec,std_mc(2,:),'o')
plot(tvec,std_ego(:,2))

pause(0.1)


%% caclulate probability of collision at each timestep
disp('calculating probability and energy of collision ')
Pcol = zeros(1,length(Q));
Ecol = zeros(1,length(Q));
Pcolgaus = zeros(1,length(Q));
Ecolgaus = zeros(1,length(Q));
Pcolmc = zeros(1,length(Q));
Ecolmc = zeros(1,length(Q));

ftprint = 4;
mo = 1500;
mc = 1500;

ni = NaN(1,Ntstep);

for n = 1:Ntstep
    
    %find cells that are within footprint in this case it is just
    %neighboring cells
    robs = repmat(mu_obs(n,:),[N 1])+randn(N,2)*chol(S_obs{n});
    is_to_check = find(Q{n});
    ni(n) = length(is_to_check);
    
    for i = is_to_check'
           
        
           sub_i = ind2subnd([m m],i);
           
           xl = xgrid(i,1)+[0 h(1)]+[-ftprint(1) ftprint(1)];
           vl = xgrid(i,2)+[0 h(2)];
         
           
           obs_moments = multivdtmom_unormalized([0 2],[xl(1);0],[xl(2);30],mu_obs(n,:),S_obs{n});
   
           Pcol(n) = Pcol(n)+Q{n}(i)*prod(h)*obs_moments(1);

           Ecol(n) = Ecol(n)+Q{n}(i)*(prod(h)*obs_moments(3)-...
                        2*(vl(2)^2-vl(1)^2)/2*h(1)*obs_moments(2)+(vl(2)^3-vl(1)^3)/3*h(1)*obs_moments(1));     
        
    end

    
    for i = 1:N
        
        xl = p{n}(1,i)+[-ftprint(1) ftprint(1)];
        vl = p{n}(2,i);
        
        obs_moments = multivdtmom_unormalized([0 2],[xl(1);0],[xl(2);50],mu_obs(n,:),S_obs{n});
        
        Pcolmc(n) = Pcolmc(n)+obs_moments(1);
        
        Ecolmc(n) = Ecolmc(n)+(obs_moments(3)-2*vl*obs_moments(2)+vl^2*obs_moments(1));
        
    end
    
    
    rego = repmat(mu_ego(n,:),[N 1])+randn(N,2)*chol(S_ego{n});
    
    Lcol = abs(robs(:,1)-rego(:,1))<ftprint;
    
    Pcolgaus(n) = sum(Lcol)/N;
    Ecolgaus(n) = sum((robs(Lcol,2)-rego(Lcol,2)).^2)/N;
    
  
end

disp(['max number of cells used: ',num2str(max(ni))])
disp(['mean number of cells used: ',num2str(mean(ni))])

energy_mult = mc^2*mo/2/(mc+mo)^2;
Ecol =  energy_mult*Ecol;
Ecolmc = energy_mult*Ecolmc/N;
Ecolgaus = energy_mult*Ecolgaus;
Pcolmc = Pcolmc/N;

figure(1)
subplot(2,2,1)
plot((0:length(Q)-1)*dt*T,Pcol)
hold on
plot((0:length(Q)-1)*dt*T,Pcolmc,'o')
plot((0:length(Q)-1)*dt*T,Pcolgaus)
ylabel('P_{col}')
xlabel('time')
grid on
subplot(2,2,2)
plot((0:length(Q)-1)*dt*T,Pcol-Pcolmc)
hold on
plot(0,0)
plot((0:length(Q)-1)*dt*T,Pcolgaus-Pcolmc)
xlabel('time')
ylabel('Error in P_{col}')
grid on
subplot(2,2,3)
plot((0:length(Q)-1)*dt*T,Ecol)
hold on
plot((0:length(Q)-1)*dt*T,Ecolmc,'o')
plot((0:length(Q)-1)*dt*T,Ecolgaus)

xlabel('time')
ylabel('E_{col}')
legend({'finite volume','montecarlo','gaussian'})
grid on
subplot(2,2,4)
plot((0:length(Q)-1)*dt*T,Ecol-Ecolmc)
hold on
plot(0,0)
plot((0:length(Q)-1)*dt*T,Ecolgaus-Ecolmc)
grid on
xlabel('time')
ylabel('Error in E_{col}')

