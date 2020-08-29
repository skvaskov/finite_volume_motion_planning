clear
close

%velocity and yawrate
u = 10;

%time horizon
T = 2;
ts = 0.01;
Tvec = 0:ts:T;
Nsteps = length(Tvec)-1;

%initial condition
mu0 = 4;
S0 = 1;
d = 1;

N = 100;
grid_lower_bounds = 0;
grid_upper_bounds = 20;
h = (grid_upper_bounds-grid_lower_bounds)./N;
xGrid = sub_to_x_val(ind2subnd(N,1:prod(N)),h,grid_lower_bounds);

%shift to make mu0 a grid point
[~,closestI] = min(sum((xGrid-mu0).^2,2));
xGrid = xGrid+(mu0-xGrid(closestI,:));
grid_lower_bounds = grid_lower_bounds+(mu0-xGrid(closestI,:));
grid_upper_bounds = grid_upper_bounds+(mu0-xGrid(closestI,:));

%% plot initial density


Q0 = mvnpdf(xGrid,mu0,S0);

plot(xGrid,Q0);


%% define function handles for flux terms 
x = sym('x',[1 1],'real');
f = -(x(1)-u(1));
g =   0.1;
alpha = sym(zeros(3));
for i = 1:d
    for j = 1:d
        for k = 1:size(g,2)
            alpha(i,j) = alpha(i,j) + g(i,k)*g(j,k); 
        end
    end
end

B = cell(d,1);
C = cell(d,1);

for i = 1:d
    Btemp = -f(i);
    Ctemp = sym(0);
    for j = 1:d
       Btemp = Btemp + 1/2*(diff(alpha(i,j),x(j))); 
       Ctemp = Ctemp + 1/2*alpha(i,j);
    end
    if isempty(symvar(Btemp))
        B{i} = @(xin) double(Btemp)*ones(size(xin,1),1);
    else
        B{i} = matlabFunction(Btemp,'vars',{x});
    end
    if isempty(symvar(Ctemp))
        C{i} = @(xin) double(Ctemp)*ones(size(xin,1),1);
    else
        C{i} = matlabFunction(Ctemp,'vars',{x});
    end
end

%% flux matrix



indI =  sub2ind([N N], (1:N)', (1:N)');
idxplus =  sub2ind([N N],(1:N)',[2:N,1]');
idxminus = sub2ind([N N],(1:N)',[N,1:N-1]');
    
A = 3*speye(prod(N));
for j = 1:d
    
    %compute grid points +/- half cell in j direction
    idxstep = zeros(size(xGrid));
    idxstep(:,j) = 1;
    xGridplus = xGrid + h(j)/2*idxstep;
    xGridminus = xGrid - h(j)/2*idxstep;
    

    
    Bplus = B{j}(xGridplus);
    Bminus = B{j}(xGridminus);
    Cplus = C{j}(xGridplus);
    Cminus = C{j}(xGridminus);
    
    wplus = h(j)*Bplus./Cplus;
    deltaplus = 1./wplus - 1./(exp(wplus)-1);
    
    wminus = h(j)*Bminus./Cminus;
    deltaminus = 1./wminus - 1./(exp(wminus)-1);
    
    
    A(idxplus)= A(idxplus)-2*ts/h(j)*((1-deltaplus).*Bplus+1/h(j)*Cplus);
    A(idxminus) = A(idxminus)-2*ts/h(j)*(-deltaminus.*Bminus+1/h(j)*Cminus);
     
    A(indI) = A(indI) + 2*ts/h(j)*(-deltaplus.*Bplus+1/h(j)*Cplus) + 2*ts/h(j)*((1-deltaminus).*Bminus+1/h(j)*Cminus);
end



%% step forward
Q = NaN(prod(N),Nsteps+1);
Q(:,1) = Q0;

rt = NaN(1,Nsteps);

for k = 1:Nsteps
    b = 4*Q(:,k)-Q(:,max(1,k-1));
    
    tic
    Q(:,k+1) = (A\b);
    rt(k) = toc;
    disp(['Iter ',num2str(k),'/',num2str(Nsteps),'; ',num2str(rt(k)),' s'])
    
end










%% monte carlo
Nsamples = 1e3;
X = cell(1,Nsamples);

w = sym('w',[size(g,2) 1],'real');
dxdt = matlabFunction(f+g*w,'vars',{x,w});


tic
for j = 1:Nsamples
    X{j} = NaN(d,Nsteps+1);
    X{j}(:,1) = mu0 + randn([1 d])*chol(S0);
 
    for k = 1:Nsteps
        wtemp = randn([size(g,2) 1]);
        X{j}(:,k+1) = X{j}(:,k) + ts*dxdt( X{j}(:,k)',wtemp);  
        
    end
end

Xcat = cat(1,X{:});

%% make gif
fn = figure(1);
cla
ylim([0 max(max(Q))*1.2])
xlim([grid_lower_bounds,grid_upper_bounds])

axis tight manual % this ensures that getframe() returns a consistent size
filename = 'cooperchang1D.gif';

histogram(Xcat(:,1),'normalization','pdf');
plot(xGrid,Q(:,1),'Color','k','LineWidth',1.0);

for k = 1:Nsteps
    
      % Capture the plot as an image 
      frame = getframe(fn); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      
      % Write to the GIF File 
      if k == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',ts); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',ts); 
      end 
      cla
      histogram(Xcat(:,k+1),'normalization','pdf');
      hold on
      plot(xGrid,Q(:,k+1),'Color','k','LineWidth',1.0);
end

