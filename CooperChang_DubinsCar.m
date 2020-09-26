clear
close

%velocity and yawrate
u = [5,1];

%time horizon
T = 2;
ts = 0.1;
Tvec = 0:ts:T;
Nsteps = length(Tvec)-1;

%initial condition
mu0 = [0;0;0]';
S0 = 0.01*eye(3);

N = [100,100,100];
grid_lower_bounds = [-5 -5 -pi/4];
grid_upper_bounds = [15 15 pi];
h = (grid_upper_bounds-grid_lower_bounds)./N;
xGrid = sub_to_x_val(ind2subnd(N,1:prod(N)),h,grid_lower_bounds);

%shift to make mu0 a grid point
[~,closestI] = min(sum((xGrid-mu0).^2,2));
xGrid = xGrid+(mu0-xGrid(closestI,:));
grid_lower_bounds = grid_lower_bounds+(mu0-xGrid(closestI,:));
grid_upper_bounds = grid_upper_bounds+(mu0-xGrid(closestI,:));

%% plot initial density
xGrid2D = sub_to_x_val(ind2subnd(N(1:2),1:prod(N(1:2))),h(1:2),grid_lower_bounds(1:2));
XX = reshape(xGrid2D(:,1),[N(1) N(2)]);
YY = reshape(xGrid2D(:,2),[N(1) N(2)]);

Q0 = mvnpdf(xGrid,mu0,S0);
Q02D = project_onto_dimension(Q0,[1 2],N,grid_lower_bounds,grid_upper_bounds);




%% define function handles for flux terms
x = sym('x',[1 3],'real');
f = [u(1)*cos(x(3));u(1)*sin(x(3));u(2)];
g = [0.1*cos(x(3)),0;0.1*sin(x(3)),0;0,0.01];
%g = zeros(3);
alpha = sym(zeros(3));
for i = 1:3
    for j = 1:3
        for k = 1:size(g,2)
            alpha(i,j) = alpha(i,j) + g(i,k)*g(j,k);
        end
    end
end

B = cell(3,1);
C = cell(3,1);

for i = 1:3
    Btemp = -f(i);
    Ctemp = sym(0);
    for j = 1:3
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
AIplus = cell(3,1);
AIminus = cell(3,1);

Ilist = (1:prod(N))';

subGrid = ind2subnd(N,1:prod(N));


idx_mat = sub2indnd([prod(N) prod(N)],[Ilist,Ilist]);
tic
% AI = 3*ones(prod(N),1);
% IdxJ = (1:prod(N))';
% IdxI = (1:prod(N))';

A = 3*speye(prod(N));
for j = 1:3
    
    %compute grid points +/- half cell in j direction
    idxstep = zeros(size(xGrid));
    idxstep(:,j) = 1;
    xGridplus = xGrid + h(j)/2*idxstep;
    xGridminus = xGrid - h(j)/2*idxstep;
    
    
    %     idxplus = sub2indnd(N,min(subGrid+idxstep,repmat(N,[prod(N) 1])));
    %     idxminus = sub2indnd(N,max(subGrid-idxstep,repmat(ones(1,3),[prod(N) 1])));
    Lp = subGrid(:,j)<N(j);
    Lm = subGrid(:,j)>1;
    
    idxplus = sub2indnd(N,subGrid(Lp,:)+idxstep(Lp,:));
    idxminus = sub2indnd(N,subGrid(Lm,:)-idxstep(Lm,:));
    
    
    Bplus = B{j}(xGridplus);
    Bminus = B{j}(xGridminus);
    Cplus = C{j}(xGridplus);
    Cminus = C{j}(xGridminus);
    
    wplus = h(j)*Bplus./(Cplus+1e-12);
    deltaplus = 1./wplus - 1./(exp(wplus)-1);
    
    deltaplus(wplus==0.0) =  0.5;
    
    wminus = h(j)*Bminus./(Cminus+1e-12);
    deltaminus = 1./wminus - 1./(exp(wminus)-1);
    
    deltaminus(wplus==0.0) =  0.5;
    
    %      AIplus{j} = -2*ts/h(j)*((1-deltaplus).*Bplus+1/h(j)*Cplus);
    %      AIminus{j} = -2*ts/h(j)*(-deltaminus.*Bminus+1/h(j)*Cminus);
    %
    
    idxplus_mat =  sub2indnd([prod(N) prod(N)],[Ilist(Lp),idxplus]);
    A(idxplus_mat) = A(idxplus_mat) -2*ts/h(j)*((1-deltaplus(Lp)).*Bplus(Lp)+1/h(j)*Cplus(Lp));
    
    
    idxminus_mat =  sub2indnd([prod(N) prod(N)],[Ilist(Lm),idxminus]);
    A(idxminus_mat) =A(idxminus_mat) -2*ts/h(j)*(-deltaminus(Lm).*Bminus(Lm)+1/h(j)*Cminus(Lm));
    
    
    A(idx_mat) = A(idx_mat) + 2*ts/h(j)*(-deltaplus.*Bplus+1/h(j)*Cplus) + 2*ts/h(j)*((1-deltaminus).*Bminus+1/h(j)*Cminus);
    if any(any(isnan(A)) )
        error('NaN Elements in A matrix')
    end
end
disp(['Built A matrix: ',num2str(toc),' s'])
% AIplus = cat(1,AIplus{:});
% idxplus = cat(1,idxplus{:});
% AIminus = cat(1,AIminus{:});
% idxminus = cat(1,idxminus{:});
%
% A = sparse(repmat( (1:prod(N))',[7 1]),[(1:prod(N))';idxplus;idxminus],[AI;AIplus;AIminus],prod(N),prod(N));
figure(1)
cla
spy(A)
pause
%% step forward
Q = NaN(prod(N),Nsteps+1);
Q(:,1) = Q0;

rt = NaN(1,Nsteps);
figure(1)
cla
mesh(XX,YY,reshape(Q02D,[N(1) N(2)]))
for k = 1:Nsteps
    b = 4*Q(:,k)-Q(:,max(1,k-1));
    
    tic
    %Q(:,k+1) = (A\b);
    Q(:,k+1) = pcg(A,b);
    rt(k) = toc;
    disp(['Iter ',num2str(k),'/',num2str(Nsteps),'; ',num2str(rt(k)),' s'])
    
    
    Q2D = project_onto_dimension(Q(:,k+1),[1 2],N,grid_lower_bounds,grid_upper_bounds);
    
    pause(0.1)
    mesh(XX,YY,reshape(Q2D,[N(1) N(2)]))
    
end










%% monte carlo
Nsamples = 1e3;
nw = size(g,2);
X = cell(1,Nsamples);

w = sym('w',[nw 1],'real');
dxdt = matlabFunction(f+g*w,'vars',{x,w});

tic
for j = 1:Nsamples
    X{j} = NaN(3,Nsteps+1);
    X{j}(:,1) = mu0 + randn([1 3])*chol(S0);
    for k = 1:Nsteps
        wtemp = randn([nw 1]);
        X{j}(:,k+1) = X{j}(:,k) + ts*dxdt( X{j}(:,k)',wtemp);
    end
end
% Xmean = mean(X,3);
%
% P2full = zeros(3,3,Nsteps+1);
% P2diag = zeros(3,Nsteps+1);
% P3diag = zeros(3,Nsteps+1);
% P4diag = zeros(3,Nsteps+1);
%
% for k = 1:Nsteps+1
%     for j = 1:Nsamples
%        P2full(:,:,k) =  P2full(:,:,k) + 1/Nsamples * (X(:,k,j)-Xmean(:,k))*(X(:,k,j)-Xmean(:,k))';
%        P3diag(:,k) = P3diag(:,k) + 1/Nsamples * (X(:,k,j)-Xmean(:,k)).^3;
%        P4diag(:,k) = P4diag(:,k) + 1/Nsamples * (X(:,k,j)-Xmean(:,k)).^4;
%
%     end
%     P2diag(:,k) = diag(P2full(:,:,k));
% end
% toc

Xcat = cat(2,X{:});
figure(1)
hold on
plot3(Xcat(1,:),Xcat(2,:),zeros(size(Xcat(1,:))),'k.')
