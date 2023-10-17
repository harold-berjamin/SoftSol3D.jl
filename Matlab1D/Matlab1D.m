% Chocka.m implements a (MUSCL)-Osher finite volume scheme for plane
% shear-wave propagation in viscoelastic soft solids. The material is
% assumed incompressible, with nonlinear behaviour and Fung-Simo QLV rheology.
% By Harold Berjamin, NUI Galway, 2022.

% clear;

%% settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mechanical parameters from Tripathi et al. (2019)
rho = 1e3; % 1.0e3
mu  = 2.684e3; % 2.684e3
c0  = sqrt(mu/rho);  % linear shear-wave speed
g1  = 0.0434; % 0.0434
g2  = 0.0466; % 0.0466
g3  = 0.2213; % 0.2213
w1 = 2*pi * 1e1; % 2*pi * 1e1
w2 = 2*pi * 1e2; % 2*pi * 1e2
w3 = 2*pi * 1e3; % 2*pi * 1e3
bet = 4.4;   % 4.4
law = 'pol'; % constitutive law, among lin, pol

Om = 82.35; % frequency
src = @(t) 2e0 * (sin(Om*t) - 0.5*sin(2*Om*t)) .* ((t>=0) - (t>=2*pi/Om)); % forcing signal

% numerics
Nx = 1e3; % number of cells
Co = 0.95; % Courant number
x  = linspace(-0.2, 0.2, Nx+1); % abscissas
tend = 0.08; % final time 8e-2 to 24e-2
slim = @(a,b) mm(a,b); % limiter, among zero, lincomb, mm, MC - cf. end

% graphics
visu = 5e2; % plot update every 'visu' iterations

%% initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grid
dx = x(2) - x(1); % cell size
x  = [x(1)-2*dx, x(1)-dx, x, x(end)+dx, x(end)+2*dx]; % abscissas
is = find(abs(x)<dx/2); % x=0

% constitutive law, speed of sound
[W, sig, S12, S22, c] = getConstitutive(law, mu, rho, bet);

% fields
t = 0;
n = 0;
Q = [0*x; 0*x; 0*x; 0*x; 0*x; 0*x; 0*x; 0*x];
Qnew = Q;
fp05 = Q;
fm05 = Q;
qEvol = Q;

% CFL condition
velo = 0*x;
for i=1:Nx+5
    velo(i) = c(Q(:,i));
end
cmax = max(velo);
dt = Co * dx / cmax;

% Gaussian quadrature params
s1 = 0.1*(5-sqrt(15));
s2 = 0.5;
s3 = 0.1*(5+sqrt(15));
W1 = 5/18;
W2 = 8/18;
W3 = 5/18;

% graphics
figure(1);
clf;

subplot(1,2,1);
hnum1 = plot(x, Q(1,:), 'b.-', 'DisplayName','num'); hold on
xlim([x(1) x(end)]);
% ylim([-1.5, 1.5] * 0.5*Ampl/c0^2);
xline(x(3),'k:');
xline(x(end-2),'k:');
yline(0, 'k-'); 
xlabel('x (m)');
ylabel('\gamma');

subplot(1,2,2);
hnum2 = plot(x, Q(2,:), 'b.-', 'DisplayName','num'); hold on
xlim([x(1) x(end)]);
% ylim([-1.5, 1.5] * 0.5*Ampl/c0);
xline(x(3),'k:');
xline(x(end-2),'k:');
yline(0, 'k-'); 
xlabel('x (m)');
ylabel('v (m/s)');

ht = sgtitle(strcat('t = ',num2str(t)));

% system
F  = @(q) [-q(2); -sig(q)/rho; 0; 0; 0; 0; 0; 0]; % flux function
A  = @(q) [0, -1, 0, 0, 0, 0, 0, 0; ...
    -c(q)^2, 0, 1/rho, q(1)/rho, 1/rho, q(1)/rho, 1/rho, q(1)/rho; ...
    0, 0, 0, 0, 0, 0, 0, 0; ...
    0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0]; % df/dq
Aa = @(q) [c(q), 0, -1/rho/c(q), -q(1)/rho/c(q), -1/rho/c(q), -q(1)/rho/c(q), -1/rho/c(q), -q(1)/rho/c(q);
    0, c(q), 0, 0, 0, 0, 0, 0; ...
    0, 0, 0, 0, 0, 0, 0, 0; ...
    0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0]; % |df/dq|

%% loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

compt = tic;

while (t + dt < tend)
    % relaxation step -----------------------------------------------------
    for i=3:Nx+3
        S12i = S12(Q(:,i));
        S22i = S22(Q(:,i));
        Q(3,i) = exp(-0.5*dt*w1)*Q(3,i) + (1-exp(-0.5*dt*w1))*g1*S12i;
        Q(4,i) = exp(-0.5*dt*w1)*Q(4,i) + (1-exp(-0.5*dt*w1))*g1*S22i;
        Q(5,i) = exp(-0.5*dt*w2)*Q(5,i) + (1-exp(-0.5*dt*w2))*g2*S12i;
        Q(6,i) = exp(-0.5*dt*w2)*Q(6,i) + (1-exp(-0.5*dt*w2))*g2*S22i;
        Q(7,i) = exp(-0.5*dt*w3)*Q(7,i) + (1-exp(-0.5*dt*w3))*g3*S12i;
        Q(8,i) = exp(-0.5*dt*w3)*Q(8,i) + (1-exp(-0.5*dt*w3))*g3*S22i;
    end
    
    % propagation step ----------------------------------------------------
    
    % shifts
    Qp1 = circshift(Q, [0 -1]);
    Qm1 = circshift(Q, [0  1]);
    
    % MUSCL-Osher flux (2nd order MUSCL, gauss-legendre 3)
    % evolved, extroplated values
    Delt = slim(Q-Qm1, Qp1-Q);
    qR = Q + 0.5*Delt;
    qL = Q - 0.5*Delt;
    for i=1:Nx+5
        qEvol(:,i) = 0.5*dt/dx*(F(qL(:,i))-F(qR(:,i)));
    end
    qMp05 = qR + qEvol;
    qPp05 = circshift(qL+qEvol, [0 -1]);
    qDiff = qPp05 - qMp05;
    % Osher flux computation
    for i=2:Nx+3
        Aa1 = Aa(qMp05(:,i) + s1*qDiff(:,i));
        Aa2 = Aa(qMp05(:,i) + s2*qDiff(:,i));
        Aa3 = Aa(qMp05(:,i) + s3*qDiff(:,i));
        Aamoy = W1*Aa1 + W2*Aa2 + W3*Aa3;
        fp05(:,i) = 0.5*(F(qMp05(:,i))+F(qPp05(:,i))) - 0.5*Aamoy*qDiff(:,i);
    end
    
    % update (flux differences)
    fm05 = circshift(fp05, [0  1]);
    Qnew = Q - dt/dx*(fp05-fm05);

    % inflow/outflow BC
    Qnew(:,1)  = Qnew(:,3);
    Qnew(:,2)  = Qnew(:,3);
    Qnew(:,Nx+4) = Qnew(:,Nx+3);
    Qnew(:,Nx+5) = Qnew(:,Nx+3);
    
    % relaxation step -----------------------------------------------------
    for i=3:Nx+3
        S12i = S12(Qnew(:,i));
        S22i = S22(Qnew(:,i));
        Qnew(3,i) = exp(-0.5*dt*w1)*Qnew(3,i) + (1-exp(-0.5*dt*w1))*g1*S12i;
        Qnew(4,i) = exp(-0.5*dt*w1)*Qnew(4,i) + (1-exp(-0.5*dt*w1))*g1*S22i;
        Qnew(5,i) = exp(-0.5*dt*w2)*Qnew(5,i) + (1-exp(-0.5*dt*w2))*g2*S12i;
        Qnew(6,i) = exp(-0.5*dt*w2)*Qnew(6,i) + (1-exp(-0.5*dt*w2))*g2*S22i;
        Qnew(7,i) = exp(-0.5*dt*w3)*Qnew(7,i) + (1-exp(-0.5*dt*w3))*g3*S12i;
        Qnew(8,i) = exp(-0.5*dt*w3)*Qnew(8,i) + (1-exp(-0.5*dt*w3))*g3*S22i;
    end
    
    % src
    Qnew(2,is) = Qnew(2,is) + dt*src(t+dt)/dx;
    
    % update
    Q = Qnew;
    t = t + dt;
    n = n + 1;
    
    for i=3:Nx+3
        velo(i) = c(Q(:,i));
    end
    cmax = max(velo);
    dt = Co * dx / cmax;
    
    % graphics
    if rem(n,visu)==0
        set(hnum1, 'YData', Q(1,:));
        set(hnum2, 'YData', Q(2,:));
        set(ht, 'String', strcat('t = ',num2str(t)));
        drawnow;
    end
end

disp(strcat('Computational time (s):', num2str(toc(compt))));

% final plot
if visu>0
    % analytical solution (linear)
    % FFT analysis and synthesis
    % sampling
    T = 2*pi/Om;
    Ne = 5e2;
    fe = Ne/T;
    te = (0:Ne-1)/fe;
    se = src(te);
    % zero-padding
    se = [se, zeros(1,2^(nextpow2(10*Ne))-Ne)];
    Nez = length(se);
    % fft
    she = fftshift( fft(se, Nez)/fe );
    wez = 2*pi*fe * (-Nez/2:(Nez-1)/2)/Nez;
    cI = c0 * sqrt( 1 - g1*w1./(w1+1i*wez) - g2*w2./(w2+1i*wez) - g3*w3./(w3+1i*wez) );
    p = 0.5*she.*fe./cI;
    vth = linspace(1,Nx+5);
    for i=1:Nx+5
        % ifft
        shei = exp(1i*wez.*(t - abs(x(i)-x(is))./cI));
        shhe = ifft( ifftshift(p.*shei) );
        ghhe = ifft( ifftshift(p.*shei./cI) );
        vth(i) = real( shhe(1) );
        gth(i) = -sign(x(i)) * real( ghhe(1) );
    end
    % plot
    subplot(1,2,1);
    set(hnum1, 'YData', Q(1,:));
    hold on
    plot(x,-sign(x).*src(t-abs(x)/c0)/(2*c0^2),'k--');
    plot(x,gth,'k-');
    subplot(1,2,2);
    set(hnum2, 'YData', Q(2,:));
    hold on
    plot(x,src(t-abs(x)/c0)/(2*c0),'k--');
    plot(x,vth,'k-');
    set(ht, 'String', strcat('t = ',num2str(t)));
    drawnow;
end

%% functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% models
function [W, sig, S12, S22, c] = getConstitutive(law, mu, rho, bet)
    c0 = sqrt(mu/rho);
    if strcmp(law,'lin')
        W   = @(q) 0.5*mu*q(1)^2;
        sig = @(q) mu*q(1) - q(3) - q(5) - q(7);
        S12 = @(q) mu*q(1);
        S22 = @(q) 0*q(1);
        c = @(q) c0;
    elseif strcmp(law,'pol')
        W   = @(q) 0.5*mu*q(1)^2*(1+bet/3*q(1)^2);
        sig = @(q)  mu*(q(1)+2*bet/3*q(1)^3) - q(1)*(q(4)+q(6)+q(8)) - q(3) - q(5) - q(7);
        S12 = @(q)  mu*(q(1)+2*bet/3*q(1)^3)*(1+q(1)^2/3);
        S22 = @(q) -mu*(q(1)+2*bet/3*q(1)^3)*q(1)/3;
        c = @(q) sqrt( (mu*(1+2*bet*q(1)^2) - (q(4)+q(6)+q(8)))/rho );
    else
        disp('unknown law!'); return
    end
end

% limiters

function res = lincomb(a,b,k)
    res = 0.5*(1+k)*a + 0.5*(1-k)*b;
end

function res = zero(a,b)
    res = 0*a;
end

function res = mm(a,b)
    s = size(a);
    res = zeros(s);
    nl = s(1);
    for i=1:nl
        res(i,:) = 0.5*(sign(a(i,:))+sign(b(i,:))).*min([abs(a(i,:));abs(b(i,:))]);
    end
end

function res = MC(a,b)
    s = size(a);
    res = zeros(s);
    nl = s(1);
    for i=1:nl
        res(i,:) = 0.5*(sign(a(i,:))+sign(b(i,:))).*min([2*abs(a(i,:));0.5*abs(a(i,:)+b(i,:));2*abs(b(i,:))]);
    end
end