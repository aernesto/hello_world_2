%lkjhiu ;lkj kj;kij;
%adlkfj;alkdjfa c
% compbiasvITI.m
%
% computes the bias as a function of the difference between the 
%

h = 0.1;    % threshold
tau = 100;   % time constant of plasticity 
b = 0.01;      % rate of plasiticty change
p0 = 3;     % plasticity target value
s = 0.005;

N = 180;% 2000;    % number of grid points
Nsim = 1;  % number of sims per point
ind = [1:N];
dt = 0.1;  % time step
T = 1300;    % end time
dx = 2*pi/N;    % space step
x = linspace(-pi,pi-dx,N)';  % space grid
nt = round(T/dt)+1; % number of time points
timey = linspace(0,T,nt);   % time vector
cmu = zeros(nt,1);
I=nan(N,nt); 
ITIs = linspace(100,500,21);
bias = zeros(21,1);
bhi = bias; blo = bias; bs = zeros(Nsim,1);
tvec = linspace(0,500,21);
U=nan(N,nt);
U(:,1)=0;
P=nan(N,nt);
P(:,1)=1;
tmod=zeros(1,nt);

x1 = pi/2;
x2=0;

Delay1=300;
Delay2=500;

evidenceMod=1%2;
DecayMod=10;

for k=1%:21,
    ITI=ITIs(k);
for l=1:Nsim, %U = zeros(N,1); P = ones(N,1);
for j = 1:nt-1,
    I(:,j)=0; tmod(j) = j*dt;
	if 10<tmod(j) & tmod(j)<60, I(:,j) = evidenceMod*cos(x-x1); end % was 2*cos(x-x1)
%          if Delay1+10<tmod(j) & tmod(j)<Delay1+60, I(:,j) = -5; end %changed from 5
    if Delay1+10+ITI<tmod(j) & tmod(j)<Delay1+40+ITI, I(:,j) =  2*cos(x-x1); end
    
    fu = heaviside(U(:,j)-h);
    f1c = dx*P(:,j)'.*cos(x')*fu;    f1s = dx*P(:,j)'.*sin(x')*fu;
    nos = randn*cos(x)+randn*sin(x);
  
    U(:,j+1) = U(:,j) + dt*(DecayMod*-U(:,j)+f1c*cos(x)+f1s*sin(x)+I(:,j))+s*sqrt(dt)*nos;
    P(:,j+1) = P(:,j) + dt*(1-P(:,j)+b*fu.*(p0-P(:,j)))/tau;
    if j*dt==Delay1+60+ITI+Delay2 
        [junk,mi] = max(U(:,j+1)); bs(l) = x(mi);
        break
    end
  
end
end
    bias(k)=sum(bs)/Nsim; bhi(k)=bias(k)+std(bs); blo(k)=bias(k)-std(bs);
end

 figure
 colormap(hot)
 imagesc([0,T*10],[-180,180],U);
 hold on
 
 figure
 [~,TEST]=max(U);
 plot(x(TEST)*180/pi);
 bs*180/pi
% degbias = 180*bias'/pi; deghi=180*bhi'/pi; deglo=180*blo'/pi; degvec=ITIs;
% % hold on, plot(x1vec,bhi,'r','linewidth',2);
% % hold on, plot(x1vec,blo,'r','linewidth',2);
% clear stdpatch
% vertx=[degvec'; degvec(end:-1:1)']; vertx(:,2)=[deghi'; deglo(end:-1:1)'];
% stdpatch.Vertices=vertx;
% stdpatch.Faces=[1:42];
% stdpatch.FaceColor=[250 128 114]/255;
% stdpatch.EdgeColor=[250 128 114]/255;
% figure(1), hold on, patch(stdpatch);
% plot(degvec,degbias,'r-o','markersize',12,'linewidth',4);
% axis([100 500 0 15])
% set(gca,'xtick',[0:100:500]);
% set(gca,'fontsize',30);

% figure(2), hold on, plot(degvec,deghi-degbias,'r-o','markersize',12,'linewidth',4);
% set(gca,'xtick',[-180:90:180]);
% set(gca,'fontsize',30);