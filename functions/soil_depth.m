function [soil_depth] = soil_depth(Z,x_cellsize,endtime,DEM_name,plots_on)
% Non-linear hillslope evolution code for generating soil depths (from Josh
% Roering)
f = waitbar(0,'Computing Soil Depth','Name',['Soil Depth for ' DEM_name]);

% Predefined parameters
z=Z;
zrock=z-1;
zchan=zeros(size(z));
% dx=x_cellsize/3.28;
dx=x_cellsize;
E=0;
interval=1000;
waitbar_update_interval=3; % seconds

% Determine sizes of matrices and break code if they don't match
[m,n]=size(z);
[mzr,nzr]=size(zrock);
[mc,nc]=size(zchan);

x=repmat(0:dx:(n-1)*dx,m,1);
y=repmat((0:dx:(m-1)*dx)',1,n);


if m~=mzr && n~=nzr && mc~=mzr && nc~=nzr
    'z, zrock, and zchan matrices not equal'
    return
end

% Parameters Flags
depthK=0;
slopeprod=0;

% FireInt=fire frequency; burn=time that burn  parameters persist
fireint=9e9;
burn=0;

% Background transport parameters
backalpha=0.007;
backbeta=-2;
backSc=1.25;
backdt=1;

% Fire transport parameters
firealpha=0.045;
firebeta=-5;
fireSc=1.03;
firedt=0.25;

% Soil production parameters
Po=0.0003;
mu=3;
prps=2.0;

% Set K and Sc values if not using K=f(depth) model
if depthK~=1
    K(1:m,1:n)=0.005;
    Sc=1.25;
    dt=1;
end

% Additional parameters and initialization
criteria=(burn/2)/fireint;
criteria=0;
grac=0.025;
time=0;
cycle=1;
% disp('Simulating Hillslope Evolution: 000000 years');

zhill=(zchan<0.9);

% Add initial surfaces to output matrix
elev(1:m,1:n,1)=z;
soil(1:m,1:n,1)=z-zrock;
[dzdx,dzdy]=gradient(z,dx);
grado(1:m,1:n,1)=sqrt((dzdx.^2)+(dzdy.^2));
curv(1:m,1:n,1)=4.*del2(z,dx);
times(1)=0;

tocstore=0;
timestore=0;
tic;
while time <= (endtime+(interval/2))
    zstore=z;
    
    % Calculate transport parameters depending on fire cycle and soil depth
    if depthK==1
        frac=abs((time/fireint)-round(time/fireint));
        if frac < criteria
            K=(firealpha/firebeta)*(exp(firebeta.*(z-zrock))-1);
            Sc=fireSc;
            dt=firedt;
        else
            K=(backalpha/backbeta)*(exp(backbeta.*(z-zrock))-1);
            Sc=backSc;
            dt=backdt;
        end
    end
    
    % Calculate gradient
    [dzdx,dzdy]=gradient(z,dx);
    grad=sqrt((dzdx.^2)+(dzdy.^2));
    
    % Identifies points where gradients are near or exceed Sc
    gradless=grad<(Sc-(grac*Sc));
    gradmore=gradless<0.9;
    
    % Calculate cosine vector and lower bedrock surface along interior nodes
    if slopeprod==1
        cosin=cos(atan(grad));
        zrock=zhill.*(zrock-(((dt*Po)./cosin).*exp(-1*mu.*(z-zrock).*cosin)));
    else
        zrock=zhill.*(zrock-((dt*Po)*exp(-1*mu.*(z-zrock))));
    end
    
    % Erode hillslope nodes according to nonlinear transport
    den=1-((grad./Sc).^2);
    [d2zdx,d2zdxy]=gradient(dzdx,dx);
    [d2zdyx,d2zdy]=gradient(dzdy,dx);
    lap=4.*del2(z,dx);
    num2=2*(((dzdx.^2).*d2zdx)+((dzdy.^2).*d2zdy)+(2.*dzdy.*dzdx.*d2zdxy));
    z=z+((gradless.*zhill.*(dt/prps).*K).*((lap./den)+(num2./(Sc*Sc.*den.*den))));
    cellrestore=find(isnan(z)==1);
    z(cellrestore)=zstore(cellrestore);
    
    
    
    % If erosion exceeds soil thickness, set surface elevation equal to bedrock
    z=(zhill.*gradless.*((zrock.*(z<=zrock))+(z.*(z>zrock))))+zhill.*gradmore.*zrock+(zchan.*z);
    cellrestore=find(isnan(z)==1);
    z(cellrestore)=zstore(cellrestore);
    
    
    % Output interval done here
%     if time >= (interval*cycle)
%         cycle=cycle+1;
%         elev(1:m,1:n,cycle)=z;
%         soil(1:m,1:n,cycle)=(0./zhill)+(z-zrock);
%         [dzdx2,dzdy2]=gradient(z,2);
%         grado(1:m,1:n,cycle)=(0./zhill)+sqrt((dzdx2.^2)+(dzdy2.^2));
%         curv(1:m,1:n,cycle)=(0./zhill)+(4.*del2(z,dx));
%         times(cycle)=time;
%         output=[num2str(time)];
%         disp(output);
%     end
    
    % Drop the boundaries by one increment and tally time clock
    z=(zchan.*(z-(E*dt)))+zhill.*z;
    cellrestore=find(isnan(z)==1);
    z(cellrestore)=zstore(cellrestore);
    
    time=time+dt;
    
%     if rem(time,5)==0 && toc/time*(endtime-time)/60>=0
    if toc-tocstore>=waitbar_update_interval
        if (toc-tocstore)/(time-timestore)*((endtime+(interval/2))-time)/60 >= 60
        waitbar(time/endtime,f,['Computing Soil Depth (Time Remaining: ' ...
            num2str((toc-tocstore)/(time-timestore)*((endtime+(interval/2))-time)/3600,'%.2f') ' hours)']);
        else
           waitbar(time/endtime,f,['Computing Soil Depth (Time Remaining: ' ...
            num2str((toc-tocstore)/(time-timestore)*((endtime+(interval/2))-time)/60,'%.2f') ' minutes)']);
        end
        tocstore=toc;
        timestore=time;
    end
end

delete(f)

% Outputs
soil_depth=(z-zrock);
soil_depth=double(soil_depth);

% Plots
if plots_on==1
    histogram(soil_depth)
    figure;
    surf(x,y,z,soil_depth);
    colorbar;
end
end

