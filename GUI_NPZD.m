function [] = GUI_NPZD()
% Graphical User Interface (GUI) for a simple 1D NPZD model based on Kuhn 
% et al. (2015, doi: 10.1016/j.pocean.2015.07.004).
%
% The model is representative of a location at 50 degree N in the North
% Atlantic Ocean and the mixed layer evolution from this location is im-
% posed. The model is run for 2 years, but only the 2nd year is shown in 
% the auto-generated plots. The satellite-observed surface phytoplankton
% evolution is shown for comparison in the surface property plot.
% 
% The GUI allows the user to modify 4 paramters, the latitude of solar 
% forcing, the initial nutrient concentration, the maximum phytoplankton
% growth rate, and the maximum zooplankton grazing rate. These parameters
% can be adjusted by either using the slider to adjust values or by typing 
% them directly into the corresponding editbox. If the number entered is 
% outside the range of the slider, the number will be reset.
%
% The buttons below the four sliders allow the user to run the model and
% quit the GUI. The evolution of state variables at the surface layer will
% be written into a Matlab file along with essential meta-information.
%

% GUI definition
S.fh = figure('units','pixels',...
              'position',[300 300 300 500],...
              'menubar','none',...
              'name','NPZD inputs',...
              'numbertitle','off',...
              'resize','off');
S.tx(1) = uicontrol('style','text',...
                 'unit','pix',...
                 'position',[20 450 260 30],...
                 'string','Latitude (deg N) of solar forcing',...
                 'fontsize',16);
S.sl(1) = uicontrol('style','slide',...
                 'unit','pix',...
                 'position',[20 400 260 30],...
                 'min',20,'max',70,'val',50);             
S.ed(1) = uicontrol('style','edit',...
                 'unit','pix',...
                 'position',[20 430 260 30],...
                 'fontsize',16,...
                 'string','50');
S.tx(2) = uicontrol('style','text',...
                 'unit','pix',...
                 'position',[20 370 260 30],...
                 'string','Initial N (mmol N /m3)',...
                 'fontsize',16);
S.sl(2) = uicontrol('style','slide',...
                 'unit','pix',...
                 'position',[20 320 260 30],...
                 'min',1,'max',15,'val',5);             
S.ed(2) = uicontrol('style','edit',...
                 'unit','pix',...
                 'position',[20 350 260 30],...
                 'fontsize',16,...
                 'string','5'); 
S.tx(3) = uicontrol('style','text',...
                 'unit','pix',...
                 'position',[20 290 260 30],...
                 'string','Maximum P growth rate (1/d)',...
                 'fontsize',16);
S.sl(3) = uicontrol('style','slide',...
                 'unit','pix',...
                 'position',[20 240 260 30],...
                 'min',0.01,'max',2,'val',0.15);             
S.ed(3) = uicontrol('style','edit',...
                 'unit','pix',...
                 'position',[20 270 260 30],...
                 'fontsize',16,...
                 'string','0.15'); 
S.tx(4) = uicontrol('style','text',...
                 'unit','pix',...
                 'position',[20 210 260 30],...
                 'string','Maximum Z grazing rate (1/d)',...
                 'fontsize',16);
S.sl(4) = uicontrol('style','slide',...
                 'unit','pix',...
                 'position',[20 160 260 30],...
                 'min',0.01,'max',2,'val',1.0);             
S.ed(4) = uicontrol('style','edit',...
                 'unit','pix',...
                 'position',[20 190 260 30],...
                 'fontsize',16,...
                 'string','1.0');  
S.pbrun = uicontrol('style','push',...
                 'unit','pix',...
                 'position',[20 90 260 60],...
                 'string','Run the model and save its output',...
                 'fontsize',16,...
                 'callback',{@pbrun_call,S},...
                 'backgroundc',[0 0.6 0.5],...
                 'busyaction','cancel',...% So multiple pushes don't stack.
                 'interrupt','off'); 
S.pbquit = uicontrol('style','push',...
                 'unit','pix',...
                 'position',[20 20 260 60],...
                 'string','Quit',...
                 'fontsize',16,...
                 'callback',{@pbquit_call},...
                 'backgroundc',[0.8 0.09 0.11],...
                 'busyaction','cancel',...% So multiple pushes don't stack.
                 'interrupt','off');             
set([S.ed,S.sl],'call',{@ed_call,S});  % Shared Callback.
end

function [] = pbrun_call(varargin)
% Callback for pushbutton.
[h,S] = varargin{[1,3]};  % Get calling handle and structure.
%h = varargin{1}; % Get the caller's handle.
col = get(h,'backg');  % Get the background color of the figure.
set(h,'str','RUNNING AND SAVING...','fontsize',16,'backg',[1 .6 .6]) % Change color of button. 
% The pause (or drawnow) is necessary to make button changes appear.
% To see what this means, try doing this with the pause commented out.
drawnow;  % FLUSH the event queue.
% Here is where you put whatever function calls or processes that the
% pushbutton is supposed to activate. 
% This is where the NPZD model is called. First, get any values that were 
% reset by the sliders or directly, and display them on the command line.
p=nan*ones(4,1);
for i=1:length(p)
  p(i) = str2double(get(S.ed(i),'string'));
end
disp('You have chosen the following values:')
for i=1:4
  disp([' ' get(S.tx(i),'string') ':  ' get(S.ed(i),'string')])
end
% Next, call the model.
NPZD_main(p)
%
set(h,'str','Run the model','fontsize',16,'backg',col)  % Now reset the button features.
end

function [] = pbquit_call(varargin)
% Callback for pushbutton.
h = varargin{1}; % Get the caller's handle.
col = get(h,'backg');  % Get the background color of the figure.
set(h,'str','QUITTING...','fontsize',16,'backg',[1 .6 .6]) % Change color of button. 
% The pause (or drawnow) is necessary to make button changes appear.
% To see what this means, try doing this with the pause commented out.
drawnow;  % FLUSH the event queue.
% Here is where you put whatever function calls or processes that the
% pushbutton is supposed to activate. 
% In this case we want to quit the model:
close all % closing all figure windows
end



function [] = ed_call(varargin)
% Callback for the edit box and slider.
[h,S] = varargin{[1,3]};  % Get calling handle and structure.

switch h  % Who called?
    case S.ed(1)
        L = get(S.sl(1),{'min','max','value'});  % Get the slider's info.
        E = str2double(get(h,'string'));  % Numerical edit string.
        if E >= L{1} && E <= L{2}
            set(S.sl(1),'value',E)  % E falls within range of slider.
        else
            set(h,'string',L{3}) % User tried to set slider out of range. 
        end
    case S.sl(1)
        set(S.ed(1),'string',get(h,'value')) % Set edit to current slider.
    
    case S.ed(2)
        L = get(S.sl(2),{'min','max','value'});  % Get the slider's info.
        E = str2double(get(h,'string'));  % Numerical edit string.
        if E >= L{1} && E <= L{2}
            set(S.sl(2),'value',E)  % E falls within range of slider.
        else
            set(h,'string',L{3}) % User tried to set slider out of range. 
        end
    case S.sl(2)
        set(S.ed(2),'string',get(h,'value')) % Set edit to current slider.
        
    case S.ed(3)
        L = get(S.sl(3),{'min','max','value'});  % Get the slider's info.
        E = str2double(get(h,'string'));  % Numerical edit string.
        if E >= L{1} && E <= L{2}
            set(S.sl(3),'value',E)  % E falls within range of slider.
        else
            set(h,'string',L{3}) % User tried to set slider out of range. 
        end
    case S.sl(3)
        set(S.ed(3),'string',get(h,'value')) % Set edit to current slider.
        
    case S.ed(4)
        L = get(S.sl(4),{'min','max','value'});  % Get the slider's info.
        E = str2double(get(h,'string'));  % Numerical edit string.
        if E >= L{1} && E <= L{2}
            set(S.sl(4),'value',E)  % E falls within range of slider.
        else
            set(h,'string',L{4}) % User tried to set slider out of range. 
        end
    case S.sl(4)
        set(S.ed(4),'string',get(h,'value')) % Set edit to current slider.
        
    otherwise
        % Do nothing, or whatever.
end
end


function NPZD_main(p)

% This script runs the NPZD model for a single location.  Chlorophyll and
% temperature data are loaded.
%
% Biological Parameters
param_str={'att','alpha','mumax','kN', 'gmax','kP','beta','lPD',...
    'lZN','lZD','lDN','wD','wP'};
% parameter list for bin NA5 form Angela
param = [0.1 0.03 0.15 0.815 1.026 0.5 0.575 0.007 0.02 0.332 0.019 2.438 0.591];
param(3) = p(3); % maximum phytoplankton growth rate 
param(5) = p(4); % maximum zooplankton grazing rate 

%-------------------------------------------------------------------------

load input_data % loads mld, phy, mld_time, temp, temp_time
latitude = p(1);

%-------------------------------------------------------------------------
% Settings related to the temporal resolution
time_res = [2 0.05]; % Time resolution: [nyears dt]

%-------------------------------------------------------------------------
% Settings related to the vertical resolution 
% Note that these correspond to the temperature climatology. Don't change.
Z_res=[5 300]; % depth of each vertical layer and total depth
n = Z_res(2)/Z_res(1); % number of vertical layers 

%-------------------------------------------------------------------------
% Initial conditions
init=[p(2) 0.1 0.01 0.1]; % initial values for N-P-Z-D state variables
T0=ones(n,1)*init;  % initial conditions apply equally for all depths.

%-------------------------------------------------------------------------
% Model integration in time
[State]=npzd_driver(T0,mld,mld_time,temp,temp_time,time_res,Z_res,param,latitude);

%-------------------------------------------------------------------------
% Plotting

nyears = time_res(1); 
dt = time_res(2); 
nsteps=round(365*nyears/dt);
time = (1:nsteps)*dt;
dz = Z_res(1);  
Nz = Z_res(2)/dz;
z = ((0:Nz-1)+0.5)*dz; % mean depth of each layer, in meters

% For surface plot, plot last year in comparison with obs.
ind_sly = 365*(nyears-1)/dt; % index for start of last year

% Calculate daily averages for the last model year (surface only)
var_name = {'Nut' 'Phy' 'Zoo' 'Det'}; 
for i =1:4
  x = squeeze(State(1,i,ind_sly:end));
  x_int = nan(1,365);
  x_int(1) = mean(x(1:1/dt));
  for j=2:365
    x_int(j) = mean(x((j-1)/dt:j/dt));
  end
  eval([var_name{i} '= x_int;'])
end

phy_obs = interp1(mld_time,phy,1:365); % interpolate satellite-based phytoplankton

% plots
xticks=1:30.5:365;

figure 
set (gcf,'position',[902   591   735   359]);
axes('position',[0.1 0.1 0.4 0.85])
hold on
plot(phy_obs,'k','linewidth',3);
plot(Phy(1,:)','color',[0 0.6 0.5],'linewidth',3);
plot(Zoo(1,:)','color',[0.13 0.37 0.66],'linewidth',3);
legend('P observations','Model P','Model Z','location','northwest')
legend('boxoff')
set(gca,'fontsize',12,'ylim',[0 0.6], 'xlim',[1 365],'box','on')
datetick('x','mmm','keeplimits','keepticks')
ylabel('Surface concentrations (mmol N m^{-3})')
%
axes('position',[0.55 0.1 0.4 0.85])
hold on
plot(Nut(1,:)','color',[0.8 0.09 0.11],'linewidth',3);%
plot(Det(1,:)','color',[0.8 0.11 0.55],'linewidth',3);%
legend('Model N','Model D')
legend('boxoff')
set(gca,'fontsize',12, 'xlim',[1 365],'box','on')
datetick('x','mmm','keeplimits','keepticks')

figure
set (gcf,'position',[1000   564   568  774]);
set (gcf,'position',[1000   564   350  774]);
var_str = {'N' 'P' 'Z' 'D'};
colormap(roma_inv13)
for var=1:4
    St=squeeze(State(:,var,:)); 
    subplot(4,1,var); 
    pcolor(time,-z,St); shading flat; hcb = colorbar; hold on
    ylabel('Depth (m)','fontsize',14); 
    text(130,-180,[var_str{var} ' (mmol N m^{-3})'],'fontsize',14,'fontweight','bold','color','w');
    set(gca,'fontsize',14,'layer','top','tickdir','out','xlim',[0 365]);
    if var<4
        set(gca,'ylabel',[],'xlabel',[],'yticklabel',[],'xticklabel',[])
    end
    if var==4
        datetick('x','mmm','keeplimits','keepticks')
    end
    set(hcb,'location','south','color',[1 1 1])
    hold on; plot([mld_time(2:end-1); mld_time(2:end-1)+365],repmat(-mld(2:end-1),2,1),'w-','linewidth',2); hold off
end
set(gcf,'PaperPosition',[0.25 2.5 5 9])

fn = ['NPZD_run_' datestr(now,30)];
save(fn,'param','param_str','latitude','Nut','Phy','Zoo','Det','phy_obs')
disp(['Saved output to ' fn '.mat']) 
disp('  ')

end


function [State]=npzd_driver(T0,mld,mld_time,temp,temp_time,time_res,Z_res,bio_param,latitude)

% This function is the main driver to run the model
%
% This function is called by:
%                         -- GUI_NPZD
%
% This function calls:
%                         -- physics (physical model)
%                         -- surface_light (determination of incoming light)
%                         -- biology (biological model)
% Inputs:
%       T0 (initial conditions)
%       mld (mixed layer depth vector)
%       temp (temperature array)
%       time_res (1 by 2 vector with temporal resolution: simulation length
%                 in years and time step)
%       Z_res(1 by 2 vector with vertical resolution: total depth and bin depth)
%       bio_param (biological parameters)
%       latitude
%
% Outputs:
%       State (model results, 4 state variables)
%
%-------------------------------------------------------------------------
%------------------------------------------------------------------------
%

tic

% Define variables
nvar=4;
% variable 1 ='nutrients'; 
% variable 2 ='phytoplankton'; 
% variable 3 ='zooplankton'; 
% variable 4 ='detritus'; 

% Define  resolution
nyears = time_res(1);               % number of years of run
dt = time_res(2);                   % time step (days)
nsteps = round(365*nyears/dt);      % number of time steps

dz = Z_res(1);                      % number of levels
n = Z_res(2)/dz;                    % vertical grid spacing (in m)

% Model run
State = nan*ones(n,nvar,nsteps); T=T0;
time = (1:nsteps)*dt;

for k=1:nsteps
    
    MLD  = round(interp1(mld_time,mld,mod(time(k),365))); % Mixed layer depth in meters
    temp_prof = interp1(temp_time,temp',mod(time(k),365))'; % interpolate temp profile in time
    
    T = physics(T,dt,dz,MLD); % Vertically diffuse
    Io= surface_light(latitude,k*dt); % Incoming light
    T = biology(T,dz,dt,bio_param,temp_prof,Io); % Runs biological model 
                                
  if isnan(T)
       display(MLD);
       display(TEMP);
       display(k);
       error('NaN results')
  end   
    
% Store state all variables
    State(:,:,k)=T;          % NPZD                                   
end
 
%--------------------------------------------------------------------------

elapsedtime=toc;
fprintf('Elapsed time is: %3.2f seconds\n', elapsedtime)
end

function T=physics(T,dt,dz,MLD)

n=size(T,1);                                               % number of levels

% Values of diffusivities above and below mld
nu1=0.0088*60*60*24;                                       % upper layer diffusivities (Converted from m2/s to m2/d)
nu2= 0.000014*60*60*24;  

% Define vector of diffusivities for all (n+1) cell boundaries
Zcell=dz*[0:n]';                                           % Define cell boundaries
Nu=nu1*ones(n+1,1); Nu(Zcell>MLD)=nu2;

% Define Crank-Nicolson matrices 
M=diag(-[Nu(1:n)+Nu(2:end)])+diag(Nu(2:end-1),-1)+diag(Nu(2:end-1),1);
M(1,1)=-Nu(2);                                             % no flux top BC
M(end,end)=-Nu(end-1);                                     % no flux bottom BC
M1= eye(n)-dt/2/dz^2*M;
M2= eye(n)+dt/2/dz^2*M;

% Update the state variables
T=M1\(M2*T);
end

function Io = surface_light(latitude,time)
%
% This function calculates the hourly radiation from the formula 
% Io = SC * ECC .* C_theta * PARfrac * cloud 
% Where SC = solar constant
%       ECC = eccentricity correction (between 0.9 and 1.1 approx.)
%       C_theta = cosine of the solar zenith angle, as function of
%                 latitude, time of the year and hour of the day
%       PARfrac = photosynthetically active radiation fraction 0.43 of total
%                 incoming light
%       cloud = fraction of light available after cloud attenuation
%

% Constants
SC = 1366.1;  % solar constant W m^{-2}
PHI = latitude*2*pi/360; % latitude angle 
PARfrac = 0.43; % PAR fraction
cloud = 0.3; 

% Eccentricity correction
%  calculated using the Spencer's fourier series.
gam = 2*pi*time./365;  % time of the year in radians
ECC = 1.000110 + (0.034221 * cos(gam)) + (0.001280 * sin (gam)) + ...
    (0.000719*cos(2*gam) + (0.000077 * sin(2*gam))); % Eccentricity: correction factor of Earth's orbit (~0.9 - ~1.1)

% Calculate declination angle
delta = 0.006918 - 0.399912 * cos(gam) + 0.070257 * sin (gam) - ...
    0.006758 * cos(2* gam) + 0.000907 * sin (2*gam) - 0.002697 * cos(3*gam) +...
    0.001480*sin(3*gam);

% Calculate hour angle
%  The hour angle is equal to zero at local noon and increases in magnitude
%  by pi/12 or 15 degrees for every hour before and after noon.
hour = 24*mod(time,1);
h = pi/12 * (12 - hour);

% Calculate cosine of the solar zenith angle (C_theta)
C_theta = sin(PHI) * sin(delta)+ cos(PHI) * cos(delta) .* cos(h);

% Calculate hourly solar radiation
Io = SC * ECC.* C_theta * PARfrac * (1-cloud); % PAR after attenuation by clouds
night = Io <0;% Set night hours to zero
Io(night)=0;
end

function [T]=biology(T,dz,dt,bio_param,temp_prof,Io) 
%
% This function contains the biological model.
%

%------------------------------------------------------------------------- 
% Biological state variables
inut=1; iphy=2; izoo=3; idet=4;

% Biological parameters
att   = bio_param(1); % light attenuation (att), units are per meter
alpha = bio_param(2); % initial slope of P-I curve (alpha) mmol N (W m-2)-1 d-1
mu0   = bio_param(3); % maximum growth rate (mumax) at T =0 C , units are per day
kN    = bio_param(4); % half-saturation of nutrients uptake(kN), units are mmol N m-3
gmax  = bio_param(5); % maximum grazing rate (gmax), units are per day
kp    = bio_param(6); % half-saturation of grazing(kp), units are mmol N m-3
beta  = bio_param(7); % assimilation efficiency (beta), dimensionless
lPD   = bio_param(8); % mortality of phytoplantkon (lPD), per day
lZN   = bio_param(9); % excretion of zooplankton (lZN), per day
lZD   = bio_param(10);% mortality of zooplankton (lZD), per day
lDN   = bio_param(11);% remineralization rate(lDN), per day
wD    = bio_param(12);% sinking speed of detritus (wD), meters per day
wP    = bio_param(13);% sinking speed of phytoplantkon (wP), meters per day
%-------------------------------------------------------------------------

Nz=size(T,1);
dz2 = 0.5*dz;      % half a grid box
kz = 1:Nz;
z = (kz-1)*dz+dz2; % mean depth of each layer, in meters

%-------------------------------------------------------------------------
% Temperature dependency
%-------------------------------------------------------------------------
mumax = (1.066.^temp_prof).*mu0; % temperature-dependend maximum growth rate (Eppley,1972)

%-------------------------------------------------------------------------
% Light
%-------------------------------------------------------------------------
I = Io * exp(-att*z'); % light attenuation with depth according to Beer's Law

%-------------------------------------------------------------------------
% Phytoplankton Growth
%-------------------------------------------------------------------------
PHI = alpha.*I./sqrt(mumax.*mumax+alpha.*alpha.*I.*I); % Light limitation
cff = dt*mumax.*PHI.*T(:,iphy)./(kN+T(:,inut)); % Growth rate * dt
T(:,inut) = T(:,inut)./(1+cff); % Updated nutrients
T(:,iphy) = T(:,iphy)+cff.*T(:,inut); % Updated phytoplankton

%--------------------------------------------------------------------------
% Zooplankton Grazing
%--------------------------------------------------------------------------
cff = dt*gmax*T(:,izoo).*T(:,iphy)./(kp+T(:,iphy).*T(:,iphy));
T(:,iphy) = T(:,iphy)./(1+cff); % Phytoplankton - grazing
T(:,izoo) = T(:,izoo)+cff.*T(:,iphy).*beta; % Zooplankton + Assimilated grazing fraction
non_assim = cff.*T(:,iphy).*(1-beta); % Unassimilated grazing fraction

%--------------------------------------------------------------------------
% Zooplankton Loss Terms
%--------------------------------------------------------------------------
% Mortality and (1-beta) go to Detritus
cff = dt*lZD*max(T(:,izoo)-.001,0.0).^2; % Mortality of zooplankton 
T(:,izoo) = T(:,izoo) - cff; % Zooplankton - Mortality
T(:,idet) = T(:,idet) + cff + non_assim; % Detritus + Zoo Mortality + Unassimilated Grazing

% Excretion goes to Nutrients
cff = dt*lZN*max(T(:,izoo)-.001,0.0); % Metabolic losses
T(:,izoo) = T(:,izoo) - cff;
T(:,inut) = T(:,inut) + cff;

%--------------------------------------------------------------------------
% Phytoplankton Loss Terms
%--------------------------------------------------------------------------
% Mortality of phytoplantkton goes to Detritus
cff = dt*lPD*max(T(:,iphy)-.001,0.0); % same units of phyto
T(:,iphy) = T(:,iphy) - cff; % Updated phytoplankton
T(:,idet) = T(:,idet) + cff; % Updated detritus

%--------------------------------------------------------------------------
% Remineralization
%--------------------------------------------------------------------------
cff = dt*lDN*max(T(:,idet)-.001,0.0);
T(:,idet) = T(:,idet) - cff;
T(:,inut) = T(:,inut) + cff;

%--------------------------------------------------------------------------
% Sinking of phytoplankton and detritus.
%--------------------------------------------------------------------------
phyMin = 0.001;
cff = [phyMin; T(1:end-1,iphy)];
T(end,inut) = T(end,inut) + dt*wP*(T(end,iphy)-phyMin)/dz; 
T(:,iphy)   = T(:,iphy) - dt*wP*(T(:,iphy)-cff)/dz; 

cff = [phyMin; T(1:end-1,idet)];
T(end,inut) = T(end,inut) + dt*wD*(T(end,idet)-phyMin)/dz;
T(:,idet)   = T(:,idet) - dt*wD*(T(:,idet)-cff)/dz;
end


