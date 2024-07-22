%% define layer thicknesses
ideltas = [100 50 50 100 50 50 100]'; % give your layer thicknesses
z = sum(ideltas(:)/1000);

%% for radial parabolic GRIN
n1 = 1.565; % maximum index
n2 = 1.555; % minimum index
grinprofile = 'spherical'; % design of the grin profile
discretization = 6; % the discrete steps in change of index
shapeparameters = [2 3 z]; % [aperture, outer diameter, thickness]

%% print configuration
qty = 9;
spacings = [12 15];
offcenter = [0 0];

liftdistance = 5; % mm
liftspeed = 2;    % mm/s
retractspeed = 3; % mm/s
holdtime = 10;    % s

%% other critical settings
numbersoflayerexposures = 5;
refresh = 1;
zinterval = 20;

%% calculation
GRINcalculation = calculateGRIN(ideltas,z,n1,n2,...
    grinprofile,discretization,shapeparameters,...
    qty,spacings,offcenter,liftdistance,...
    liftspeed,retractspeed,holdtime,...
    numbersoflayerexposures,refresh,zinterval);