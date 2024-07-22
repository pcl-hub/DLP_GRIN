function conversionprofileprediction = zconversionprofile_matrixcalculation...
                        (formulation,printer,numberofLayers,...
                         layerThicknesses,exposureTimes,grayscaleValues,...
                         zinterval,refresh,layersofinterest,stage,output)
%%
% this function calculates z conversion profile for one pixel
% Instead of a layer-by-layer loop, matrix calculation is used

%% initialize

f = formulation;

ac = f.resinProperties.criticalConversion;
criteffdose = f.resinProperties.criticalEffectiveDose;

Dp = f.workingCurve.depthofPenetration; % resin depth of penetration
Ec = f.workingCurve.criticalEnergy; % resin critical energy

Eoi = f.cureKinetics.oxygenInhibition.oxygenInhibitionEnergy; % oxygen inhibition energy
ck = f.cureKinetics.cureKineticsModel.fitting; % cure kinetics cfit
w = f.cureKinetics.cureKineticsModel.doseRateDependency.average; % dose rate factor of cure kinetics;

gsf = printer.grayscalefunction; 

% Input printing parameters
numlayers = numberofLayers; % number of layers
deltas = layerThicknesses; % layer thicknesses
taus = exposureTimes; % exposure times
gsvs = grayscaleValues; % grayscale values;

% deltas, taus, and gsvs be n-by-1 vectors
if ~(length(deltas)==numlayers); deltas = deltas(1)*ones(numlayers,1); end
if ~(length(taus)==numlayers); taus = taus(1)*ones(numlayers,1); end
if ~(length(gsvs)==numlayers); gsvs = gsvs(1)*ones(numlayers,1); end

% if specific layers are of interest
lois = sort(layersofinterest);

%% calculate conversion profile

dz = zinterval; % sampling points interval along z
z = flip(dz:dz:sum(deltas(:)))'; % sampling z positions for conversion prediction

zsminus = cumsum(deltas,1,'reverse') - deltas;
zmat = z*ones(1,numlayers)-ones(size(z))*zsminus'; % sampling z positions of all stages

p = zmat>0; % real points

I0 = feval(gsf,gsvs); % incident light intensities of all printing stages
I0(I0<0) = 0; % make sure I is not negative

I0mat = ones(size(z))*I0'; % light intensity matrice for all sampled z positions and all stages
taumat = ones(size(z))*taus'; % exposure time matrice for all sampled z positions and all stages

I1mat = zeros(size(zmat));
I1mat(p) = I0mat(p).*exp(-zmat(p)/Dp);

dosemat = zeros(size(zmat));
dosemat(p) = I1mat(p).*taumat(p); % dose matrice for all sampled points of all stages

I2mat = zeros(size(zmat));
I2mat(I1mat>0) = I1mat(I1mat>0).^(w-1);

if logical(refresh)     % under this setting, the exposure history of the points
    p = dosemat>Eoi;    % below gelation is lost when resin refreshes
    effdosemat = zeros(size(zmat));    
    effdosemat(p) = (dosemat(p)-Eoi).*I2mat(p); % when effective dose is larger than
    p = effdosemat>criteffdose;                 % the critical effective dose, it prints
    p = logical(cumsum(cumsum(p,1,'reverse'),2)); % status of print
    dosemat(~p) = 0; % set the dose of all non-print points to 0 due to refreshing setting
end

if ~isempty(lois) % if only the conversions of specific layers are of interest
    zmat0 = [zeros(size(z)) zmat(:,1:end-1)];
    pois = any((zmat(:,lois)>0 & zmat0(:,lois)<=0),2);

    zmat = zmat(pois,:);
    I2mat = I2mat(pois,:);
    dosemat = dosemat(pois,:);
end

cumdosemat = cumsum(dosemat,2);

oioc = cumdosemat>Eoi; % oxygen inhibition overcame
oioc0 = [zeros(size(oioc(:,1))) double(oioc)];
oioc1 = logical(diff(oioc0,1,2));

dosemat(oioc1) = cumdosemat(oioc1)-Eoi;

effdosemat = zeros(size(dosemat)); % effective doses
effdosemat(oioc) = dosemat(oioc).*I2mat(oioc);

cumeffdosemat = cumsum(effdosemat,2); % accumulated effective doses

if isempty(stage); stage = numlayers; end

% the input stage needs to be empty, 'all', 'final',
% or any number or combination of numbers between 1 and numlayers
switch stage 
    case 'all'
        stage = 1:numlayers;
    case 'final'
        stage = numlayers;
end

z = zmat(:,stage);
cumeffdoses = cumeffdosemat(:,stage);

p = cumeffdoses>0;

% the lines below converts accumulated effective
% doses to conversion values using ck (cure kinetics)
conversions = zeros(size(cumeffdoses));
conversions(p) = feval(ck,cumeffdoses(p)); 
conversions(conversions<0) = 0;

%%
z = z/1000; 

if isempty(output); output = 'full'; end

switch output
    case 'simple'
        conversionprofileprediction = struct('z',z, 'conversionProfile',conversions);

    otherwise
        conversionprofileprediction = struct('formulation',formulation,...
            'printer',printer,'printParameters', table(deltas, taus, gsvs,...
            'VariableNames',["layerThickness","exposureTime","grayscaleValue"]),...
            'refresh',refresh,...
            'prediction', struct('layersofInterest',lois,'stage',stage,...
            'z',z,'effectiveDoses',cumeffdoses,'conversionProfile',conversions),...        
            'units',struct('z','mm','layerThickness','\{mu}m','exposureTime',...
            's','lighIntensity','mW/cm^{2}','grayscaleValue', [],...
            'effectiveDose',['(mW/cm^{2})^{' num2str(w) '} s'],'conversion',[]));
end