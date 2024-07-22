function calculation = calculateprintparameters(formulation, printer, ...
                        conversionmatrix, ideltas, itaus, igsvs, ...
                        numbersoflayerexposures, refresh, zinterval)
%%
% conversionmatrix is a 2D numerical matrix.
% each column of the matrix represents a conversion profile along z.
% the number of columns (Dim. 2) corresponds to the number of
% different z conversion profiles.

% prefix i is for input
% prefix c is for current
% no prefix for predicted parameters

totallayerthickness = sum(ideltas,1);
targetthickness = zinterval*size(conversionmatrix,1);

if totallayerthickness ~= targetthickness
    error(sprintf(['The layer thicknesses do not '...
                    'add up to the target thickness.\n'...
                    'Please enter proper layer thicknesses'...
                    ' and try again.']));
end

%%
DCT = 'Distrib_Computing_Toolbox';
if license('test', DCT); pool = gcp; else; pool = []; end
clear DCT

%%
numlayers = size(ideltas,1);
numprofiles = size(conversionmatrix,2);

s = [length(itaus) length(igsvs) numprofiles];
maxnumexposure = min([s(1) s(3)]);

numexpos = numbersoflayerexposures;
if isempty(numexpos); numexpos = s(1); end    % the number of exposures can
numexpos = numexpos(:)';                      % be 1, 2, or 3, or larger; 
p = numexpos<1 | numexpos>3;                  % if it is empty, < 1, or > 3, 
numexpos(p) = s(1);                           % the number will be determined
numexpos = unique(numexpos);                  % by the algorithm.

maxnumexpos = numexpos;
p = maxnumexpos == s(1);
maxnumexpos(p) = maxnumexposure;

m = length(numexpos);

numexposures = zeros(numlayers,m);

deltas = ideltas;
taus = zeros(numlayers,numprofiles,m);
gsvs = zeros(numlayers,numprofiles,m);

minsses = zeros(numlayers,numprofiles,m);

z = flip(zinterval:zinterval:sum(deltas(:)))';
zsminus = cumsum(deltas,1,'reverse') - deltas;
zmat = z*ones(1,numlayers)-ones(size(z))*zsminus';

p = zmat>0;
p0 = [false(size(p(:,1))) p(:,1:end-1)];
p = p & ~p0;

for j = 1:m

    for layernum = numlayers:-1:1
    
        sses = zeros(s);
        p1 = p(:,layernum);

        cdeltas = deltas(layernum:end,1);
        ctaus = taus(layernum:end,:,j);
        cgsvs = gsvs(layernum:end,:,j);        
    
        targetconversions = conversionmatrix(p1,:);    

        if isempty(pool)
            for i = 1:numprofiles
                currentParams = [cdeltas ctaus(:,i) cgsvs(:,i)];
                sses(:,:,i) = calculateSSE(formulation, printer,...
                    targetconversions(:,i), currentParams,...
                    itaus, igsvs, refresh, zinterval);
            end
        else
            parfor i = 1:numprofiles
                currentParams = [cdeltas ctaus(:,i) cgsvs(:,i)];
                sses(:,:,i) = calculateSSE(formulation, printer,...
                    targetconversions(:,i), currentParams,...
                    itaus, igsvs, refresh, zinterval);
            end
        end
        
        ssespermute = permute(sses,[1 3 2]);
        [sses2d, inds] = min(ssespermute,[],3);
    
        [tinds,minsse] = findminimum(sses2d,numexpos(j));
        
        ginds = inds(((1:s(3))-1)*s(1)+tinds);
        
        taus(layernum,:,j) = itaus(tinds);
        gsvs(layernum,:,j) = igsvs(ginds);
        
        minsses(layernum,:,j) = minsse;

        pt = taus(layernum,:,j)>0;
        numexposures(layernum,j) = numel(unique(taus(layernum,pt,j)));
    end
end

predictedconversions = zeros([size(conversionmatrix) m]);
predictions(m) = struct;

for j = 1:m
    if isempty(pool)
        for i = 1:numprofiles
            prediction1 = zconversionprofile_matrixcalculation...
                            (formulation, printer,numlayers,...
                            deltas,taus(:,i,j),gsvs(:,i,j),...
                            zinterval,refresh,[],[],'simple');
            predictedconversions(:,i,j) = prediction1.conversionProfile;
        end
    else
        parfor i = 1:numprofiles
            prediction1 = zconversionprofile_matrixcalculation...
                            (formulation, printer,numlayers,...
                            deltas,taus(:,i,j),gsvs(:,i,j),...
                            zinterval,refresh,[],[],'simple');
            predictedconversions(:,i,j) = prediction1.conversionProfile;
        end
    end

    predictions(j).maxNumberofExposuresPerLayer = maxnumexpos(j);
    predictions(j).calculatedParameters = table(numexposures(:,j),...
        deltas*ones(1,numprofiles),taus(:,:,j),gsvs(:,:,j),...
        'VariableNames',["numExposures","layerThickness",...
        "exposureTime","grayscaleValue"]);
    predictions(j).predictedConversionMatrix = table(z,...
        predictedconversions(:,:,j),'VariableNames',["z","conversion"]);
    predictions(j).sse = minsses(:,:,j);

end

delete(pool);

z = z/1000;

calculation  = struct('formulation',formulation,'printer',printer,...
    'refresh',logical(refresh),'targetConversionMatrix',...
    table(z,conversionmatrix,'VariableNames',["z","conversion"]),...
    'predictions',predictions,'units',struct('z','mm','layerThickness',...
    '\{mu}m','exposureTime','s','grayscaleValue',[],'conversion',[]));


function sse = calculateSSE(formulation, printer, targetconversionprofile,...
                            currentParams, itaus, igsvs, refresh, zinterval)

tcp = targetconversionprofile;
sse = zeros(length(itaus), length(igsvs));

for i = 1:length(itaus)
    for j = 1:length(igsvs)
        currentParams(1,[2 3]) = [itaus(i) igsvs(j)];
        numlayers = size(currentParams,1);
        cdeltas = currentParams(:,1);
        ctaus = currentParams(:,2);
        cgsvs = currentParams(:,3);

        calculation = zconversionprofile_matrixcalculation...
            (formulation, printer,numlayers,cdeltas,ctaus,cgsvs,...
             zinterval,refresh,1,'final','simple');

        cp = calculation.conversionProfile;
        sse(i,j) = (cp-tcp)'*(cp-tcp);
    end
end


function [tauinds,minsse] = findminimum(sse,numexpo)

rowcombs = nchoosek(1:size(sse,1),numexpo);
ssse = zeros(size(rowcombs,1),1);

for i = 1:size(rowcombs,1)
    sse1 = sse(rowcombs(i,:),:);
    minsse1 = min(sse1,[],1);
    ssse(i) = sum(minsse1(:));
end

[~,i] = min(ssse,[],1);
sse1 = sse(rowcombs(i,:),:);
[minsse,tauinds1] = min(sse1,[],1);

tauinds = rowcombs(i,tauinds1);