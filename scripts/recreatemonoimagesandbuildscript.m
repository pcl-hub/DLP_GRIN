function recreatemonoimagesandbuildscript(GRINcalculation, savetopath)

images = GRINcalculation.GRIN.baseImages;
imsize = size(images); if length(imsize) == 2; imsize(3) = 1;end
images = reshape(images,[imsize(1)*imsize(2), imsize(3)]);

[fullsize, pixelsize, qty, spacings, offcenter, liftdistance,liftspeed,...
    retractspeed, holdtime] = extractconfiguration(GRINcalculation);

maxtau = 1000; % assumes maximum single exposure time possible is 1000s

%%
for p = 1:length(GRINcalculation.calculations.predictions)
    
    prediction = GRINcalculation.calculations.predictions(p);

    numexposure = prediction.maxNumberofExposuresPerLayer;

    calculatedParameters = prediction.calculatedParameters;

    % numexposures = calculatedParameters.numExposures;
    
    deltas= calculatedParameters.layerThickness;
    taus = calculatedParameters.exposureTime;
    gsvs = calculatedParameters.grayscaleValue;
    
    bppositions = (cumsum(deltas,1)- deltas(1,1))/1000; 
    bppositions(end,:) = bppositions(end,:) + 0.05;
    % 0.05 is the thickness difference 
    % btw. a glass slide and a FEP film
    
    [sortedtaus,inds] = sort(taus,2,'descend');
    stagetaus = sortedtaus-[sortedtaus(:,2:end) zeros(size(taus(:,1)))];
    
    layercounts = sum(ceil(stagetaus/maxtau),2);
    totallayercount = sum(layercounts);

    liftdistances = zeros(totallayercount,1);
    liftdistances(cumsum(layercounts,1),1) = liftdistance;
    
    %%
    fpath = [savetopath 'prediction' num2str(p) '\files\'];
    mkdir(fpath); monobuildscriptfileheader(fpath,totallayercount); 

    fileID = fopen([fpath 'buildscript.ini'],'a');

    layercount = 0;

    for i = 1:size(taus,1)
    
        gsvs1 = ones(size(images(:,1)))*gsvs(i,inds(i,:));
        
        images1 = cumsum(images(:,inds(i,:)).*gsvs1,2);

        images1 = reshape(images1,imsize);

        m = ceil(stagetaus(i,:)/maxtau);
        
        ni = find(m);

        for j = ni % or flip(ni)

            m1 = m(j);
            bpposition1 = bppositions(i,j);
            exposuretime1 = ones(m1,1)*maxtau;
            exposuretime1(m1,1) = stagetaus(i,j)-(m1-1)*maxtau;

            image1 = images1(:,:,j);
            fullimage1 = placetheimages(image1,...
                qty,spacings,offcenter,fullsize,pixelsize);
            
            for j1 = 1:m1
                layercount = layercount+1;
                
                filename1 = monoimagename(layercount);
                fullfilename1 = [fpath filename1 '.png'];
                imwrite(uint8(fullimage1),fullfilename1,'png');
            
                fprintf(fileID,'%1$.4f,%2$ 10s,%3$ .1f,%4$ .1f,%5$ .1g,%6$ .1f\n',...
                 bpposition1, filename1,exposuretime1(j1,1),liftdistances(layercount,1),0,0);
            end

        end
    end
    
    fclose(fileID);

    %%
    interlayertime = (liftdistance/liftspeed+...
        liftdistance/retractspeed+holdtime+1)*(size(taus,1)-1);
    totalprinttime = (sum(stagetaus(:))+interlayertime)/3600;
    lastlayerprinttime = sum(stagetaus(end,:))/60;
    
    disp(['Prediction #' num2str(p) ':']);
    disp(['max number of exposures per layer: ' num2str(numexposure)]);
    disp(['total number of exposures is ' num2str(layercount) ';']);
    disp(['total print time will be ' num2str(totalprinttime) ' hours; ']);
    disp(['last layer will take ' num2str(lastlayerprinttime) ' minutes.']);

end



function [fullsize, pixelsize, qty, spacings, offcenter, liftdistance,liftspeed,...
    retractspeed, holdtime] = extractconfiguration(GRINcalculation)

fullsize = GRINcalculation.printer.resolution;
pixelsize = GRINcalculation.printer.pixelSize;

qty = GRINcalculation.printConfiguration.locationInfo.quantity;
spacings = GRINcalculation.printConfiguration.locationInfo.spacings;
offcenter = GRINcalculation.printConfiguration.locationInfo.offcenter;

liftdistance = GRINcalculation.printConfiguration.interlayerParameters.liftDistance;
liftspeed = GRINcalculation.printConfiguration.interlayerParameters.liftSpeed;
retractspeed = GRINcalculation.printConfiguration.interlayerParameters.retractSpeed;
holdtime = GRINcalculation.printConfiguration.interlayerParameters.holdTime;


function monobuildscriptfileheader(fpath,totallayercount)

fileID = fopen([fpath 'buildscript.ini'],'w');
fprintf(fileID, '%1$20s\n','Machine = MONO3_64UM');
fclose(fileID);

strnumber = floor(log10(totallayercount)+20);
strformat = ['%1$' num2str(strnumber) 's\n'];
strlc = num2str(totallayercount);

fileID = fopen([fpath 'buildscript.ini'],'a');
fprintf(fileID, '%1$21s\n','Slice thickness = 100');
fprintf(fileID, strformat,['number of slices = ' strlc]);
fprintf(fileID, '%1$22s\n','illumination time = 20');
fprintf(fileID, '%1$29s\n','number of override slices = 0');
fprintf(fileID, '%1$31s\n','override illumination time = 20');
fprintf(fileID, '%1$25s\n','support burn in time = 20');
fclose(fileID);


function filename = monoimagename(layernumber)
numdigits = floor(log10(layernumber)+1);
monoprefix = {'S00000' 'S0000' 'S000' 'S00' 'S0'};
filename = [monoprefix{numdigits} num2str(layernumber) '_P1'];