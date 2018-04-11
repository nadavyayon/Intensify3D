function Intensify3D(DirectoryOf16bit,MinIm,MaxIm,STDNumber,MaxTissueIntensity,FilterSize,HasBackground,NormalizeType,Threads,ImageToSelectFrom,handles) 
tic % start Stopwatch
%% initiation
warning('off');
Dialog('Getting Started....',handles);
if Threads
    if ~isempty(gcp('nocreate'));
%         delete(gcp); 
    else
        parpool(Threads);
    end   
end
cd(DirectoryOf16bit);
% delete old files
Folder16bitInfo = dir('*.tif');
if exist('NormalizedBackground') == 0
    mkdir('NormalizedBackground'); 
else
    cd('NormalizedBackground');
    toDelete = dir;
    for k = 3 : length(toDelete)
        FileName = toDelete(k).name;
        delete(FileName);
    end
    if exist('SupportImages')==0 && HasBackground
        mkdir('SupportImages');
    else
        if HasBackground 
             cd('SupportImages');
             toDelete = dir;
             for k = 3 : length(toDelete)
            FileName = toDelete(k).name;
            delete(FileName);
        end
        end
    end    
end
cd(DirectoryOf16bit);
%get image size
FirstImageFileIndex = find([Folder16bitInfo.bytes]>10*1000,1,'first');
Ti = Tiff([Folder16bitInfo(FirstImageFileIndex+MinIm).name],'r');
ImMd = Ti.read();
ImD = double(ImMd);
ImSize = size(ImD);
Ti.close();
FilterXY = [FilterSize FilterSize];

%% init background estimation
if HasBackground
    cd('NormalizedBackground');
    mkdir('SupportImages');
    if HasBackground<3
        [Inn , ~, ~] = neighbourND( 1:ImSize(1)*ImSize(2), ImSize, [1 1] );
        Inn((Inn==0)) = 1; % define edges as 1
        BackroundThreshold = [];
    else 
        Inn = [];
        % measure statistics of 10 readomly selected images 
        SelectedImages = 1:ceil(length(Folder16bitInfo)/20):length(Folder16bitInfo);
        counter = 0;
        for im = SelectedImages
            counter = counter + 1;
            Ti = Tiff([DirectoryOf16bit,'' filesep '',Folder16bitInfo(FirstImageFileIndex+im-1).name],'r');
            ImD = double(Ti.read());
            Ti.close();
            ImD_Downsampled = medfilt2(ImD(1:10:end,1:10:end));
            ThresholdValue(counter) = SimpleThrehold(ImD_Downsampled);
        end
        BackroundThreshold = quantile(ThresholdValue,1-(STDNumber+ 0.001)/10.1);        
    end
end
clc
%% estimate quantile from selcted Image
Ti = Tiff([DirectoryOf16bit,'' filesep '',Folder16bitInfo(FirstImageFileIndex+ImageToSelectFrom-1).name],'r');
ImMd = Ti.read();
% Convert to double 
ImD = double(ImMd);
% find quantile 
NoOfQuantiles = 10000; 
ImageQuantiles = quantile(ImD(:),NoOfQuantiles);
[~,UserDefinedQuantile] = min(abs(ImageQuantiles-MaxTissueIntensity)); % find quantile 
UserDefinedQuantile = UserDefinedQuantile/NoOfQuantiles;
clc
      
%% Normalize along the XY plane
% image statistics
EstimatedMaxTissueIntensitySeriers = zeros(length(MinIm:MaxIm),1); % for upper precentile
SemiQuantiles = zeros(length(MinIm:MaxIm),NoOfQuantiles); % for Semi-Qantile 
LowerUpperQuantiles = zeros(length(MinIm:MaxIm),2);% for Contrast Stretch
if HasBackground
    % Estimate Background
    Dialog('Estimating Background Area....',handles); 
    parfor_progress(MaxIm+1-MinIm);
    if Threads
        parfor i = MinIm:MaxIm         
            warning('off','all');
            % init
            ImSizePar = ImSize;
            Folder16bitInfoPar = Folder16bitInfo;
            Ti = Tiff([DirectoryOf16bit,'' filesep '',Folder16bitInfoPar(FirstImageFileIndex+i-1).name],'r');
            ImD = double(Ti.read());
            Ti.close();
            % find quantile value for each image individually 
            ImageQuantiles = quantile(ImD(:),NoOfQuantiles);
            if UserDefinedQuantile==1  % enable the user to define a value larger than the brightest pixel
                EstimatedMaxTissueIntensity = max(ImD(:)); 
            else
                EstimatedMaxTissueIntensity = ImageQuantiles(UserDefinedQuantile*NoOfQuantiles)
            end
            FindBackground(i,DirectoryOf16bit,HasBackground,ImD,EstimatedMaxTissueIntensity,Inn,BackroundThreshold,ImSizePar,FilterXY,STDNumber)
            parfor_progress;
        end
    else
        for i = MinIm:MaxIm  
            warning('off','all');
            % init
            ImSizePar = ImSize;
            Folder16bitInfoPar = Folder16bitInfo;
            Ti = Tiff([DirectoryOf16bit,'' filesep '',Folder16bitInfoPar(FirstImageFileIndex+i-1).name],'r');
            ImD = double(Ti.read());
            Ti.close();

            % find quantile value for each image individually 
            ImageQuantiles = quantile(ImD(:),NoOfQuantiles);
            if UserDefinedQuantile==1  % enable the user to define a value larger than the brightest pixel
                EstimatedMaxTissueIntensity = max(ImD(:)); 
            else
                EstimatedMaxTissueIntensity = ImageQuantiles(UserDefinedQuantile*NoOfQuantiles)
            end     
            FindBackground(i,DirectoryOf16bit,HasBackground,ImD,EstimatedMaxTissueIntensity,Inn,BackroundThreshold,ImSizePar,FilterXY,STDNumber)
            parfor_progress;
        end
    end
    warning('off','all');
    % Z filter background estimation
    if (MaxIm-MinIm > 7)
    Dialog('Smoothing Background Estimation..... ',handles);
    if Threads
        parfor i = MinIm:MaxIm
            if (i<(MinIm+2)); ii = [0 0 0 1 2]; end
            if (i>=(MinIm+2))&&(i<=(MaxIm-2)); ii = [-2 -1 0 1 2]; end
            if (i>(MaxIm-2)); ii = [-2 -1 0 0 0]; end
            Ts1 = Tiff([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'SupportImages' filesep 'Supp_Image_', sprintf('%3.4d',i+ii(1))],'r');
            ImS1 = Ts1.read();
            Ts2 = Tiff([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'SupportImages' filesep 'Supp_Image_', sprintf('%3.4d',i+ii(2))],'r');
            ImS2 = Ts2.read();
            Ts = Tiff([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'SupportImages' filesep 'Supp_Image_', sprintf('%3.4d',i+ii(3))],'r+');
            ImS = Ts.read();          
            Ts4 = Tiff([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'SupportImages' filesep 'Supp_Image_', sprintf('%3.4d',i+ii(4))],'r');
            ImS4 = Ts4.read();
            Ts5 = Tiff([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'SupportImages' filesep 'Supp_Image_', sprintf('%3.4d',i+ii(5))],'r');
            ImS5 = Ts5.read();
            ImS = (ImS1+ImS2+ImS+ImS4+ImS5)>3;
            Ts1.close(); Ts2.close(); Ts.close();Ts4.close(); Ts5.close();
            SaveTiffInCode([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'SupportImages' filesep 'Tissue_Area_', sprintf('%3.4d',i)],1,'w',ImS)
        end
    else 
        for i = MinIm:MaxIm       
            if  (i<(MinIm+2)); ii = [0 0 0 1 2]; end
            if (i>=(MinIm+2))&&(i<=(MaxIm-2)); ii = [-2 -1 0 1 2]; end
            if (i>(MaxIm-2)); ii = [-2 -1 0 0 0]; end    
            Ts1 = Tiff([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'SupportImages' filesep 'Supp_Image_', sprintf('%3.4d',i+ii(1))],'r');
            ImS1 = Ts1.read();
            Ts2 = Tiff([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'SupportImages' filesep 'Supp_Image_', sprintf('%3.4d',i+ii(2))],'r');
            ImS2 = Ts2.read();
            Ts = Tiff([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'SupportImages' filesep 'Supp_Image_', sprintf('%3.4d',i+ii(3))],'r+');
            ImS = Ts.read();          
            Ts4 = Tiff([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'SupportImages' filesep 'Supp_Image_', sprintf('%3.4d',i+ii(4))],'r');
            ImS4 = Ts4.read();
            Ts5 = Tiff([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'SupportImages' filesep 'Supp_Image_', sprintf('%3.4d',i+ii(5))],'r');
            ImS5 = Ts5.read();
            ImS = (ImS1+ImS2+ImS+ImS4+ImS5)>3;
            Ts1.close(); Ts2.close(); Ts.close();Ts4.close(); Ts5.close();
            SaveTiffInCode([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'SupportImages' filesep 'Tissue_Area_', sprintf('%3.4d',i)],1,'w',ImS)
        end
    end
    for i = MinIm:MaxIm; delete([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'SupportImages' filesep 'Supp_Image_', sprintf('%3.4d',i)]); end
    %close(h)
    end
    
% Normalize background detection
    Dialog('Correcting Background across XY plane....',handles);
    parfor_progress(MaxIm+1-MinIm)
    if Threads 
        parfor i = MinIm:MaxIm
            % Load Images
            FilterXYPar = FilterXY;
            Folder16bitInfoPar = Folder16bitInfo;
            Ti = Tiff([DirectoryOf16bit,'' filesep '',Folder16bitInfoPar(FirstImageFileIndex+i-1).name],'r');
            ImD = double(Ti.read());
            Ti.close();
            % support
            Ts = Tiff([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'SupportImages' filesep 'Tissue_Area_', sprintf('%3.4d',i)],'r');
            ImS = double(Ts.read());
            ImSupp = Ts.read();
            Ts.close();        
            % find quantile value for each image individually 
            ImageQuantiles = quantile(ImD(:),NoOfQuantiles);
            if UserDefinedQuantile==1  % enable the user to define a value larger than the brightest pixel
                 EstimatedMaxTissueIntensity = max(ImD(:)); 
            else
                 EstimatedMaxTissueIntensity = ImageQuantiles(UserDefinedQuantile*NoOfQuantiles)
            end
            % Remove_Peaks and zeros 
            ImD_Smooth = ImD;
            Maxes = (ImD > EstimatedMaxTissueIntensity);
            TissuePixelvalues = ImD(logical(ImS.*(ImD_Smooth<quantile(ImD_Smooth(:),0.9)).*(ImD_Smooth>quantile(ImD_Smooth(:),0.1))));
            if size(TissuePixelvalues(:),1)>100     
               RandValues = emprand(TissuePixelvalues,sum(sum(Maxes)),1); % generate signal-replacing values
               ImD_Smooth(Maxes)= RandValues;   
            end     
            % Smooth along XY
            ImD_Smooth(not(ImS)) = EstimatedMaxTissueIntensity/2;%quantile(ImD(TissuePixelIndices),UserDefinedQuantile);
            ImD_Smooth = Savitzky_Golay(ImD_Smooth,FilterXYPar(1),FilterXYPar(2),1);
            ImD_Smooth = Savitzky_Golay(ImD_Smooth,FilterXYPar(1),FilterXYPar(2),1); % 2nd time
            ImD_Smooth_Norm = ImD_Smooth/max(ImD_Smooth(:));   
            NormXY = (ImD./ImD_Smooth_Norm);
            SaveTiffInCode([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'Norm_',Folder16bitInfo(FirstImageFileIndex+i-1).name],16,'w',NormXY);
            NewBackgroundEstimateImage = NormXY;
            NewMBI = min(NewBackgroundEstimateImage(Maxes));            
            % log Statistics 
            % for upper precentile
            EstimatedMaxTissueIntensitySeriers(i) = EstimatedMaxTissueIntensity;
            % for Semi-Qantile
            ImLin = ImSupp(:)>0;
            if UserDefinedQuantile ~= 1 
               UpperQuantile = NewMBI; 
            else
               UpperQuantile = quantile(NormXY(ImLin),0.99);
            end
            Maxes = (NormXY(:) > UpperQuantile); 
            ImageQuantiles = quantile(NormXY(logical(ImLin.*not(Maxes))),NoOfQuantiles);
            SemiQuantiles(i,:) =  ImageQuantiles; 
            % for Contrast Stretch
            LowerValue = quantile(NormXY(ImLin),0.1);
            LowerUpperQuantiles(i,:) = [LowerValue UpperQuantile];
            % for upper precentile
            EstimatedMaxTissueIntensitySeriers(i) = UpperQuantile;
            clc
            parfor_progress;      
        end
    else 
        for i = MinIm:MaxIm
            % Load Images
            FilterXYPar = FilterXY;
            Folder16bitInfoPar = Folder16bitInfo;
            Ti = Tiff([DirectoryOf16bit,'' filesep '',Folder16bitInfoPar(FirstImageFileIndex+i-1).name],'r');
            ImD = double(Ti.read());
            Ti.close();
            % support
            Ts = Tiff([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'SupportImages' filesep 'Tissue_Area_', sprintf('%3.4d',i)],'r');
            ImS = double(Ts.read());
            ImSupp = Ts.read();
            Ts.close();           
            % find quantile value for each image individually 
            ImageQuantiles = quantile(ImD(:),NoOfQuantiles);
            if UserDefinedQuantile==1  % enable the user to define a value larger than the brightest pixel
                 EstimatedMaxTissueIntensity = max(ImD(:)); 
            else
                 EstimatedMaxTissueIntensity = ImageQuantiles(UserDefinedQuantile*NoOfQuantiles)
            end       
             % Remove_Peaks and zeros  
             ImD_Smooth = ImD;
             Maxes = (ImD > EstimatedMaxTissueIntensity);
             TissuePixelvalues = ImD(logical(ImS.*(ImD_Smooth<quantile(ImD_Smooth(:),0.9)).*(ImD_Smooth>quantile(ImD_Smooth(:),0.1))));
             if size(TissuePixelvalues(:),1)>100 
                RandValues = emprand(TissuePixelvalues,sum(sum(Maxes)),1); % generate signal-replacing values
                ImD_Smooth(Maxes)= RandValues;   
             end     
            % Smooth along XY
            ImD_Smooth(not(ImS)) = EstimatedMaxTissueIntensity/2;%quantile(ImD(TissuePixelIndices),UserDefinedQuantile);
            ImD_Smooth = Savitzky_Golay(ImD_Smooth,FilterXYPar(1),FilterXYPar(2),1);
            ImD_Smooth = Savitzky_Golay(ImD_Smooth,FilterXYPar(1),FilterXYPar(2),1); % 2nd time
            ImD_Smooth_Norm = ImD_Smooth/max(ImD_Smooth(:));   
            NormXY = (ImD./ImD_Smooth_Norm);
            SaveTiffInCode([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'Norm_',Folder16bitInfo(FirstImageFileIndex+i-1).name],16,'w',NormXY);
            NewBackgroundEstimateImage = NormXY;
            NewMBI = min(NewBackgroundEstimateImage(Maxes));            
            % log Statistics 
            % for upper precentile
            EstimatedMaxTissueIntensitySeriers(i) = EstimatedMaxTissueIntensity;
            % for Semi-Qantile
            ImLin = ImSupp(:)>0;
            if UserDefinedQuantile ~= 1 
               UpperQuantile = NewMBI; 
            else
               UpperQuantile = quantile(NormXY(ImLin),0.99);
            end
            Maxes = (NormXY(:) > UpperQuantile); 
            ImageQuantiles = quantile(NormXY(logical(ImLin.*not(Maxes))),NoOfQuantiles);
            SemiQuantiles(i,:) =  ImageQuantiles; 
            % for Contrast Stretch
            LowerValue = quantile(NormXY(ImLin),0.1);
            LowerUpperQuantiles(i,:) = [LowerValue UpperQuantile];
            % for upper precentile
            EstimatedMaxTissueIntensitySeriers(i) = UpperQuantile;
            clc
            parfor_progress;     
        end
    end
else % Normalize no background detection 
    Dialog('Correcting Background Across The XY plane...',handles);
    parfor_progress(MaxIm+1-MinIm);
    if Threads
        parfor i = MinIm:MaxIm  
            warning('off','all');
            % init
            FilterXYPar = FilterXY;
            Folder16bitInfoPar = Folder16bitInfo;
            Ti = Tiff([DirectoryOf16bit,'' filesep '',Folder16bitInfoPar(FirstImageFileIndex+i-1).name],'r');
            ImD = double(Ti.read());
            Ti.close();     
            % find quantile value for each image individually 
            ImageQuantiles = quantile(ImD(:),NoOfQuantiles);
            if UserDefinedQuantile==1  % enable the user to define a value larger than the brightest pixel
                EstimatedMaxTissueIntensity = MaxTissueIntensity; 
            else
                EstimatedMaxTissueIntensity = ImageQuantiles(UserDefinedQuantile*NoOfQuantiles)
            end
            TissuePixelIndices = (ImD<EstimatedMaxTissueIntensity);   
            % Remove_Peaks and zero pixels     
            ImD_Smooth = ImD;
            Maxes = (ImD > EstimatedMaxTissueIntensity);
            TissuePixelvalues = ImD(logical(TissuePixelIndices.*(ImD_Smooth<quantile(ImD_Smooth(:),0.9)).*(ImD_Smooth>quantile(ImD_Smooth(:),0.1))));
            RandValues = emprand(TissuePixelvalues,sum(sum(Maxes)),1); % generate signal-replacing values
            ImD_Smooth(Maxes)= RandValues;
            % Smooth along XY
            ImD_Smooth(not(TissuePixelIndices)) = EstimatedMaxTissueIntensity/2;%quantile(ImD(TissuePixelIndices),UserDefinedQuantile);
            ImD_Smooth = Savitzky_Golay(ImD_Smooth,FilterXYPar(1),FilterXYPar(2),1);
            ImD_Smooth = Savitzky_Golay(ImD_Smooth,FilterXYPar(1),FilterXYPar(2),1); % 2nd time
            ImD_Smooth_Norm = ImD_Smooth/max(ImD_Smooth(:));
            NormXY = (ImD./ImD_Smooth_Norm);
            SaveTiffInCode([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'Norm_',Folder16bitInfo(FirstImageFileIndex+i-1).name],16,'w',NormXY);
            NewBackgroundEstimateImage = NormXY;
            NewMBI = min(NewBackgroundEstimateImage(Maxes));
            % for Semi-Qantile
            ImSupp = ones(size(ImD,1),size(ImD,2));
            ImLin = ImSupp(:)>0;
            if UserDefinedQuantile ~= 1 
               UpperQuantile = NewMBI; 
            else
               UpperQuantile = quantile(NormXY(ImLin),0.99);
            end
            Maxes = (NormXY(:) > UpperQuantile); 
            ImageQuantiles = quantile(NormXY(logical(ImLin.*not(Maxes))),NoOfQuantiles);
            SemiQuantiles(i,:) =  ImageQuantiles; 
            % for Contrast Stretch
            LowerValue = quantile(NormXY(ImLin),0.1);
            LowerUpperQuantiles(i,:) = [LowerValue UpperQuantile]; 
            % for upper precentile
            EstimatedMaxTissueIntensitySeriers(i) = UpperQuantile;
            clc
            parfor_progress;
        end
    else
        for i = MinIm:MaxIm  
            warning('off','all');
            % init
            FilterXYPar = FilterXY;
            Folder16bitInfoPar = Folder16bitInfo;
            Ti = Tiff([DirectoryOf16bit,'' filesep '',Folder16bitInfoPar(FirstImageFileIndex+i-1).name],'r');
            ImD = double(Ti.read());
            Ti.close(); 
             %     find quantile value for each image individually 
            ImageQuantiles = quantile(ImD(:),NoOfQuantiles);
            if UserDefinedQuantile==1  % enable the user to define a value larger than the brightest pixel
                EstimatedMaxTissueIntensity = MaxTissueIntensity; 
            else
                EstimatedMaxTissueIntensity = ImageQuantiles(UserDefinedQuantile*NoOfQuantiles);
            end
            TissuePixelIndices = (ImD<EstimatedMaxTissueIntensity);
            % Remove_Peaks and zero pixels     
            ImD_Smooth = ImD;
            Maxes = (ImD > EstimatedMaxTissueIntensity);
            TissuePixelvalues = ImD(logical(TissuePixelIndices.*(ImD_Smooth<quantile(ImD_Smooth(:),0.9)).*(ImD_Smooth>quantile(ImD_Smooth(:),0.1))));
            RandValues = emprand(TissuePixelvalues,sum(sum(Maxes)),1); % generate signal-replacing values
            ImD_Smooth(Maxes)= RandValues;
            % Smooth along XY
            ImD_Smooth(not(TissuePixelIndices)) = EstimatedMaxTissueIntensity/2;%quantile(ImD(TissuePixelIndices),UserDefinedQuantile);
            ImD_Smooth = Savitzky_Golay(ImD_Smooth,FilterXYPar(1),FilterXYPar(2),1);
            ImD_Smooth = Savitzky_Golay(ImD_Smooth,FilterXYPar(1),FilterXYPar(2),1); % 2nd time
            ImD_Smooth_Norm = ImD_Smooth/max(ImD_Smooth(:));
            NormXY = (ImD./ImD_Smooth_Norm);
            SaveTiffInCode([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'Norm_',Folder16bitInfo(FirstImageFileIndex+i-1).name],16,'w',NormXY);       
            % now that the image is corrected we estimate what the new MBI is.
            NewBackgroundEstimateImage = NormXY;
            NewMBI = min(NewBackgroundEstimateImage(Maxes));
            % for Semi-Qantile
            ImSupp = ones(size(ImD,1),size(ImD,2));
            ImLin = ImSupp(:)>0;
            if UserDefinedQuantile ~= 1 
                UpperQuantile = NewMBI;
            else
                UpperQuantile = quantile(NormXY(ImLin),0.99);
            end
            Maxes = (NormXY(:) > UpperQuantile); 
            ImageQuantiles = quantile(NormXY(logical(ImLin.*not(Maxes))),NoOfQuantiles);
            SemiQuantiles(i,:) =  ImageQuantiles; 
            % for Contrast Stretch
            LowerValue = quantile(NormXY(ImLin),0.1);
            LowerUpperQuantiles(i,:) = [LowerValue UpperQuantile];
            % for upper precentile
            EstimatedMaxTissueIntensitySeriers(i) = UpperQuantile;
            clc
            parfor_progress;      
        end
    end
end
%% Normalize Images along Z
clc
if NormalizeType == 1 % quantile    
     SemiQuantilesMean = quantile(SemiQuantiles(MinIm:MaxIm,:),0.98);
     Dialog('Quantile Normalization...',handles);
     parfor_progress(MaxIm+1-MinIm); % initiate parallel progress tracking
     if Threads 
         parfor i = MinIm:MaxIm
             warning('off','all');
             Folder16bitInfoPar = Folder16bitInfo;
             ImSizePar = ImSize;
             Ti = Tiff([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'Norm_',Folder16bitInfoPar(FirstImageFileIndex+i-1).name],'r+');
             ImD = double(Ti.read());
             % Convert to double 
             if HasBackground
                 Ts = Tiff([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'SupportImages' filesep 'Tissue_Area_', sprintf('%3.4d',i)],'r');
                 ImS = Ts.read(); 
                 Ts.close();
             else 
                 ImS = ImD>=0;
             end
             ImSLin = ImS(:);
             % Remove_Peaks (Somas)
             if UserDefinedQuantile ~= 1 
                  UpperQuantile = EstimatedMaxTissueIntensitySeriers(i); 
             else
                  UpperQuantile = quantile(ImD(ImSLin),0.99);
             end
             LowerValue = quantile(ImD(ImSLin),0.5);
             Maxes = (ImD > UpperQuantile);
             NormXYZ = ImD;
             NormXYZ(Maxes) = (ImD(Maxes) - LowerValue)*( (SemiQuantilesMean(end) - SemiQuantilesMean(NoOfQuantiles/2)) / (UpperQuantile - LowerValue) ) + SemiQuantilesMean(NoOfQuantiles/2); % contrast sretching of high intensity pixels
             ImSLin = logical(ImSLin.*not(Maxes(:))) ;
             ImD_N = ImD.*double(ImS);
             vq = interp1(1:NoOfQuantiles,SemiQuantilesMean,linspace(1,NoOfQuantiles,sum(ImSLin)));
             [~,OrderdImD] = sort(ImD_N(ImSLin));
             QunatiledImd = ImD_N(ImSLin);
             QunatiledImd(OrderdImD) = vq;
             NormXYZ(ImSLin) = QunatiledImd;
             NormXYZ = reshape(NormXYZ,ImSizePar(1),ImSizePar(2));
             NormXYZ(Maxes) = (max(NormXYZ(ImSLin))/min(NormXYZ(Maxes)))*NormXYZ(Maxes); 
             SaveTiffInCode([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'Norm_',Folder16bitInfo(FirstImageFileIndex+i-1).name],16,'w',NormXYZ);
             parfor_progress;
         end
     else
         for i = MinIm:MaxIm
             warning('off','all');
             Folder16bitInfoPar = Folder16bitInfo;
             ImSizePar = ImSize;
             Ti = Tiff([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'Norm_',Folder16bitInfoPar(FirstImageFileIndex+i-1).name],'r+');
             ImD = double(Ti.read());
             % Convert to double  
             if HasBackground
                 Ts = Tiff([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'SupportImages' filesep 'Tissue_Area_', sprintf('%3.4d',i)],'r');
                 ImS = Ts.read(); 
                 Ts.close();
             else 
                 ImS = ImD>=0;
             end
             ImSLin = ImS(:);
             % Remove_Peaks (Somas)
             if UserDefinedQuantile ~= 1 
                 UpperQuantile = EstimatedMaxTissueIntensitySeriers(i); 
             else
                 UpperQuantile = quantile(ImD(ImSLin),0.99);
             end
             LowerValue = quantile(ImD(ImSLin),0.5);
             Maxes = (ImD > UpperQuantile);
             NormXYZ = ImD;
             NormXYZ(Maxes) = (ImD(Maxes) - LowerValue)*( (SemiQuantilesMean(end) - SemiQuantilesMean(NoOfQuantiles/2)) / (UpperQuantile - LowerValue) ) + SemiQuantilesMean(NoOfQuantiles/2); % contrast sretching of high intensity pixels
             ImSLin = logical(ImSLin.*not(Maxes(:))) ;
             ImD_N = ImD.*double(ImS);
             vq = interp1(1:NoOfQuantiles,SemiQuantilesMean,linspace(1,NoOfQuantiles,sum(ImSLin)));
             [~,OrderdImD] = sort(ImD_N(ImSLin));
             QunatiledImd = ImD_N(ImSLin);
             QunatiledImd(OrderdImD) = vq;
             NormXYZ(ImSLin) = QunatiledImd;
             NormXYZ = reshape(NormXYZ,ImSizePar(1),ImSizePar(2));
             NormXYZ(Maxes) = (max(NormXYZ(ImSLin))/min(NormXYZ(Maxes)))*NormXYZ(Maxes); 
             SaveTiffInCode([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'Norm_',Folder16bitInfo(FirstImageFileIndex+i-1).name],16,'w',NormXYZ);
%              parfor_progress;
         end
     end
%      parfor_progress;
end  
if NormalizeType == 2 % contrast stretch 
   LowerUpperQuantilesMean = quantile(LowerUpperQuantiles(MinIm:MaxIm,:),0.98);
   Dialog('Contrast Stretch Normalization...',handles);
   parfor_progress(MaxIm+1-MinIm); % initiate parallel progress tracking
   if Threads
       parfor i = MinIm:MaxIm
           warning('off','all');
           ImSizePar = ImSize;
           Folder16bitInfoPar = Folder16bitInfo;
           Ti = Tiff([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'Norm_',Folder16bitInfoPar(FirstImageFileIndex+i-1).name],'r+');
           ImD = double(Ti.read());
           % Convert to double 
           NormXYZ = ImD;
           if HasBackground
               Ts = Tiff([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'SupportImages' filesep 'Tissue_Area_', sprintf('%3.4d',i)],'r');
               ImS = Ts.read(); 
               Ts.close();
           else
               ImS = ImD>=0;
           end
           ImSLin = ImS(:);
           if UserDefinedQuantile ~= 1 
               UpperQuantile = EstimatedMaxTissueIntensitySeriers(i); 
           else
               UpperQuantile = quantile(ImD(ImSLin),0.95);
           end
           LowerQuantile = quantile(ImD(ImSLin),0.1);
           NormXYZ(ImSLin) = (ImD(ImSLin) - LowerQuantile)*( (LowerUpperQuantilesMean(2) - LowerUpperQuantilesMean(1)) / (UpperQuantile - LowerQuantile) ) + LowerUpperQuantilesMean(1); % contrast sretching
           SaveTiffInCode([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'Norm_',Folder16bitInfo(FirstImageFileIndex+i-1).name],16,'w',NormXYZ);
           parfor_progress;
       end
   else
       for i = MinIm:MaxIm
           warning('off','all');
           ImSizePar = ImSize;
           Folder16bitInfoPar = Folder16bitInfo;
           Ti = Tiff([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'Norm_',Folder16bitInfoPar(FirstImageFileIndex+i-1).name],'r+');
           ImD = double(Ti.read());
           % Convert to double 
           NormXYZ = ImD;
           if HasBackground
               Ts = Tiff([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'SupportImages' filesep 'Tissue_Area_', sprintf('%3.4d',i)],'r');
               ImS = Ts.read(); 
               Ts.close();
           else 
               ImS = ImD>=0;
           end
           ImSLin = ImS(:);
           if UserDefinedQuantile ~= 1 
               UpperQuantile = EstimatedMaxTissueIntensitySeriers(i); 
           else
               UpperQuantile = quantile(ImD(ImSLin),0.95);
           end
           LowerQuantile = quantile(ImD(ImSLin),0.1);
           NormXYZ(ImSLin) = (ImD(ImSLin) - LowerQuantile)*( (LowerUpperQuantilesMean(2) - LowerUpperQuantilesMean(1)) / (UpperQuantile - LowerQuantile) ) + LowerUpperQuantilesMean(1); % contrast sretching
           SaveTiffInCode([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'Norm_',Folder16bitInfo(FirstImageFileIndex+i-1).name],16,'w',NormXYZ);
           parfor_progress;
        end
     end
end    
 if NormalizeType == 3 % Upper Percentile
     EstimatedMaxTissueIntensitySeriersNorm = EstimatedMaxTissueIntensitySeriers/max(EstimatedMaxTissueIntensitySeriers);
     Dialog('Upper Percentile Normalization...',handles);
     parfor_progress(MaxIm+1-MinIm); % initiate parallel progress tracking
     if Threads
         parfor i = MinIm:MaxIm
        warning('off','all');
        Folder16bitInfoPar = Folder16bitInfo;
        Ti = Tiff([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'Norm_',Folder16bitInfoPar(FirstImageFileIndex+i-1).name],'r+');
        NormXY = double(Ti.read());         
        NormXYZ = (NormXY./EstimatedMaxTissueIntensitySeriersNorm(i)); %Normalize
        SaveTiffInCode([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'Norm_',Folder16bitInfo(FirstImageFileIndex+i-1).name],16,'w',NormXYZ);
        parfor_progress;
         end
     else
        for i = MinIm:MaxIm
        warning('off','all');
        Folder16bitInfoPar = Folder16bitInfo;
        Ti = Tiff([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'Norm_',Folder16bitInfoPar(FirstImageFileIndex+i-1).name],'r+');
        NormXY = double(Ti.read());
        NormXYZ = (NormXY./EstimatedMaxTissueIntensitySeriersNorm(i)); %Normalize
        SaveTiffInCode([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'Norm_',Folder16bitInfo(FirstImageFileIndex+i-1).name],16,'w',NormXYZ);
        parfor_progress;
        end
     end
 end    
set(handles.Execute,'FontSize',20);
set(handles.Execute,'ForegroundColor','black');
set(handles.Execute,'String','Start');
msgbox(['Done!, in ',num2str(toc),' seconds']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% auxillary functions
% find background 
function FindBackground(i,DirectoryOf16bit,HasBackground,ImD,EstimatedMaxTissueIntensity,Inn,BackroundThreshold,ImSizePar,FilterXY,STDNumber)
    SignalPixels = (ImD>EstimatedMaxTissueIntensity);
    ImD(SignalPixels) = EstimatedMaxTissueIntensity; % remove spikes in flourescence
    if HasBackground == 1 % find background with EM
        PCAscore = ImagePCAInCode(ImD,Inn);
        [~,BackgroundPixelIndices] = SeperateTissueBackgroundEMInCode(ImD,PCAscore,STDNumber); 
    end
    if HasBackground == 2 % find background with adaptive K means
        PCAscore = ImagePCAInCode(ImD,Inn);
        silh = [];
        idx = [];
        DoensampledData = PCAscore(1:1000:end,:);
        for k = 2:4
             idx(:,k-1) = kmeans(DoensampledData,k,'Distance','sqEuclidean');
             [silh(:,k-1)] = silhouette(DoensampledData,idx(:,k-1),'sqEuclidean');
        end   
        [~,BestK] = max(mean(abs(silh)));
        BestK = BestK+1;
        BestKInd = cell(BestK,1);
        MeanK = [];
        STDK = [];
        for k = 1:BestK;
            BestKInd{k} = find(idx(:,BestK-1) == k);
            MeanK(k,:) = mean(DoensampledData(BestKInd{k},:));
            STDK(k,:) = std(DoensampledData(BestKInd{k},:));
        end
        [BackgroundMeanValue,BackgroundMean] = min(MeanK(:,1));
        BackgroundPixelIndices = find(PCAscore(:,1) < (BackgroundMeanValue+STDK(BackgroundMean,1)*(19.5-STDNumber*2)));
    end   
    if HasBackground == 3 % bacground brighness     
        BackgroundPixelIndices = find(ImD<BackroundThreshold);
    end     
    SuppImage = ones(ImSizePar(1),ImSizePar(2));
    SuppImage(BackgroundPixelIndices) = 0;
    warning('off','all');
    SuppFilterXY = round(FilterXY(1)/4)+ mod(round(FilterXY(1)/4),2)-1;
    SuppImageSmooth = Savitzky_Golay(SuppImage,SuppFilterXY,SuppFilterXY,1);
    SuppImageSmooth(SuppImageSmooth<0.5) = 0; 
    SuppImageSmooth(SuppImageSmooth>=0.5) = 1;    
    SaveTiffInCode([DirectoryOf16bit,'' filesep 'NormalizedBackground' filesep 'SupportImages' filesep 'Supp_Image_', sprintf('%3.4d',i)],1,'w',logical(SuppImageSmooth));
    clc
end
% simape Threhold
function ThresholdValue = SimpleThrehold(ImD_Downsampled)
            counter = 0;
            PassThresholdpixels = [];
            ThresholdJumps = round(min(ImD_Downsampled(:))):2:round(max(ImD_Downsampled(:)));
            for ii = ThresholdJumps
                counter = counter+1;
                PassThresholdpixels(counter) = sum(ImD_Downsampled(:)>ii);
            end
            GetSlope = movingslope(movingslope(PassThresholdpixels,round(length(ThresholdJumps)/100)),2);
%             plot(ThresholdJumps(PeakIndex(MaxIndex):end),[ChangeInPassThreshold])
            [Value PeakIndex] = findpeaks(GetSlope);
            [MaxValue MaxIndex] = max(Value);
            ChangeInPassThreshold = GetSlope(PeakIndex(MaxIndex):end);
            Candidates = find(islocalmin([ChangeInPassThreshold]));
            
         
            
            ThresholdIndex = Candidates(1)+PeakIndex(MaxIndex);
            ThresholdValue = ThresholdJumps(ThresholdIndex);
end
% PCA 
function PCAscore = ImagePCAInCode(ImDBackground,Inn)
% Filter image  
    ImSd = stdfilt(ImDBackground, ones(3)); % filter
    ImMd = medfilt2(ImDBackground);   
    ImDBackgroundND = ImDBackground(Inn);
    ImMd = ImMd(:);
    ImSd = ImSd(:);
     % PCA for dimentionality reduction and rotation 
     [~,PCAscore,~,~,~,~] = pca([ImDBackgroundND ImMd ImSd] ,'NumComponents',2); % changed from 3 to 2 in version 1_3
     
end  
% Seperate TissueBackground with EM 
function [TissuePixelIndices,BackgroundPixelIndices] = SeperateTissueBackgroundEMInCode(ImD,PCAscore,STDNumber)
    ImDLinear = ImD(:);

    % test for bimodality 
    [xpdf, ~, ~] = compute_xpdf(PCAscore);
    nboot = 2000;
    [~, p_value, ~, ~] = HartigansDipSignifTest(xpdf, nboot); 
    
    % Gaussian Mixture EM
    if p_value < 0.01 %if bimodal -> perform EM
        NoOfGaussians   = 3; % background and tissue and something else
        PointForEM = round(length(ImDLinear)/100); % use only 1% for training 
        jumps = floor(length(PCAscore)/PointForEM);     
        TrainingData = PCAscore(1:jumps:end,:);

        [Labeling,~,~] = mixGaussEmInCode(TrainingData',NoOfGaussians); 
        Cluster1 = TrainingData(Labeling==1,:);
        Cluster2 = TrainingData(Labeling==2,:);
        Cluster3 = TrainingData(Labeling==3,:);
        ClusterMeans = [mean(Cluster1);mean(Cluster2);mean(Cluster3)];
        ClusterSTD = [std(Cluster1);std(Cluster2);std(Cluster3)];
       
        MiniImD =  ImD(:);
        MiniImD = MiniImD(1:jumps:end);
        ClusterIntestinyMean(1) =  mean(MiniImD((Labeling==1)));
        ClusterIntestinyMean(2) =  mean(MiniImD((Labeling==2)));
        ClusterIntestinyMean(3) =  mean(MiniImD((Labeling==3)));
%             ClusterIntestinyMean(4) =  mean(MiniImD((Labeling==4)));
        STDs = (STDNumber)^2/100; 
        [~,SmallSTDClusterIdentity] = min(sum(ClusterSTD,2));
        [~,SmallIntensityClusterIdentity] = min(ClusterIntestinyMean);
        if (SmallSTDClusterIdentity == SmallIntensityClusterIdentity) && (sum(Labeling==SmallIntensityClusterIdentity)>PointForEM/1000)
            BackgroundMean = ClusterMeans(SmallSTDClusterIdentity,:);
            BackgroundSTD = ClusterSTD(SmallSTDClusterIdentity,:);
            InElipseXYZ = ...
              ( ( (PCAscore(:,1)-BackgroundMean(1)).^2 ) / ((BackgroundSTD(1)/(STDs))^2) + ...
                ( (PCAscore(:,2)-BackgroundMean(2)).^2 ) / ((BackgroundSTD(2)/(STDs))^2) ) <= 1;             
            InElipseIndex = find(InElipseXYZ);
            NotInElipseIndex = find(not(InElipseXYZ));
            BackgroundPixelIndices = InElipseIndex;
            TissuePixelIndices = NotInElipseIndex;
        else % not able to find a gaussian distribution that has the minimal intensity background found
            TissuePixelIndices = find(ImDLinear);
            BackgroundPixelIndices = [];
        end
        
    else % not able to descriminate background found
    
        TissuePixelIndices = find(ImDLinear);
        BackgroundPixelIndices = [];
    end
%     
%         figure;hold all; 
%         plot3(PCAscore(1:2000:end,1),PCAscore(1:2000:end,2),PCAscore(1:2000:end,3),'.g')
%         plot3(ClusterMeans(:,1),ClusterMeans(:,2),ClusterMeans(:,3),'ob','MarkerFaceColor','b','MarkerSize',10)
%           plot3(PCAscore(InElipseIndex(1:2000:end),1),PCAscore(InElipseIndex(1:2000:end),2),PCAscore(InElipseIndex(1:2000:end),3),'ob')

%        ReconstructedBackImage = zeros(length(ImDLinear),1);
%        ReconstructedBackImage(BackgroundPixelIndices) = 1;
%        ReconstructedBackImage(TissuePixelIndices) = 2;
%        imagesc(reshape(ReconstructedBackImage,2560,2160));
end
% EM
function [label, model, llh] = mixGaussEmInCode(X, init)
% Perform EM algorithm for fitting the Gaussian mixture model.
% Input:
%   X: d x n data matrix
%   init: k (1 x 1) number of components or label (1 x n, 1<=label(i)<=k) or model structure
% Output:
%   label: 1 x n cluster label
%   model: trained model structure
%   llh: loglikelihood
% Written by Mo Chen (sth4nth@gmail.com).
%% init
% fprintf('EM for Gaussian mixture: running ... \n');
tol = 1e-6;
maxiter = 1000;
llh = -inf(1,maxiter);
R = initialization(X,init);
for iter = 2:maxiter
    [~,label(1,:)] = max(R,[],2);
    R = R(:,unique(label));   % remove empty clusters
    model = maximization(X,R);
    [R, llh(iter)] = expectation(X,model);
    if abs(llh(iter)-llh(iter-1)) < tol*abs(llh(iter)); break; end;
end
llh = llh(2:iter);

function R = initialization(X, init)
n = size(X,2);
if isstruct(init)  % init with a model
    R  = expectation(X,init);
elseif numel(init) == 1  % random init k
    k = init;
    label = ceil(k*rand(1,n));
    R = full(sparse(1:n,label,1,n,k,n));
elseif all(size(init)==[1,n])  % init with labels
    label = init;
    k = max(label);
    R = full(sparse(1:n,label,1,n,k,n));
else
    error('ERROR: init is not valid.');
end
end
function [R, llh] = expectation(X, model)
mu = model.mu;
Sigma = model.Sigma;
w = model.w;

n = size(X,2);
k = size(mu,2);
R = zeros(n,k);
for i = 1:k
    R(:,i) = loggausspdf(X,mu(:,i),Sigma(:,:,i));
end
R = bsxfun(@plus,R,log(w));
T = logsumexp(R,2);
llh = sum(T)/n; % loglikelihood
R = exp(bsxfun(@minus,R,T));
end
function model = maximization(X, R)
[d,n] = size(X);
k = size(R,2);
nk = sum(R,1);
w = nk/n;
mu = bsxfun(@times, X*R, 1./nk);

Sigma = zeros(d,d,k);
r = sqrt(R);
for i = 1:k
    Xo = bsxfun(@minus,X,mu(:,i));
    Xo = bsxfun(@times,Xo,r(:,i)');
    Sigma(:,:,i) = Xo*Xo'/nk(i)+eye(d)*(1e-6);
end

model.mu = mu;
model.Sigma = Sigma;
model.w = w;
end
function y = loggausspdf(X, mu, Sigma)
d = size(X,1);
X = bsxfun(@minus,X,mu);
[U,p]= chol(Sigma);
if p ~= 0
    error('ERROR: Sigma is not PD.');
end
Q = U'\X;
q = dot(Q,Q,1);  % quadratic term (M distance)
c = d*log(2*pi)+2*sum(log(diag(U)));   % normalization constant
y = -(c+q)/2;
end

function s = logsumexp(X, dim)
% Compute log(sum(exp(X),dim)) while avoiding numerical underflow.
%   By default dim = 1 (columns).
% Written by Mo Chen (sth4nth@gmail.com).
if nargin == 1, 
    % Determine which dimension sum will use
    dim = find(size(X)~=1,1);
    if isempty(dim), dim = 1; end
end

% subtract the largest in each dim
y = max(X,[],dim);
s = y+log(sum(exp(bsxfun(@minus,X,y)),dim));   % TODO: use log1p
i = isinf(y);
if any(i(:))
    s(i) = y(i);
end
end


end
% Calculate the linear indices for neighboring points in a matrix 
function [Iadj , Radj, Nfound ] = neighbourND( index, sizeA, res )
% function  [Iadj , Radj, Nfound] = neighbour3D( index,  sizeA, res )
% Calculate the linear indices for neighboring points in a matrix 
% Second output is and array of distances based on an input resolution vector
% This resolution vector defaults to ones(1,ndims)
% The output Nfound reports the number of neighbours found in within the
% matrix. For 2D we expect up to 8, for 3D up to 26 etc...
% 
% Example 1:
% A is a 128x128x16 image data matrix with a spatial resolution of
% 0.1x 0.25x .375 mm^3 
% to get the neighbouring point linear indices for point 456 we do
% sizeA = [128 128 16]
% [ Iadj , Radj, Nfound] = neighbourND( 456, sizeA, [ .10 .25 .375] )
%
% NEW: now index can be a column array with linear indices
% Output Iadj will be Nx8 (2D) or Nx26 (3D) etc and Radj will be 
% a row array 1x8 or 1x26 etc... 
%
% Example 2:
% create points near the center of a 144x192x16 matrix
% spatial resolution .3 x .3x 5 mm^3
% idx = (-6:1:6)+((144*192*3)+144*96+76)
%[ Iadj , Radj, Nfound] = neighbourND( idx , [144,192, 32] , [.3, 0.3, 5])
% Results in 11x26 matrix Iadj, 
% 26 distances in Radj and Nfound is 26
%
% The neighbour indices outside the matrix will be zero!
% when a single index is entered the outside points are still removed so a
% point in a 3D matrix at the edge can sill return 17 neighbours or even less
% when it is a corner.
%==============================================

%==============================================
% Ronald Ouwerkerk 2010 NIH/NIDDK 
% New version: Now handles arrays of indices
% This script is made available on Matlab file exchange by the author 
% for use by other Matlab programmers.
% This script is not intended for commercial use.
% If used for published work a reference or acknowledgement is greatly 
% appreciated.
% The function was tested for several 1D(col and row), 2D, 3D and 4D cases
% I cannot be sure that it really works for all dimensionalities. 
% Let me know if you find a bug (and feel free to squash it for me)
%==============================================

%% Set defaults and process input parameters
% first two are arbitary values for a demo
if nargin <1
    % default index [7,6,2]
    index = 128*128+ 128*5+7;
end

if nargin < 2
    % default size 128x128xN with N big enough for the index
    i3 = floor( index /128/128);
    disp( 'Demo mode')
    sizeA =[128, 128, i3+2];
end

% Get dimensionality
ndimA = length( sizeA );

%Set default resolution to isotropic distances
if nargin < 3
    res =ones(1, length( sizeA) );
else
    if length(res) < ndimA;
        errstr = sprintf('\nError in %s.\n The length of the resolution array (%d) must equal the number of matrix dimensions (%d)\n', ...
                                         mfilename,                                           length(res)  ,                                                         ndimA  );
        disp(errstr)
        help( mfilename)
        return
    else
        % reduce the resolution array, last digit is probably slice
        % thickness, irrelevant if we have one slice only
        res = res( 1:ndimA );
    end
end

%% explicit version of ind2sub 
% ind2sub requires multiple output arguments, one for each dimension
ilin = index(:);
np = length( ilin );
imat = ones( np, ndimA);

for di = ndimA:-1:2    
    blocksize = prod( sizeA( 1:(di-1)  ) );
    ndi = 1+ floor( ( ilin-1) / blocksize );
    ilin = ilin- (ndi -1) *blocksize;
    imat(:,di) = ndi;
end
imat(:,1) = ilin;

%% Find the indices of neighbours
% Get all the index permutations for neighbours ( -1, +1) over all
% dimensions. The total number of neighbours should be three  to the power Ndim
% minus one if we discard the original point itself

% initialize the shift index array
nneighb = 3^ndimA;
nbi = zeros( nneighb, ndimA);

di = ndimA;
while ( di ) 
    N = 3^(di-1);
    ni = 1:N;
    while( ni(end) < nneighb+1 )
        for val=[-1, 0, 1]
              nbi( ni ,di ) = val;
              ni = ni+ N;
        end
    end
    di = di-1;
end

%% Create distance matrix
d = ones(nneighb, 1) * res;
d = d.*abs( nbi );
% create a row vector with distances
dvec = sqrt( sum( d.^2, 2))';
% Get index to exclude the original point: distance = 0
notorig = logical( dvec > 0 );

%% Add the input index array to nbi to get all neighbours
% set up the array for neighbour indices
nd = length( index);
Iadj = zeros( nd, nneighb );
kdo = notorig(ones(nd,1), : ); 

for di = 1:ndimA
    indices = imat( :, di );
    shifts = nbi( :, di )';
    neighbindices = indices( :, ones( 1,nneighb)) +shifts( ones(nd, 1), : ) ;
    maxmat = sizeA( di );
    % set up mask matrix to keep indices within limits and excllude the original point
    s = logical( neighbindices <= maxmat );
    s =logical( neighbindices > 0 ) & s;
    kdo = kdo & s;
    % Calculate the linear index
    if di == 1       
        Iadj( kdo ) =  neighbindices( kdo );
    else
        blocksize = prod( sizeA( 1:(di-1)  ) );
        m = neighbindices-1;
        Iadj(kdo )  = Iadj(kdo )+ m(kdo)*blocksize;
    end
end

%% Select only the sensible points for the neighbour index and distances matrices
% Remove columns that have no valid indices anywhere at all (e.g. origin)
% for shorter index lists with  all points near the edges more may be
% removed.
if nd == 1
    allkdo = any( kdo, 1);
    Iadj = Iadj( :, allkdo);
    Radj = dvec( allkdo );
    Nfound = length(  find( allkdo ) );
else
    Nfound = nneighb-1;
    Radj = dvec;
    iself = (Radj == 0);
    Iadj = Iadj(:,~iself);
    Radj = Radj(~iself);
end
end
% Filter using Savitzky-Golay filtering.
function doublySmoothedImage = Savitzky_Golay(imageArray,WindowV,WindowH,PolynomialOrder)
% Filter using Savitzky-Golay filtering.
% By Image Analyst

% i = 900
% workspace;	% Make sure the workspace panel is showing.
% fontSize = 14;
% Read in standard MATLAB gray scale demo image.
% imageArray = imread('football.jpg');
% imageArray = imread(FolderInfo(i+4).name);
% imageArray = double(imageArray);
[~, ~, ~] = size(imageArray);
% subplot(2, 2, 1);
% imshow(imageArray, [0 255]);
% title('Original Grayscale Image', 'FontSize', fontSize);
% set(gcf, 'Position', get(0,'Screensize')); % Enlarge figure to full screen.
% set(gcf,'name','Savitzky-Golay Filter Demo by ImageAnalyst','numbertitle','off') 
% % Apply the Savitzky-Golay filter.
% First apply it in the vertical (row) direction.
k = PolynomialOrder; % Order of the polynomial
windowSize = WindowV;
verticallySmoothedImage = sgolayfilt(imageArray, k, windowSize, [], 1);
% subplot(2, 2, 2);
% imshow(verticallySmoothedImage, [0 255]);
% title('Savitzky-Golay filtered in the vertical direction only', 'FontSize', fontSize);
% Apply the Savitzky-Golay filter.
% now apply it in the horizintal (row) direction.
% k = PolynomialOrder; % Order of the polynomial
windowSize = WindowH;
%horizontallySmoothedImage = sgolayfilt(imageArray, k, windowSize, [], 2);
% subplot(2, 2, 3);
% imshow(horizontallySmoothedImage, [0 255]);
% title('Savitzky-Golay filtered in the horizontal direction only', 'FontSize', fontSize);
doublySmoothedImage = sgolayfilt(verticallySmoothedImage, k, windowSize, [], 2);
% subplot(2, 2, 4);
% imshow(doublySmoothedImage, [0 255]);
% title('Savitzky-Golay filtered in both directions', 'FontSize', fontSize);
end
% Progress monitor (progress bar) that works with parfor.
function percent = parfor_progress(N)
%PARFOR_PROGRESS Progress monitor (progress bar) that works with parfor.
%   PARFOR_PROGRESS works by creating a file called parfor_progress.txt in
%   your working directory, and then keeping track of the parfor loop's
%   progress within that file. This workaround is necessary because parfor
%   workers cannot communicate with one another so there is no simple way
%   to know which iterations have finished and which haven't.
%
%   PARFOR_PROGRESS(N) initializes the progress monitor for a set of N
%   upcoming calculations.
%
%   PARFOR_PROGRESS updates the progress inside your parfor loop and
%   displays an updated progress bar.
%
%   PARFOR_PROGRESS(0) deletes parfor_progress.txt and finalizes progress
%   bar.
%
%   To suppress output from any of these functions, just ask for a return
%   variable from the function calls, like PERCENT = PARFOR_PROGRESS which
%   returns the percentage of completion.
%
%   Example:
%
%      N = 100;
%      parfor_progress(N);
%      parfor i=1:N
%         pause(rand); % Replace with real code
%         parfor_progress;
%      end
%      parfor_progress(0);
%
%   See also PARFOR.

% By Jeremy Scheff - jdscheff@gmail.com - http://www.jeremyscheff.com/

narginchk(0, 1);

if nargin < 1
    N = -1;
end

percent = 0;
w = 50; % Width of progress bar

if N > 0
    f = fopen('parfor_progress.txt', 'w');
    if f<0
        error('Do you have write permissions for %s?', pwd);
    end
    fprintf(f, '%d\n', N); % Save N at the top of progress.txt
    fclose(f);
    
    if nargout == 0
        disp(['  0%[>', repmat(' ', 1, w), ']']);
    end
elseif N == 0
    delete('parfor_progress.txt');
    percent = 100;
    
    if nargout == 0
        disp([repmat(char(8), 1, (w+9)), char(10), '100%[', repmat('=', 1, w+1), ']']);
    end
else
    if ~exist('parfor_progress.txt', 'file')
        error('parfor_progress.txt not found. Run PARFOR_PROGRESS(N) before PARFOR_PROGRESS to initialize parfor_progress.txt.');
    end
    
    f = fopen('parfor_progress.txt', 'a');
    fprintf(f, '1\n');
    fclose(f);
    
    f = fopen('parfor_progress.txt', 'r');
    progress = fscanf(f, '%d');
    fclose(f);
    percent = (length(progress)-1)/progress(1)*100;
    
    if nargout == 0
        perc = sprintf('%3.0f%%', percent); % 4 characters wide, percentage
        disp([repmat(char(8), 1, (w+9)), char(10), perc, '[', repmat('=', 1, round(percent*w/100)), '>', repmat(' ', 1, w - round(percent*w/100)), ']']);
    end
end
end
% compute_xpdf
function [x2, n, b] = compute_xpdf(x)
  x2 = reshape(x, 1, numel(x));
  [n, b] = hist(x2, 40);
  x2 = sort(x2);
  x2 = interp1 (1:length(x2), x2, 1:1000:length(x2));
end
% calculates Hartigan's DIP statistic and its significance for the empirical p.d.f  XPDF
function [dip, p_value, xlow,xup]=HartigansDipSignifTest(xpdf,nboot)

%  function		[dip,p_value,xlow,xup]=HartigansDipSignifTest(xpdf,nboot)
%
% calculates Hartigan's DIP statistic and its significance for the empirical p.d.f  XPDF (vector of sample values)
% This routine calls the matlab routine 'HartigansDipTest' that actually calculates the DIP
% NBOOT is the user-supplied sample size of boot-strap
% Code by F. Mechler (27 August 2002)

% calculate the DIP statistic from the empirical pdf
[dip,xlow,xup, ~, ~, ~, ~, ~] = HartigansDipTest(xpdf);
N=length(xpdf);

% calculate a bootstrap sample of size NBOOT of the dip statistic for a uniform pdf of sample size N (the same as empirical pdf)
boot_dip=[];
for i=1:nboot
   unifpdfboot=sort(unifrnd(0,1,1,N));
   [unif_dip]=HartigansDipTest(unifpdfboot);
   boot_dip=[boot_dip; unif_dip];
end;
boot_dip=sort(boot_dip);
p_value=sum(dip<boot_dip)/nboot;

% % Plot Boot-strap sample and the DIP statistic of the empirical pdf
% figure(1); clf;
% [hy,hx]=hist(boot_dip); 
% bar(hx,hy,'k'); hold on;
% plot([dip dip],[0 max(hy)*1.1],'r:');
end
% does the dip calculation for an ordered vector XPDF
function	[dip,xl,xu, ifault, gcm, lcm, mn, mj] = HartigansDipTest(xpdf)

% function	[dip,xl,xu, ifault, gcm, lcm, mn, mj]=HartigansDipTest(xpdf)
%
% This is a direct translation by F. Mechler (August 27 2002)
% into MATLAB from the original FORTRAN code of Hartigan's Subroutine DIPTST algorithm 
% Ref: Algorithm AS 217 APPL. STATIST. (1985) Vol. 34. No.3 pg 322-325
%
% Appended by F. Mechler (September 2 2002) to deal with a perfectly unimodal input
% This check the original Hartigan algorithm omitted, which leads to an infinite cycle
%
% HartigansDipTest, like DIPTST, does the dip calculation for an ordered vector XPDF using
% the greatest convex minorant (gcm) and the least concave majorant (lcm),
% skipping through the data using the change points of these distributions.
% It returns the 'DIP' statistic, and 7 more optional results, which include
% the modal interval (XL,XU), ann error flag IFAULT (>0 flags an error)
% as well as the minorant and majorant fits GCM, LCM, and the corresponding support indices MN, and MJ

% sort X in increasing order in column vector
x=sort(xpdf(:));
N=length(x);
mn=zeros(size(x));
mj=zeros(size(x));
lcm=zeros(size(x));
gcm=zeros(size(x));
ifault=0;

% Check that N is positive
if (N<=0) 
   ifault=1;
   fprintf(1,'\nHartigansDipTest.    InputError :  ifault=%d\n',ifault);
   return;
end;

% Check if N is one
if (N==1)
   xl=x(1);
   xu=x(N);
   dip=0.0;
   ifault=2;
   fprintf(1,'\nHartigansDipTest.    InputError :  ifault=%d\n',ifault);
   return;
end;

if (N>1)
   % Check that X is sorted
   if (x ~= sort(x))
      ifault=3;
      fprintf(1,'\nHartigansDipTest.    InputError :  ifault=%d\n',ifault);
      return;
   end;
   % Check for all values of X identical OR for case 1<N<4
   if ~((x(N)>x(1)) && (4<=N))
      xl=x(1);
      xu=x(N);
      dip=0.0;
      ifault=4;
      fprintf(1,'\nHartigansDipTest.    InputError :  ifault=%d\n',ifault);
      return;
   end;
end;

% Check if X is perfectly unimodal
% Hartigan's original DIPTST algorithm did not check for this condition
% and DIPTST runs into infinite cycle for a unimodal input
% The condition that the input is unimodal is equivalent to having 
% at most 1 sign change in the second derivative of the input p.d.f.
xsign=-sign(diff(diff(x)));
% This condition check below works even 
% if the unimodal p.d.f. has its mode in the very first or last point of the input 
% because then the boolean argument is Empty Matrix, and ANY returns 1 for an Empty Matrix
posi=find(xsign>0);
negi=find(xsign<0);
if isempty(posi) || isempty(negi) || all(posi<min(negi))
   % A unimodal function is its own best unimodal approximation, with a zero corresponding dip
   xl=x(1);
   xu=x(N);
   dip=0.0;
   ifault=5;
	%fprintf(1,'\n  The input is a perfectly UNIMODAL input function\n');
   return;
end;

% LOW  contains the index of the current estimate of the lower end of the modal interval
% HIGH contains the index of the current estimate of the upper end of the modal interval
fn=N;
low=1;
high=N;
dip=1/fn;
% xl=x(low);
% xu=x(high);

% establish the indices over which combination is necessary for the convex minorant fit
mn(1)=1;
for j=2:N
   mn(j)=j-1;
   % here is the beginning of a while loop
   mnj=mn(j);
   mnmnj=mn(mnj);
   a=mnj-mnmnj;
   b=j-mnj;
   while ~( (mnj==1) || ((x(j)-x(mnj))*a < (x(mnj)-x(mnmnj))*b))
      mn(j)=mnmnj;
      mnj=mn(j);
      mnmnj=mn(mnj);
      a=mnj-mnmnj;
      b=j-mnj;
   end;   % here is the end of the while loop
end; % end  for j=2:N

% establish the indices over which combination is necessary for the concave majorant fit
mj(N)=N;
na=N-1;
for jk=1:na
   k=N-jk;
   mj(k)=k+1;
   % here is the beginning of a while loop
   mjk=mj(k);
   mjmjk=mj(mjk);
   a=mjk-mjmjk;
   b=k-mjk;
   while ~( (mjk==N) || ((x(k)-x(mjk))*a < (x(mjk)-x(mjmjk))*b))
      mj(k)=mjmjk;
      mjk=mj(k);
      mjmjk=mj(mjk);
      a=mjk-mjmjk;
      b=k-mjk;
   end;   % here is the end of the while loop
end; % end  for jk=1:na

itarate_flag = 1;

% start the cycling of great RECYCLE
while itarate_flag 

% collect the change points for the GCM from HIGH to LOW
% CODE BREAK POINT 40
ic=1;
gcm(1)=high;
igcm1=gcm(ic);
ic=ic+1;
gcm(ic)=mn(igcm1);
while(gcm(ic) > low)
   igcm1=gcm(ic);
   ic=ic+1;
   gcm(ic)=mn(igcm1);
end;
icx=ic;

% collect the change points for the LCM from LOW to HIGH
ic=1;
lcm(1)=low;
lcm1=lcm(ic);
ic=ic+1;
lcm(ic)=mj(lcm1);
while(lcm(ic) < high)
   lcm1=lcm(ic);
   ic=ic+1;
   lcm(ic)=mj(lcm1);
end;
icv=ic;

% ICX, IX, IG are counters for the convex minorant
% ICV, IV, IH are counters for the concave majorant
ig=icx;
ih=icv;

% find the largest distance greater than 'DIP' between the GCM and the LCM from low to high
ix=icx-1;
iv=2;
d=0.0;

% Either GOTO CODE BREAK POINT 65 OR ELSE GOTO CODE BREAK POINT 50;
if ~(icx~=2 || icv~=2)
   d=1.0/fn;
else
   iterate_BP50=1;
   while iterate_BP50
		% CODE BREAK POINT 50
		igcmx=gcm(ix);
      lcmiv=lcm(iv);
      if ~(igcmx > lcmiv)
         % if the next point of either the GCM or LCM is from the LCM then calculate distance here
         % OTHERWISE, GOTO BREAK POINT 55
         lcmiv1=lcm(iv-1);
         a=lcmiv-lcmiv1;
         b=igcmx-lcmiv1-1;
         dx=(x(igcmx)-x(lcmiv1))*a/(fn*(x(lcmiv)-x(lcmiv1)))-b/fn;
         ix=ix-1;
         if(dx < d) 
            goto60 = 1; 
         else
            d=dx;
            ig=ix+1;
            ih=iv;
            goto60 = 1;
         end;
      else
         % if the next point of either the GCM or LCM is from the GCM then calculate distance here
         % CODE BREAK POINT 55
         lcmiv=lcm(iv);
         igcm=gcm(ix);
         igcm1=gcm(ix+1);
         a=lcmiv-igcm1+1;
         b=igcm-igcm1;
         dx=a/fn-((x(lcmiv)-x(igcm1))*b)/(fn*(x(igcm)-x(igcm1)));
         iv=iv+1;
         if ~(dx < d) 
            d=dx;
            ig=ix+1;
            ih=iv-1;
         end;
         goto60 = 1;
      end;
      
      if goto60
         % CODE BREAK POINT 60
         if (ix < 1) 
             ix=1;
         end;
         if (iv > icv) 
             iv=icv; 
         end;
         iterate_BP50 = (gcm(ix) ~= lcm(iv)); 
      end;
   end; % End of WHILE iterate_BP50
end; % End of ELSE (IF ~(icx~=2 | icv~=2)) i.e., either GOTO CODE BREAK POINT 65 OR ELSE GOTO CODE BREAK POINT 50

% CODE BREAK POINT 65
itarate_flag = ~(d < dip);
if itarate_flag
% if itarate_flag is true, then continue calculations and the great iteration cycle
% if itarate_flag is NOT true, then stop calculations here, and break out of great iteration cycle to BREAK POINT 100
   
% calculate the DIPs for the corrent LOW and HIGH

% the DIP for the convex minorant
dl=0.0;
% if not true, go to CODE BREAK POINT 80
if (ig ~= icx)
   icxa=icx-1;
   for j=ig:icxa
      temp=1.0/fn;
   	jb=gcm(j+1);
      je=gcm(j);
      % if not true either, go to CODE BREAK POINT 74
      if ~(je-jb <= 1)
         if~(x(je)==x(jb))
            a=(je-jb);
            const=a/(fn*(x(je)-x(jb)));
            for jr=jb:je
               b=jr-jb+1;
               t=b/fn-(x(jr)-x(jb))*const;
               if (t>temp) 
                   temp=t; 
               end;
            end;
         end;
      end;
      % CODE BREAK POINT 74
      if (dl < temp) 
          dl=temp; 
      end;
   end;
end;

% the DIP for the concave majorant
% CODE BREAK POINT 80
du=0.0;
% if not true, go to CODE BREAK POINT 90
if ~(ih==icv)
   icva=icv-1;
   for k=ih:icva
      temp=1.0/fn;
      kb=lcm(k);
      ke=lcm(k+1);
      % if not true either, go to CODE BREAK POINT 86
      if ~(ke-kb <= 1)
         if ~(x(ke)==x(kb))
            a=ke-kb;
            const=a/(fn*(x(ke)-x(kb)));
            for kr=kb:ke
               b=kr-kb-1;
               t=(x(kr)-x(kb))*const-b/fn;
               if (t>temp) 
                   temp=t; 
               end;
            end;
         end;
      end;
      % CODE BREAK POINT 86
      if (du < temp) 
          du=temp; 
      end;
   end;
end;

% determine the current maximum
% CODE BREAK POINT 90
dipnew=dl;
if (du > dl) 
    dipnew=du; 
end;
if (dip < dipnew) 
    dip=dipnew; 
end;
low=gcm(ig);
high=lcm(ih);      

end; % end of IF(itarate_flag) CODE from BREAK POINT 65

% return to CODE BREAK POINT 40 or break out of great RECYCLE;
end; % end of WHILE of great RECYCLE

% CODE BREAK POINT 100
dip=0.5*dip;
xl=x(low);
xu=x(high);

end
% generate random datasets from empirical distribution
function xr = emprand(dist,varargin)
%EMPRAND Generates random numbers from empirical distribution of data.                   
% This is useful when you do not know the distribution type (i.e. normal or
% uniform), but you have the data and you want to generate random 
% numbers form the data. The idea is to first construct cumulative distribution
% function (cdf) from the given data. Then generate uniform random number and
% interpolate from cdf. 
%
% USAGE:
%         xr = EMPRAND(dist)        - one random number  
%         xr = EMPRAND(dist,m)      - m-by-m random numbers
%         xr = EMPRAND(dist,m,n)    - m-by-n random numbers
%                                             
% INPUT:
%    dist - vector of distribution i.e. data values                                   
%       m - generates m-by-m matrix of random numbers  
%       n - generates m-by-n matrix of random numbers
%       
% OUTPUT:
%    xr - generated random numbers                                                                       
%        
% EXAMPLES:
% % Generate 1000 normal random numbers
% mu = 0; sigma = 1; nr = 1000;
% givenDist = mu + sigma * randn(nr,1);
% generatedDist = emprand(givenDist,nr,1);
% %
% % % Plot histogram to check given and generated distribution
% [n,xout] = hist(givenDist);
% hist(givenDist);
% hold on
% hist(generatedDist,xout)
% %
% Plot cdf to check given and generated distribution
% figure
% x = sort(givenDist(:));      % Given distribution
% p = 1:length(x);
% p = p./length(x);
% plot(x,p,'color','r');      
% hold on
%
% xr = sort(generatedDist(:)); % Generated distribution
% pr = 1:length(xr);
% pr = pr./length(xr);
% 
% plot(xr,pr,'color','b');
% xlabel('x')
% ylabel('cdf')
% legend('Given Dist.','Generated Dist.')
% title('1000 random numbers generated from given normal distribution of data');
% 
% HISTORY:
% version 1.0.0, Release 05-Jul-2005: Initial release
% version 1.1.0, Release 16-Oct-2007: Some bug fixes and improvement of help text
%    1. Can handle NaN values in dist
%    2. Extraplolate for out of range
%    3. Calling function EMPCDF is included within this function
%
% See also: 

% Author: Durga Lal Shrestha
% UNESCO-IHE Institute for Water Education, Delft, The Netherlands
% eMail: durgals@hotmail.com
% Website: http://www.hi.ihe.nl/durgalal/index.htm
% Copyright 2004-2007 Durga Lal Shrestha.
% $First created: 05-Jul-2005
% $Revision: 1.1.0 $ $Date: 16-Oct-2007 21:47:47 $

% ***********************************************************************
%% INPUT ARGUMENTS CHECK

narginchk(1,3);
if ~isvector(dist)
    error('Invalid data size: input data must be vector')
end
if nargin == 2 
    m = varargin{1};
    n = m;
elseif nargin == 3
    m = varargin{1};
    n = varargin{2};
else
    m = 1;
    n = 1;
end

%% COMPUTATION
x = dist(:);
% Remove missing observations indicated by NaN's.
t = ~isnan(x);
x = x(t);

% Compute empirical cumulative distribution function (cdf)
xlen = length(x);
x = sort(x);
p = 1:xlen;
p = p./xlen;   

% Generate uniform random number between 0 and 1
ur =  rand(m,n);

% Interpolate ur from empirical cdf and extraplolate for out of range
% values.
xr = interp1(p,x,ur,[],'extrap');
end

function SaveTiffInCode(filemane,bit,mode,Image)
Ti = Tiff(filemane,mode);
Ti.setTag('ImageLength',size(Image,1));
Ti.setTag('ImageWidth',size(Image,2));
Ti.setTag('Photometric',Tiff.Photometric.MinIsBlack);
Ti.setTag('BitsPerSample',bit);
Ti.setTag('SamplesPerPixel',size(Image,3));
Ti.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
Ti.setTag('Software','MATLAB');
if bit==1; Ti.write(Image); end
if bit==8; Ti.write(uint8(Image)); end
if bit==16; Ti.write(uint16(Image)); end

Ti.close();

end

function Dialog(message,handles) 
    set(handles.Execute,'FontSize',9);
    set(handles.Execute,'ForegroundColor','red');
    set(handles.Execute,'String',message);
    drawnow
end