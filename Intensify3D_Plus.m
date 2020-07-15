
function Intensify3D_Plus(Folder,NSTD,MSFS,MBI,HasBackground,NormalizeType,Threads,ImageToMeasure,FolderToMeasure,handles,UserSelectedThrehold,TissueSmoothing,ExistingSupportFilesFolder,XYNorm,Between_Stack_Normalization)
% for normalizaing your entire imaging exp in 3D and across different biological samples
% to operate see User_GUI_Intensify3DPlus_v_1_0
% all files must be in tif format and of the same xy size 
% Folder - Parent folder with subdirectories
% NSTD - Number of standard deviations to use for automatic threshold detection
% MSFS - Spatial filter size
% MBI - Maximum Background Intensity
% HasBackground - does the image have imaging media background
% NormalizeType - in-sample normalization type
% Threads - how many cores to use
% ImageToMeasure - Representative image to measure general parameters all
% images will be normalized according to the intensity quantiles of this
% image
% FolderToMeasure - Folder containing "ImageToMeasure"
% handles - Graphical handles 
% UserSelectedThrehold - manual threshold for background detection
% TissueSmoothing - Smoothing factor for background detection
% ExistingSupportFilesFolder - folder that already contains good support files 
% XYNorm - Do xy normalization or only Z normalzation 
% Between_Stack_Normalization - type of normalization to use between
% samples

%% initiation
clc
% get file folders 
MotherDir = Folder;
Directories = dir(MotherDir); 
Directories = Directories(3:end,:);
Directories = Directories([Directories.isdir],:);
NoOfQuantiles = 10000;

%% initiation of parallel pool 
warning('off');
% dialog('Getting Started....',handles);
if Threads
    if ~isempty(gcp('nocreate'))
         delete(gcp); 
    else
        parpool(Threads);
    end   
end


%% delete old files
for d = 1:length(Directories)
    cd(MotherDir);

    cd(Directories(d).name)
    Folder16bitInfo = dir('*.tif');
    flag = size(dir('NormalizedBackground'));
    if ~flag(1);
        mkdir('NormalizedBackground'); 
    else
        cd('NormalizedBackground');
        toDelete = dir;
        for k = 1 : length(toDelete)
            FileName = toDelete(k).name;
            delete(FileName);
        end
        flag = size(dir('SupportImages'));
        if ~flag(1) && HasBackground
            mkdir('SupportImages');
        else
            if HasBackground && ~ExistingSupportFilesFolder
                 cd('SupportImages');
                 toDelete = dir;
                 for k = 3 : length(toDelete)
                    FileName = toDelete(k).name;
                    delete(FileName);
                 end
            end
        end    
    end
end

%% get image size and image quantiles
cd(MotherDir);
cd(Directories(FolderToMeasure).name)
Folder16bitInfo = dir('*.tif');
FirstImageFileIndex = find([Folder16bitInfo.bytes]>10*1000,1,'first');
Ti = Tiff([Folder16bitInfo(ImageToMeasure).name],'r');
ImMd = Ti.read();
ImD = double(ImMd);
ImSize = size(ImD);
Ti.close();
FilterXY = [MSFS MSFS];        
quantiles = quantile(ImD(:),NoOfQuantiles);
largerthan = find(quantiles>MBI);
if ~isempty(largerthan)
    UserSelectedQuantile = largerthan(1)/NoOfQuantiles;
else
    UserSelectedQuantile = 1;
end

%% Normalization of individual samples one by one

SemiQuantilesMean = cell(length(Directories),1);
LowerUpperQuantilesMean = cell(length(Directories),1);
EstimatedMaxTissueIntensitySeriers = cell(length(Directories),1);
if ~ExistingSupportFilesFolder
 ExistingSupportFilesFolders = mat2cell(zeros(length(Directories),1),ones(length(Directories),1));
else
    for d = 1:length(Directories)
        ExistingSupportFilesFolders{d} = [Directories(d).folder, filesep,Directories(d).name, filesep,'NormalizedBackground', filesep,'SupportImages'];
    end
end

handles = 0;
cd(MotherDir)

   
    f = waitbar(0,'Normalizing each sample to itself..... ');
%%
for d = 1:length(Directories)
    pause(1)
    tic
        d
        waitbar(d/length(Directories),f);
        Folder = [MotherDir,'\',Directories(d).name];
        cd(Folder)
        FolderInfo = dir('*.tif');
        FirstImageFileIndex = find([FolderInfo.bytes]>10*1000,1,'first');
        LastImageFileIndex = find([FolderInfo.bytes]>10*1000,1,'last');
        I = imread([Folder,'\',FolderInfo(FirstImageFileIndex+ImageToMeasure-1).name]);
        MaxBsckgroundGuess = quantile(I(:),UserSelectedQuantile);
        MaxIm(d) = LastImageFileIndex-FirstImageFileIndex+1;
        delete([Folder,'\workspace.mat']); 
        MinIm = 1;
        [NewLowerUpperQuantiles{d},NewSemiQuantiles{d}]= Intensify3D_Core(Folder,MinIm,MaxIm(d)...
            ,NSTD,MaxBsckgroundGuess,MSFS,HasBackground,NormalizeType,Threads,ImageToMeasure,handles,...
            UserSelectedThrehold,TissueSmoothing,ExistingSupportFilesFolders{d},XYNorm) ;    
    toc
end

%% Get Image statistics from all samples post normalization
AllSemiQuantilesMean = [];
AllLowerUpperQuantilesMean = [];
AllEstimatedMaxTissueIntensitySeriers = [];
for d = 1:length(Directories) 
    AllSemiQuantilesMean = [AllSemiQuantilesMean; NewSemiQuantiles{d}]; 
    AllLowerUpperQuantilesMean = [AllLowerUpperQuantilesMean; NewLowerUpperQuantiles{d}]; 
    AllEstimatedMaxTissueIntensitySeriers = [AllEstimatedMaxTissueIntensitySeriers; NewLowerUpperQuantiles{d}(2)];     
end
MedianAllSemiQuantilesMean = median(AllSemiQuantilesMean);
MedianAllLowerUpperQuantilesMean = median(AllLowerUpperQuantilesMean);
NormAllEstimatedMaxTissueIntensitySeriers = AllEstimatedMaxTissueIntensitySeriers./max(AllEstimatedMaxTissueIntensitySeriers);
cd(MotherDir)
save workspace % save parameters

%% Normalize all stacks according to same histogram
if Between_Stack_Normalization == 2 % contrast stretch 
    f = waitbar(0,'Normalizing across samples with constrast stretch');
    for d = 1:length(Directories)
        waitbar(d/length(Directories),f);
        DirectoryOf16bit = [MotherDir,'\',Directories(d).name,'\NormalizedBackground'];
        cd(DirectoryOf16bit);
        mkdir('Norm2All')
        Folder16bitInfo = dir('*.tif');

       LowerUpperQuantilesMean = [MedianAllLowerUpperQuantilesMean(1) MedianAllLowerUpperQuantilesMean(2)];
        LowerQuantileStack = AllLowerUpperQuantilesMean(d,1);
        UpperQuantileStack = AllLowerUpperQuantilesMean(d,2); 
           for i = 1:length(Folder16bitInfo)
               warning('off','all');
               Ti = Tiff([Folder16bitInfo(i).folder,'' filesep',Folder16bitInfo(i).name],'r+');
               ImD = double(Ti.read());
               NormXYZ = ImD;
               if HasBackground
                   Ts = Tiff([Folder16bitInfo(i).folder, filesep 'SupportImages' filesep 'Tissue_Area_', sprintf('%3.4d',i)],'r');
                     ImS = Ts.read(); 
                     Ts.close();
               else
                   ImS = ImD>=0;
               end
               ImSLin = logical(ImS(:));  
               NormXYZ(ImSLin) = (ImD(ImSLin) - LowerQuantileStack)*( (LowerUpperQuantilesMean(2) - LowerUpperQuantilesMean(1)) / (UpperQuantileStack - LowerQuantileStack) ) + LowerUpperQuantilesMean(1); % contrast sretching          
               SaveTiffInCode([Folder16bitInfo(i).folder,filesep,'Norm2All' filesep,'All_',Folder16bitInfo(i).name],16,'w',NormXYZ);
           end

        d
    end   

end
%%
    clc
if Between_Stack_Normalization == 1 % quantile    
     SemiQuantilesMean = MedianAllSemiQuantilesMean;
     EstimatedMaxTissueIntensitySeriersCurrent = EstimatedMaxTissueIntensitySeriers{d};
     
         for i = MinIm:MaxIm
             warning('off','all');
             Folder16bitInfoPar = Folder16bitInfo;
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
             
             ImSizePar = size(ImS);
             ImSLin = ImS(:);
             % Remove_Peaks (Somas)
             if UserDefinedQuantile ~= 1 
                 UpperQuantile = EstimatedMaxTissueIntensitySeriers{i}; 
             else
                 UpperQuantile = quantile(ImD(ImSLin),0.99);
             end
             LowerValue = quantile(ImD(ImSLin),0.5);
             Maxes = (ImD > UpperQuantile);
             NormXYZ = ImD;
             NormXYZ(Maxes) = (ImD(Maxes) - LowerValue)*( (SemiQuantilesMean(end) - SemiQuantilesMean(NoOfQuantiles/2)) / (UpperQuantile - LowerValue) ) + SemiQuantilesMean(NoOfQuantiles/2); % contrast sretching of high intensity pixels
             ImSLin = logical(ImSLin.*not(Maxes(:))) ;
             ImD_N = ImD(:).*double(ImSLin);
             vq = interp1(1:NoOfQuantiles,SemiQuantilesMean,linspace(1,NoOfQuantiles,sum(ImSLin)));
             [~,OrderdImD] = sort(ImD_N(ImSLin));
             QunatiledImd = ImD_N(ImSLin);
             QunatiledImd(OrderdImD) = vq;
             NormXYZ(ImSLin) = QunatiledImd;
             NormXYZ = reshape(NormXYZ,ImSizePar(1),ImSizePar(2));
             NormXYZ(Maxes) = (max(NormXYZ(ImSLin))/min(NormXYZ(Maxes)))*NormXYZ(Maxes); 
               SaveTiffInCode([Folder16bitInfo(i).folder,filesep,'Norm2All' filesep,'All_',Folder16bitInfo(i).name],16,'w',NormXYZ);
         end
     
end  
%%
 if Between_Stack_Normalization == 3 % Upper Percentile
     EstimatedMaxTissueIntensitySeriersNorm = NormAllEstimatedMaxTissueIntensitySeriers(d);
       DirectoryOf16bit = [MotherDir,'\',Directories(d).name,'\NormalizedBackground'];
        cd(DirectoryOf16bit);
        mkdir('Norm2All')
        Folder16bitInfo = dir('*.tif');
        for i = MinIm:MaxIm
        warning('off','all');
        Folder16bitInfoPar = Folder16bitInfo;
        Ts = Tiff([Folder16bitInfo(i).folder, filesep 'SupportImages' filesep 'Tissue_Area_', sprintf('%3.4d',i)],'r');
        NormXY = double(Ti.read());
        NormXYZ = (NormXY./EstimatedMaxTissueIntensitySeriersNorm(i)); %Normalize
        SaveTiffInCode([Folder16bitInfo(i).folder,filesep,'Norm2All' filesep,'All_',Folder16bitInfo(i).name],16,'w',NormXYZ);
        end
     
 end    
close(f)
end

%%%%%%%%%%%%%%%%%%%%%%%%

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


