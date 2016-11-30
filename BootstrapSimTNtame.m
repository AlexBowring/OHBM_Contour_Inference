function BootstrapSimTNboundaries(nSubj,SvNm,nRlz)
%
% function BootstrapSim(nSubj,SvNm,nRlz)
% 
% Sim for investigating the contour infernece method
% 
% Calls SpheroidSignal.m to generate a spheroid of signal 
%
% nSubj = Number of subjects
% SvNm = File name to be saved
% nRlz = Number of toy runs
%----------------------------------------------


%------------Starting Up initialization
if (nargin<1)
  nSubj  = 60;  % Number of subjects
end
if (nargin<2)
  SvNm  = 'Normsim';  % Save name
end
if (nargin<3)
  nRlz = 5000;
end  
if exist([SvNm '.mat'], 'file')
  error('Will not overwrite sim result')
end

%------------Define parameters
%nSubj  = 120;

%nRlz = 1;

Dim    = [112 112 16];     
%Dim    = [84 84 16];     
          %-Spatial extent of signal
Smo    = [1.5 3 6];
                  %-Smoothness in terms of FHWM
Rad = 20; 
          %-Equatorial radius of spheroid
Mag = 3;
          %-Magnitude of signal
rimFWHM = 3;      %-width of padding around the image (xFWHM) to be truncated
                  % to reduce the edge effect

Thr = 1.33; 
          %- Threshold value for bootstrap procedure

nBoot = 500;
          %- Number of bootstrap samples created in bootstrap procedure
          
tau = 1/sqrt(nSubj);

StdBlk = prod(Dim([1 2])/2); % How to block up data when computing Stdev (to
                          % save memory)

%-----------Initialization of Some Variables
Smo    = ones(3,1)*Smo;     %-Smoothness. In this case, isotropic
randn('seed',sum(100*clock));
                            %-Random number generator initializaiton
nSmo   = length(Smo);       %-Number of smoothing kernel widths

wDim    = Dim + 2*ceil(rimFWHM*max(Smo(:)))*ones(1,3);  % Working image dimension


Trunc_x = {(ceil(rimFWHM*max(Smo(:)))+1):(ceil(rimFWHM*max(Smo(:)))+Dim(1))};
Trunc_y = {(ceil(rimFWHM*max(Smo(:)))+1):(ceil(rimFWHM*max(Smo(:)))+Dim(2))};
Trunc_z = {(ceil(rimFWHM*max(Smo(:)))+1):(ceil(rimFWHM*max(Smo(:)))+Dim(3))};

TrnInd = cat(2, Trunc_x, Trunc_y, Trunc_z);
   % Index for truncation

V       = prod(Dim);                          %-Image volume (voxels)

%-All noise for all subjects, allowing an outer smoothing loop
RawNoise = zeros([wDim nSubj]);

%-All signal images for all smoothing kernel widths
Sigdatamat = zeros([Dim nSmo]);                      
AC         = zeros([Dim nSmo]);
% DeltaAC    = zeros([Dim nSmo]); % not presently used

%-All subjects signal + noise images for all smoothing kernel widths
Imgdatamat = zeros([Dim nSubj]);

%- This vector stores the result for each realisation on whether 
% AC^+ < AC < AC^ for each level of smoothing- (1 if true, 0 if false) 
SubsetSuccessVector = zeros(nRlz, nSmo); 
SubsetSuccessVectorout = zeros(nRlz, nSmo); 
SubsetSuccessVectorin = zeros(nRlz, nSmo); 
SubsetSuccessVectormaxinout = zeros(nRlz, nSmo); 

%- This vector stores the threshold value 'c' for each run
ThresholdStore = zeros([nRlz nSmo]); 
ThresholdStorein = zeros([nRlz nSmo]);
ThresholdStoreout = zeros([nRlz nSmo]);
ThresholdStoremaxinout = zeros([nRlz nSmo]);

% This stores the vector SupG for each run
SupGStore = zeros(nBoot, nRlz, nSmo);
SupGoutStore = zeros(nBoot, nRlz, nSmo);
SupGinStore = zeros(nBoot, nRlz, nSmo);
SupGmaxinoutStore = zeros(nBoot, nRlz, nSmo);


%-These matrices store all the sets of interest during the bootstrap
% method for all levels of smoothing

LowerConVolStore = zeros([nRlz nSmo]);
MiddleConVolStore = zeros([nRlz nSmo]);
UpperConVolStore = zeros([nRlz nSmo]);
UpperSubsetMidVolStore = zeros([nRlz nSmo]);
MidSubsetLowerVolStore = zeros([nRlz nSmo]);

LowerConVolStoreout = zeros([nRlz nSmo]);
UpperConVolStoreout = zeros([nRlz nSmo]);
UpperSubsetMidVolStoreout = zeros([nRlz nSmo]);
MidSubsetLowerVolStoreout = zeros([nRlz nSmo]);

LowerConVolStorein = zeros([nRlz nSmo]);
UpperConVolStorein = zeros([nRlz nSmo]);
UpperSubsetMidVolStorein = zeros([nRlz nSmo]);
MidSubsetLowerVolStorein = zeros([nRlz nSmo]);

LowerConVolStoremaxinout = zeros([nRlz nSmo]);
UpperConVolStoremaxinout = zeros([nRlz nSmo]);
UpperSubsetMidVolStoremaxinout = zeros([nRlz nSmo]);
MidSubsetLowerVolStoremaxinout = zeros([nRlz nSmo]);

supG=zeros(nBoot,1);
supGout=zeros(nBoot,1);
supGin=zeros(nBoot,1);

%
% Pre-computations, common to all realisations
%
Sig = SpheroidSignal(wDim, Rad, Mag, 0); %- Signal
% Make the edge image
[a,b,c] = ndgrid(-1:1);
se = strel('arbitrary',sqrt(a.^2 + b.^2 + c.^2) <=1);
for j = 1:nSmo
  % 
  % Smooth signal 
  %
  if Smo(1,j)==0
    Sigs      = Sig;
  else
    Sigs      = zeros(wDim);
    if (length(wDim)>2)
      ss       = spm_smooth(Sig,Sigs,Smo(:,j)');
    else
      [Sigs,ss]   = spm_conv(Sig,Smo(1,j),Smo(2,j));
    end
    Sigs      = Sigs;
  end 
  tSigs    = Sigs(TrnInd{1},TrnInd{2},TrnInd{3});
  MaxtSigs = max(max(max(tSigs)));
  tSigs    = (Mag/MaxtSigs)*tSigs;
  
  Sigdatamat(:,:,:,j) = tSigs; 

  %
  % Obtaining the set AC and boundary DeltaAC
  %
  AC(:,:,:,j) = Sigdatamat(:,:,:,j)>=Thr;
  ACThrDil=imdilate(AC,se);
  ACThrEro=imerode(AC,se);
  %DeltaAC(:,:,:,j)=(ACThrDil - AC)|(AC - ACThrEro); % Not presently used
end      

%
% Realisations loop
%
for t=1:nRlz
  for j= 1:nSmo

    fprintf('Realization %5d FWHM %d mm\n',t,Smo(1,j));
    
    for i=1:nSubj
      %
      % Generate random realizations (once for all smoothings)
      %
      if j==1
        RawNoise(:,:,:,i)   = randn(wDim); %- Noise that will be added to the signal 
      end
      
      %
      % Smooth noise  
      %
      if Smo(1,j)==0
        Noises    = RawNoise(:,:,:,i);
        else
        Noises    = zeros(wDim);
            if (length(wDim)>2)
                tt       = spm_smooth(RawNoise(:,:,:,i),Noises,Smo(:,j)');
                else
                [Noises,tt] = spm_conv(RawNoise(:,:,:,i),Smo(1,j),Smo(2,j));
            end
        Noises    = Noises/sqrt(tt);
      end 
      
      %
      % Truncate to avoid edge effects
      %
      tNoises    = Noises(TrnInd{1},TrnInd{2},TrnInd{3}); 
      
      tImgs = Sigdatamat(:,:,:,j) + tNoises; %- Creates the true image of smoothed signal + smoothed noise
      
      %
      % Writing out the image for quality control
      %
      %spm_unlink(Vqc.fname);
      %Vqc = spm_create_vol(Vqc);
      %Vqc = spm_write_vol(Vqc,tRls);        
      
      %
      % Store smoothed truncated signal image in Sigdatamat and signal + noise image in Imgdatamat for applying bootstrap method 
      %
      
      Imgdatamat(:,:,:,i) = tImgs;
      
    end %========== Loop i (subjects)
    
    %
    % Obtaining the set ACHAT and the boundary DeltaACHAT
    %
    
    ACHATMeanOut = mean(Imgdatamat,4);

    ACHATStdOut = reshape(...
  mystd(reshape(Imgdatamat,[prod(Dim) nSubj]),StdBlk),...
  Dim); 
    
    ACHATTstat = ACHATMeanOut./ACHATStdOut;
    ACHAT = ACHATTstat>=Thr;
    
    % Residuals - just mean centering 
    ACHATResid = bsxfun(@minus,Imgdatamat,ACHATMeanOut);
    ACHATResid = spdiags(1./reshape(ACHATStdOut, [prod(Dim) 1]), 0,prod(Dim),prod(Dim))*reshape(ACHATResid,[prod(Dim) nSubj]);

    % Make the edge image
    ACHATTstatThrDil=imdilate(ACHAT,se);
    ACHATTstatThrEro=imerode(ACHAT,se);
    DeltaACHAT=(ACHATTstatThrDil - ACHAT)|(ACHAT - ACHATTstatThrEro);
    DeltaACHATout =logical(ACHATTstatThrDil - ACHAT);
    DeltaACHATin  =logical(ACHAT - ACHATTstatThrEro);
    
    a=whos;fprintf('Memory usage: %f GB\n',sum([a.bytes])/1024^3)
    
    
    %
    % Apply bootstrapping procedure and finding coverage probability using the true boundary delta A_c
    %
    for k=1:nBoot
    %fprintf('k=%d ',k)
    Ystar = Imgdatamat(:,:,:,randi(nSubj,[nSubj 1]));
    Ystar = bsxfun(@minus,Ystar,ACHATMeanOut);
    Ystar = spdiags(1./reshape(ACHATStdOut, [prod(Dim) 1]), 0,prod(Dim),prod(Dim))*reshape(Ystar,[prod(Dim) nSubj]);
    Ystar = reshape(Ystar, [Dim nSubj]);
  
    MuStar = sum(Ystar,4)/sqrt(nSubj);
      
    supG(k)=max(abs(MuStar(DeltaACHAT)));
      
    supGout(k)=max(abs(MuStar(DeltaACHATout)));
    supGin(k)=max(abs(MuStar(DeltaACHATin)));
      
    %------- Note: The method needs to be applied using DeltaACHAT and DeltaAC
    end
    
    supGa90=prctile(supG,90)
    
    LowerCon = ACHATTstat>=Thr-supGa90*tau;
    MiddleCon = AC(:,:,:,j);
    UpperCon = ACHATTstat>=Thr+supGa90*tau;
    
    
    MidOnUpperMAT = UpperCon.*MiddleCon;
    LowerOnMidMAT = MiddleCon.*LowerCon;
    
    UpperSubsetMidMAT = UpperCon - MidOnUpperMAT;
    MidSubsetLowerMAT = MiddleCon - LowerOnMidMAT;


    SupGStore(:,t,j) = supG; 
    ThresholdStore(t,j) = supGa90;
    
    LowerConVolStore(t,j) = sum(LowerCon(:));
    MiddleConVolStore(t,j) = sum(MiddleCon(:));
    UpperConVolStore(t,j) = sum(UpperCon(:));
    UpperSubsetMidVolStore(t,j) = sum(UpperSubsetMidMAT(:));
    MidSubsetLowerVolStore(t,j) = sum(MidSubsetLowerMAT(:)); %- Saving volumes of sets of interest
    
    if UpperSubsetMidMAT+MidSubsetLowerMAT==zeros(Dim)
      SubsetSuccessVector(t,j) = 1; 
      fprintf('success! \n');
    else 
      SubsetSuccessVector(t,j) = 0; 
      fprintf('failure! \n');
    end

    % Finding the results using the outer, dilated boundary
    supGa90out=prctile(supGout,90)
    
    LowerConout = ACHATTstat>=Thr-supGa90out*tau;
    UpperConout = ACHATTstat>=Thr+supGa90out*tau;
    
    
    MidOnUpperMATout = UpperConout.*MiddleCon;
    LowerOnMidMATout = MiddleCon.*LowerConout;
    
    UpperSubsetMidMATout = UpperConout - MidOnUpperMATout;
    MidSubsetLowerMATout = MiddleCon - LowerOnMidMATout;


    SupGoutStore(:,t,j) = supGout; 
    ThresholdStoreout(t,j) = supGa90out;

    LowerConVolStoreout(t,j) = sum(LowerConout(:));
    UpperConVolStoreout(t,j) = sum(UpperConout(:));
    UpperSubsetMidVolStoreout(t,j) = sum(UpperSubsetMidMATout(:));
    MidSubsetLowerVolStoreout(t,j) = sum(MidSubsetLowerMATout(:)); %- Saving volumes of sets of interest
    
    if UpperSubsetMidMATout+MidSubsetLowerMATout==zeros(Dim)
      SubsetSuccessVectorout(t,j) = 1; 
      fprintf('success! \n');
    else 
      SubsetSuccessVectorout(t,j) = 0; 
      fprintf('failure! \n');
    end

    % Finding the results using the inner, eroded boundary
    supGa90in=prctile(supGin,90)
    
    LowerConin = ACHATTstat>=Thr-supGa90in*tau;
    UpperConin = ACHATTstat>=Thr+supGa90in*tau;
    
    
    MidOnUpperMATin = UpperConin.*MiddleCon;
    LowerOnMidMATin = MiddleCon.*LowerConin;
    
    UpperSubsetMidMATin = UpperConin - MidOnUpperMATin;
    MidSubsetLowerMATin = MiddleCon - LowerOnMidMATin;


    SupGinStore(:,t,j) = supGin; 
    ThresholdStorein(t,j) = supGa90in;

    LowerConVolStorein(t,j) = sum(LowerConin(:));
    UpperConVolStorein(t,j) = sum(UpperConin(:));
    UpperSubsetMidVolStorein(t,j) = sum(UpperSubsetMidMATin(:));
    MidSubsetLowerVolStorein(t,j) = sum(MidSubsetLowerMATin(:)); %- Saving volumes of sets of interest
    
    if UpperSubsetMidMATin+MidSubsetLowerMATin==zeros(Dim)
      SubsetSuccessVectorin(t,j) = 1; 
      fprintf('success! \n');
    else 
      SubsetSuccessVectorin(t,j) = 0; 
      fprintf('failure! \n');
    end       
    
    % Finding the results using the max of supG calculated on the inner and
    % out boundary
    supGa90maxinout=max([supGa90in supGa90out]);
    
    LowerConmaxinout = ACHATTstat>=Thr-supGa90maxinout*tau;
    UpperConmaxinout = ACHATTstat>=Thr+supGa90maxinout*tau;
    
    
    MidOnUpperMATmaxinout = UpperConmaxinout.*MiddleCon;
    LowerOnMidMATmaxinout = MiddleCon.*LowerConmaxinout;
    
    UpperSubsetMidMATmaxinout = UpperConmaxinout - MidOnUpperMATmaxinout;
    MidSubsetLowerMATmaxinout = MiddleCon - LowerOnMidMATmaxinout;

    ThresholdStoremaxinout(t,j) = supGa90maxinout;

    LowerConVolStoremaxinout(t,j) = sum(LowerConmaxinout(:));
    UpperConVolStoremaxinout(t,j) = sum(UpperConmaxinout(:));
    UpperSubsetMidVolStoremaxinout(t,j) = sum(UpperSubsetMidMATmaxinout(:));
    MidSubsetLowerVolStoremaxinout(t,j) = sum(MidSubsetLowerMATmaxinout(:)); %- Saving volumes of sets of interest
    
    if UpperSubsetMidMATmaxinout+MidSubsetLowerMATmaxinout==zeros(Dim)
      SubsetSuccessVectormaxinout(t,j) = 1; 
      fprintf('success! \n');
    else 
      SubsetSuccessVectormaxinout(t,j) = 0; 
      fprintf('failure! \n');
    end    
    
  end  %-------- Loop j (Different levels of smoothing)
end %------ Loop t (nRlz)

PercentageSuccessVector = mean(SubsetSuccessVector, 1);
PercentageSuccessVectorout = mean(SubsetSuccessVectorout, 1);
PercentageSuccessVectorin = mean(SubsetSuccessVectorin, 1);
PercentageSuccessVectormaxinout = mean(SubsetSuccessVectormaxinout, 1);


for z = 1:nSmo
  fprintf('For %d toy runs, using data smoothed with a %d mm FWHM, For %d percent of trials the inclusion AC^+ < AC < AC^- was obtained\n', nRlz, Smo(1,z), PercentageSuccessVector(1,z)*100);
  fprintf('For %d toy runs, evaluted on the OUTER boundary, using data smoothed with a %d mm FWHM, For %d percent of trials the inclusion AC^+ < AC < AC^- was obtained\n', nRlz, Smo(1,z), PercentageSuccessVectorout(1,z)*100);
  fprintf('For %d toy runs, evaluted on the INNER boundary, using data smoothed with a %d mm FWHM, For %d percent of trials the inclusion AC^+ < AC < AC^- was obtained\n', nRlz, Smo(1,z), PercentageSuccessVectorin(1,z)*100);
  fprintf('For %d toy runs, evaluted using the max supG over the INNER and OUTER bounday, using data smoothed with a %d mm FWHM, For %d percent of trials the inclusion AC^+ < AC < AC^- was obtained\n', nRlz, Smo(1,z), PercentageSuccessVectormaxinout(1,z)*100);
end

eval(['save ' SvNm ' nSubj nRlz Dim wDim Smo Rad Mag rimFWHM Thr nBoot ThresholdStore ThresholdStorein ThresholdStoreout ThresholdStoremaxinout LowerConVolStore LowerConVolStorein LowerConVolStoreout LowerConVolStoremaxinout MiddleConVolStore UpperConVolStore UpperConVolStorein UpperConVolStoreout UpperConVolStoremaxinout UpperSubsetMidVolStore UpperSubsetMidVolStorein UpperSubsetMidVolStoreout UpperSubsetMidVolStoremaxinout MidSubsetLowerVolStore MidSubsetLowerVolStorein MidSubsetLowerVolStoreout MidSubsetLowerVolStoremaxinout SubsetSuccessVector SubsetSuccessVectorout SubsetSuccessVectorin SubsetSuccessVectormaxinout PercentageSuccessVector PercentageSuccessVectorin PercentageSuccessVectorout PercentageSuccessVectormaxinout SupGStore SupGinStore SupGoutStore']);


function y = mystd(x,blk)
% Hopefully faster, more memory efficient std - Only works for 2D arrays,
% works *across* rows (not down columns)
n = size(x,2);
nRow=size(x,1);
nBlk=nRow/blk;
y = zeros(nRow,1);
I0 = 1:blk;
for i=1:nBlk
  I = I0+(i-1)*blk;
  xbar = sum(x(I,:), 2) ./ n;
  xc = bsxfun(@minus, x(I,:), xbar);
  y(I) = sqrt(sum(xc.^2, 2) ./ (n-1));
end

return
