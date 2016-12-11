function [] = RawBrainBootstrappingFSL(String,Out)

cd(String);
[x,y,z] = ndgrid(-1:1);
se = strel('arbitrary',sqrt(x.^2 + y.^2 + z.^2) <=1);

Thr=50;  % In raw change units, mu
nBoot = 1000;

VY=spm_vol('copes.nii.gz');  % This is the file "handle" for all input
                             % images  - Ignore gzip warning
VM=spm_vol('mask.nii.gz'); % This is the handle for the mask

Mask=spm_read_vols(VM)>0;

nSubj=length(VY);
Dim=VY(1).dim;
datamat=zeros([Dim nSubj]);

tau = 1/sqrt(nSubj);

% Load up each of the 3D images into a massive 4D image
for i=1:nSubj
  % Take care to mask each image as its read
  datamat(:,:,:,i) = spm_read_vols(VY(i)).*Mask;
end

% No 2D reshaping (except with Resid)

Mean = mean(datamat,4);
Std = std(datamat,0,4); 

Mu = zeros(Dim);

Mu(Mask) = Mean(Mask);
MuThr = Mu >= Thr;

% Residuals - just mean centering 
Resid = datamat - repmat(Mean,[1 1 1 nSubj]);
Resid = reshape(Resid,[prod(Dim) nSubj]);
Resid(Mask,:) = Resid(Mask,:)./repmat(Std(Mask),[1 nSubj]);

% Make the edge image
MuThrDil=imdilate(MuThr,se);
MuThrEro=imerode(MuThr,se);
MuThrEdge=(MuThrDil-MuThr)|(MuThr-MuThrEro);
MuThrEdge=MuThrEdge & Mask; % Not necessary since there's no data
                                    % outside the mask, but just to be safe

for i=1:nBoot

  SignFlips=randi(2,[nSubj,1])*2-3; 
  % This impliments Eqn (6), the Wild Bootstrap mean 
  Ystar = reshape(reshape(Resid,[prod(Dim) nSubj]).*repmat(SignFlips',[prod(Dim) 1]), [Dim nSubj]);
  MuStar = sum(Ystar,4)/sqrt(nSubj);

  supG(i)=max(abs(MuStar(MuThrEdge)));
     
end

supGa95=prctile(supG,95)

LowerCon = Mu>=Thr-supGa95*tau*Std;
MiddleCon = MuThr;
UpperCon = Mu>=Thr+supGa95*tau*Std;


subplot(2,3,1)
imagesc(Mu(:,:,40));axis image; colorbar
subplot(2,3,2)
hist(supG,50)
abline('v',supGa95)
subplot(2,3,3)
imagesc(MuThrEdge(:,:,40));axis image; colorbar
subplot(2,3,4)
imagesc(Mu(:,:,40)>=Thr-supGa95*tau);axis image; colorbar
subplot(2,3,5)
imagesc(MuThr(:,:,40));axis image; colorbar
subplot(2,3,6)
imagesc(Mu(:,:,40)>=Thr+supGa95*tau);axis image; colorbar

cd(Out);
Vout=VY(1); % clone the first image's handle
Vout.fname = 'RawLowerConfidenceIntervalLowThresh.nii'; % crucially, change the file name!
Vout.descrip = 'Lower confidence interval!'; % Actually, put something more
                                        % informative here

Vout=spm_write_vol(Vout,LowerCon);

Vout=VY(1); % clone the first image's handle
Vout.fname = 'RawMiddleConfidenceIntervalLowThresh.nii'; % crucially, change the file name!
Vout.descrip = 'Middle confidence interval!'; % Actually, put something more
                                        % informative here

Vout=spm_write_vol(Vout,MiddleCon);

Vout=VY(1); % clone the first image's handle
Vout.fname = 'RawUpperConfidenceIntervalLowThresh.nii'; % crucially, change the file name!
Vout.descrip = 'Upper confidence interval!'; % Actually, put something more
                                        % informative here

Vout=spm_write_vol(Vout,UpperCon);

end