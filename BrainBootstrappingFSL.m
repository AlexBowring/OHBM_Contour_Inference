function [] = BrainBootstrappingFSL(String,Out)

cd(String);
[x,y,z] = ndgrid(-1:1);
se = strel('arbitrary',sqrt(x.^2 + y.^2 + z.^2) <=1);

Thr=1.0;  % In standardized effect units, mu/sigma  (not T)
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
Std = reshape(...
  mystd(reshape(datamat,[prod(Dim) nSubj]),StdBlk),...
  Dim); 

MuStdz = zeros(Dim);

MuStdz(Mask) = Mean(Mask)./Std(Mask);
MuStdzThr = MuStdz >= Thr;
fprintf('Standardized threshold %f used; equivalent to T threshold of %f\n',Thr,Thr*sqrt(nSubj));

% Residuals - just mean centering 
Resid = datamat - repmat(Mean,[1 1 1 nSubj]);
Resid = bsxfun(@minus, datamat, Mean);
Resid = reshape(Resid,[prod(Dim) nSubj]);
Resid(Mask,:) = spdiags(1./reshape(Std, [prod(Dim) 1]), 0,prod(Dim),prod(Dim))*reshape(Resid(Mask,:),[prod(Dim) nSubj]);

% Make the edge image
MuStdzThrDil=imdilate(MuStdzThr,se);
MuStdzThrEro=imerode(MuStdzThr,se);
MuStdzThrEdge=(MuStdzThrDil-MuStdzThr)|(MuStdzThr-MuStdzThrEro);
MuStdzThrEdge=MuStdzThrEdge & Mask; % Not necessary since there's no data
                                    % outside the mask, but just to be safe

for i=1:nBoot
  % This impliments the bootstrap method 
  Ystar = datamat(:,:,:,randi(nSubj,[nSubj 1]));
  Ystar = bsxfun(@minus,Ystar,Mean);
  Ystar = spdiags(1./reshape(Std, [prod(Dim) 1]), 0,prod(Dim),prod(Dim))*reshape(Ystar,[prod(Dim) nSubj]);
  Ystar = reshape(Ystar, [Dim nSubj]);
  MuStdzStar = sum(Ystar,4)/sqrt(nSubj);

  supG(i)= max(abs(MuStdzStar(MuStdzThrEdge)))      
end

supGa95=prctile(supG,95)

LowerCon = MuStdz>=Thr-supGa95*tau;
MiddleCon = MuStdzThr;
UpperCon = MuStdz>=Thr+supGa95*tau;

cd(Out);
Vout=VY(1); % clone the first image's handle
Vout.fname = 'StandardizedLowerConfidenceInterval.nii'; % crucially, change the file name!
Vout.descrip = 'Standardized Lower confidence interval!'; % Actually, put something more
                                        % informative here

Vout=spm_write_vol(Vout,LowerCon);

Vout=VY(1); % clone the first image's handle
Vout.fname = 'StandardizedMiddleConfidenceInterval.nii'; % crucially, change the file name!
Vout.descrip = 'Standardized Middle confidence interval!'; % Actually, put something more
                                        % informative here

Vout=spm_write_vol(Vout,MiddleCon);

Vout=VY(1); % clone the first image's handle
Vout.fname = 'StandardizedUpperConfidenceInterval.nii'; % crucially, change the file name!
Vout.descrip = 'Standardized Upper confidence interval!'; % Actually, put something more
                                        % informative here

Vout=spm_write_vol(Vout,UpperCon);

end

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
