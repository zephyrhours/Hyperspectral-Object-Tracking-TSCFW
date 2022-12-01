function H=hog3d(data,cellSize,hOrients,vOrients)
sz=size(data);
data=single(data);
% data=gpuArray(data);
[Gx,Gy, Gz]=gradient(data,1);
theta = atan2(Gy,Gx)+pi;
phi = atan2(Gz, sqrt((Gx).^2 + (Gy).^2))+pi/2;
magnitudes = sqrt((Gx).^2 + (Gy).^2 + (Gz).^2);
% [magnitudes, theta, phi] = imgradient3(data,'prewitt');
% theta=theta./180*pi+pi;
% phi=phi./180*pi+pi/2;
% histograms=zeros([floor(sz/cellSize),hOrients+vOrients]);
hz=[floor(sz(1:2)/cellSize),(hOrients+vOrients)*floor(sz(3)/cellSize)*4];
% H1=zeros([hz(1:2),hOrients]);
% H2=zeros([hz(1:2),vOrients]);
cols=1:cellSize:sz(3);
%gradient Quantize
thetaCell = num2cell(theta,[1,2]);
phiCell=num2cell(phi,[1,2]);
magnitudesCell=num2cell(magnitudes,[1,2]);
soft=-1;%interplot in horizental
clip=0.2;
thetaQ=cellfun(@(m,o) gradientHist(m,o,cellSize,hOrients,soft,0,clip,true),magnitudesCell, thetaCell,'UniformOutput',false);
thetaQ=cellfun(@(x) x(:), thetaQ,'UniformOutput',false );
thetaQ=squeeze(cell2mat(thetaQ));
thetaQ=arrayfun(@(x)sum(thetaQ(:,x:x+cellSize-1),2),cols,'UniformOutput',false);
thetaQ=cellfun(@(x) reshape(x,hz(1),hz(2),hOrients),thetaQ,'UniformOutput',false);
% thetaNorm=cellfun(@(x) reshape(x,hz(1),hz(2),hOrients) , (cellfun(@(x) cell2mat(arrayfun(@(y) conv2(pad(x(:,:,y).^2),[1 1;1 1],'valid')+0.01,1:hOrients,'UniformOutput',false)), thetaQ,'UniformOutput',false )),'UniformOutput',false);
% theta=cellfun(@(x,y) x./y,thetaQ,thetaNorm,'UniformOutput',false);
thetaBlock=(cellfun(@(x)cellfun(@(x) x',(arrayfun(@(y)im2colstep(double(pad(squeeze(x(:,:,y)))),[2,2]),1:hOrients,'UniformOutput',false)),'UniformOutput',false), thetaQ,'UniformOutput',false ));
thetaBlock=(cellfun(@(x) cell2mat(x),thetaBlock,'UniformOutput',false));

thetaHist=cell2mat(cellfun(@(x) normlize(x),thetaBlock,'UniformOutput',false));


phiQ=cellfun(@(m,o) gradientHist(m,o,cellSize,vOrients,soft,0,clip,false),magnitudesCell, phiCell,'UniformOutput',false);
phiQ=cellfun(@(x) x(:), phiQ,'UniformOutput',false );
phiQ=squeeze(cell2mat(phiQ));
phiQ=arrayfun(@(x)sum(phiQ(:,x:x+cellSize-1),2),cols,'UniformOutput',false);
phiQ=cellfun(@(x) reshape(x,hz(1),hz(2),vOrients),phiQ,'UniformOutput',false);
phiBlock=(cellfun(@(x)cellfun(@(x) x',(arrayfun(@(y)im2colstep(double(pad(squeeze(x(:,:,y)))),[2,2]),1:vOrients,'UniformOutput',false)),'UniformOutput',false), phiQ,'UniformOutput',false ));
phiBlock=(cellfun(@(x) cell2mat(x),phiBlock,'UniformOutput',false));

phiHist=cell2mat(cellfun(@(x) normlize(x),phiBlock,'UniformOutput',false));
H=[thetaHist phiHist];
H=reshape(H,hz);
% [U,mu,vars] = pca(H' );
% [Y,Xhat,avsq] = pcaApply(H', U, mu,11);
% H=reshape(Y',hz(1),hz(2),11);








function x=pad(x)
        x=[x x(:,end-1)];
        x=[x;x(end-1,:)];
end
function z=normlize(x)
    sz=size(x);
   y=hogNormMatrix(x);
   z=x./(repmat(y,[1,sz(2)])+eps);
   z(z>0.2)=0.2;
   y=hogNormMatrix(z);
   z=z./(repmat(y,[1,sz(2)])+eps);
end


function x = normalizeL2Hys(x)
classToUse = class(x);
x = x./(norm(x,2) + eps(classToUse)); % L2 norm
x(x > 0.2) = 0.2;                     % Clip to 0.2
x = x./(norm(x,2) + eps(classToUse)); % repeat L2 norm
end
function  y=hogNormMatrix(x)
        y=sum(abs(x).^2,2).^(1/2);   
end
end
