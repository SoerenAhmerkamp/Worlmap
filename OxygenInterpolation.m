

fname = 'L3m_20200201-20200229__GLOB_4_AV-MOD_CHL1_MO_00.nc';
ncdisp(fname)
%
vardata = ncread(fname,'CHL1_mean')';

imagesc(vardata)
caxis([0,1])

%geotiffwrite('test2.tif',[-180,-90;180,90],vardata)

%%


fname = 'A20150802015171.L3b_SNSP_CHL.x.nc';
ncdisp(fname)

A = h5read(fname,'/level-3_binned_data/chlor_a');

%

img = reshape(A.sum,5449,435);
caxis([0,1])

%geotiffwrite('test2.tif',[-180,-90;180,90],vardata)

%% Oxygen
data = m_shaperead('woa18_all_o00mn01');
load m_coasts
%%

%subplot(2,1,1)
lat = [];
lon = [];
oxy = [];
oxyInd = [];

for n = 1:length(data.ncst)
    temp = data.ncst{n};
    tempOx = [data.dbfdata{n,:}];
    tempOx(tempOx < 0) = NaN;
    
    [tempOxy, ind] = min(tempOx);
    oxyInd(end+1) = ind;
    
   % if val > -0.01
    lat(end+1) = temp(1);
    lon(end+1) = temp(2);
    oxy(end+1) = tempOxy;
    %end
end
hold on
scatter(lat,lon,100,oxy,'.')
plot(ncst(:,1),ncst(:,2),'-k')
%colormap(cmap)
%caxis([0.7, 6]), axis equal, axis tight
xlim([-180 180])
ylim([-90 90])
%colormap jet
colorbar

%% Standard Deviation 
data2 = m_shaperead('woa18_all_o00sd01');
load m_coasts
%%
subplot(2,1,2)
latStd = [];
lonStd = [];

for n = 1:length(data2.ncst)
    temp = data2.ncst{n};

    latStd(end+1) = temp(1);
    lonStd(end+1) = temp(2);
end
    
match = [];
matchStd = [];

% Check if stations are matched
for n = 1:length(lat)
    for m = 1:length(latStd)
        if lat(n) == latStd(m) && lon(n) == lonStd(m)
            match(end+1) = n; 
            matchStd(end+1) = m; 
        end
    end
end

%hold on
%plot(lat(match),lon(match),'.')
    
%plot(latStd(matchStd),lonStd(matchStd),'o')
    
%testLat = lat(match); testLon = lon(match);
%testLatStd = latStd(matchStd); testLonStd = lonStd(matchStd);

% Test lat lon
%hold on
%rval = round(rand(40,1)*100)+1;
%plot(testLat(rval),testLon(rval),'.')
%plot(testLatStd(rval),testLonStd(rval),'o')

%
%
%m_proj('WGS84');
%m_coast('color',[0 .6 0]);
%m_coast('color',[0 .6 0]);

latStd = latStd(matchStd);
lonStd = lonStd(matchStd);
stdev2 = [];
stdev = [];

for n = 1:length(matchStd)
    tempStd = data2.dbfdata{n,oxyInd(match(n))};

    if tempStd
        
    else
           tempStd = NaN;
    end

    stdev2(end+1) = tempStd;
    
    if oxy(match(n)) > 20 
        tempStd = NaN;
    end    
%
    ind = stdev > 0;
    stdev(end+1) = tempStd;
end

temp = oxy(match);

hold on
scatter(latStd,lonStd,100,stdev2,'.')
plot(ncst(:,1),ncst(:,2),'-k')
%colormap(cmap)
caxis([0, 20]), axis equal, axis tight
xlim([-180 180])
ylim([-90 90])
colormap parula

%% Standard Deviation
ind = ~isnan(stdev) & stdev >= 0;
x = latStd(ind);
y = lonStd(ind);
v = stdev(ind);

xGr = -180:1:180;
yGr = -90:1:90;

[xq,yq] = meshgrid(xGr,yGr);

mask = xq;

SN=4;
LX=0.8;
LY=0.8;

[zq eq]=divagrid(x,y,v,xq,yq); %,-SN,LX,LY
%vq = griddata(x,y,v,xq,yq);

hold on
pcolor(xq,yq,zq),shading flat
plot(ncst(:,1),ncst(:,2),'-k')
%colormap(cmap)
caxis([0, 20]), axis equal, axis tight
xlim([-180 180])
ylim([-90 90])
colormap parula

%% Oxygen Interpolation

ind = oxy >= -10; % Some values are below 0
oxy(oxy<0) = 0;
x = lat(ind);
y = lon(ind);
v = oxy(ind);

res = 1/2;


xGr = -180:res:180;
yGr = -90:res:90;

[xq,yq] = meshgrid(xGr,yGr);

hold on
siz = size(xq);
mask = zeros(siz);
for n = 1:length(k)-1
    pos = k(n)+1:k(n+1)-1;
    %fill(ncst(pos,1),ncst(pos,2),'r')
    xtem = ncst(pos,1);
    ytem = ncst(pos,2);
    %n% 9 13
    %A = 0.9;
    %xt = A*x+(1-A)*mean(xtem(1:end-1));
    %yt = A*y+(1-A)*mean(ytem(1:end-1));
    %pause
    temp = poly2mask((xtem+180)/res,(ytem+90)/res,siz(1),siz(2));
    mask = mask+temp;
end
mask(mask>1) = 0;
%mask(mask<0) = 0;
mask = imcomplement(mask);
mask(mask == 0) = NaN;
%pcolor(xq,yq,mask),shading flat

% siz = size(xq);
% mask = zeros(siz);
% 
% for n = 1:length(x)
%     
%     r = 4;
%     th = 0:pi/50:2*pi;
%     xunit = r * cos(th) + x(n);
%     yunit = r * sin(th) + y(n);
% 
%     temp = poly2mask(xunit+180,yunit+90,siz(1),siz(2));
%     mask = mask+temp;
% end
% maskTrans = mask;
% mask(mask>1) = 1;
%mask(mask == 0) = NaN;



xmax = max(xq(:));
ymax = max(yq(:));
xmin = min(xq(:));
ymin = min(yq(:));


SN=1;
LX=(xmax-xmin)/20;
LY=(ymax-ymin)/20;
   
   
[zq eq]=divagrid(x,y,v,xq,yq,SN,LX,LY,mask); %,-SN,LX,LY
%vq = griddata(x,y,v,xq,yq);

zq = fillmissing(zq,'nearest');

hold on
pcolor(xq,yq,zq),shading flat
fill(ncst(:,1),ncst(:,2),'r')
colormap(cmap), 
caxis([0, 4]), axis equal, axis tight
xlim([-180 180])
ylim([-90 90])

geotiffwrite('../logarithmicOxygen_HR.tif',[-179.5,-90;179.5,90],flipud(log10(zq)))
%geotiffwrite('../normalOxygen_LR.tif',[-179.5,-90;179.5,90],flipud(zq))



%% Standard Deviation

ind = ~isnan(stdev) & stdev >= 0;
ind2 = ~isnan(stdev2) & stdev >= 0;

x = latStd(ind);
y = lonStd(ind);
v = stdev(ind);
v2 = stdev2(ind);

%xMask = latStd(ind2);
%yMask = lonStd(ind2);

xGr = -180:1:180;
yGr = -90:1:90;

[xq,yq] = meshgrid(xGr,yGr);

% hold on
siz = size(xq);
mask = zeros(siz);

for n = 1:length(x)
    
    r = 5;
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x(n);
    yunit = r * sin(th) + y(n);

    temp = poly2mask(xunit+180,yunit+90,siz(1),siz(2));
    mask = mask+temp;
end
maskTrans = mask;
mask(mask>1) = 1;
mask(mask == 0) = NaN;


%pcolor(xq,yq,mask),shading flat

xmax = max(xq(:));
ymax = max(yq(:));
xmin = min(xq(:));
ymin = min(yq(:));


SN=1;
LX=(xmax-xmin)/20;
LY=(ymax-ymin)/20;
   
   
[zq eq]=divagrid(x,y,v2,xq,yq,SN,LX,LY,mask); %,-SN,LX,LY
%vq = griddata(x,y,v,xq,yq);

%zq = fillmissing(zq,'nearest');

hold on
pcolor(xq,yq,zq),shading flat
fill(ncst(:,1),ncst(:,2),'r')
colormap(cmap), 
caxis([0, 5]), axis equal, axis tight
xlim([-180 180])
ylim([-90 90])

%geotiffwrite('../transpMatrix.tif',[-179.5,-90;179.5,90],flipud(maskTrans))
%geotiffwrite('../normalOxygenSTD_LR.tif',[-179.5,-90;179.5,90],flipud(zq))



