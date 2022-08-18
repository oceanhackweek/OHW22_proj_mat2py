function [tr,zc,Ug,Vg,F] = read_hourly_ADCP2(ffnm)

  if nargin == 0 
    %ffnm = '~/Downloads/IMOS_ANMN-NRS_VZ_20080623_NRSYON_FV02_velocity-hourly-timeseries_END-20210623_C-20211222.nc'; 
  end

 F = nc2struct(ffnm);
 A = F.global_attributes;
 isite = strmatch('site_code',A(:,1));
 lcn = A{isite,2};

  time = F.TIME;
  zz = F.DEPTH;
  ucur = F.UCUR;
  vcur = F.VCUR;
  depth_count = F.DEPTH_count;
  inidx = F.instrument_index;
  znom = F.NOMINAL_DEPTH;
  nbins = max(unique(F.CELL_INDEX));

  dz = round(round(median(znom))/nbins*2)/2;
  %this dz is chosen to get a single velocity per depth level

  tref = min(time);
  tr = tref:1/24:max(time);
  ti = int64(round((time-tref)*24)+1);

  zr = [0:dz:round(max(zz))+dz]';
  Zc = av2(zr);
  zi = int64(round((zz-Zc(1))/dz)+1);
  Ug = nan([length(Zc) length(tr)]);
  Vg = Ug;
  zg = Ug;
  Ndat = zg;

  for ii = 1:length(inidx)
   if ~isnan(ucur(ii)) & ~isnan(zz(ii)) & zz(ii)>0 & depth_count(ii)>1
        Ug(zi(ii),ti(ii)) = ucur(ii);
        Vg(zi(ii),ti(ii)) = vcur(ii);
        zg(zi(ii),ti(ii)) = zz(ii);
        Ndat(zi(ii),ti(ii)) = depth_count(ii);
    end
  end

 % NRSYON has problems in0904 & 1206 deployments
 % find large values and set to nan
 if strcmp(lcn,'NRSYON') 
   ibg = find(abs(Vg)>1.2 | abs(Ug)>1.2 );
   Ug(ibg) = nan;
   Vg(ibg) = nan;
 end

%find top viable level
inn = find(~isnan(nanmean(Ug,1)));
nnn = length(inn);
for iz = 1:length(Zc)
  pcgd(iz) = sum(~isnan(Ug(iz,:)))/nnn;
end

% NB this extraction leaves lots of nans due to dz being same as original
% vertical bin size and motion of the tide
% average over nb to get a more robust time series (2.5m)
igd = find(pcgd>0.4);
itop = igd(1); nz = length(Zc-itop); nb = 5; nbin = ceil(nz/nb); nt = length(tr);
ug = nan(nbin,nt); vg = ug; ndat = ug; zc = nan(nbin,1);

for ii = 1:nbin
  i1 = (ii-1)*nb+itop; iav = i1:min(i1+nb-1,nz);
  if length(iav)>=3
    ug(ii,:) = nanmean(Ug(iav,:),1);
    vg(ii,:) = nanmean(Vg(iav,:),1);
    ndat(ii,:) = nansum(Ndat(iav,:),1);
    zc(ii) = nanmean(Zc(iav))';
  end
end

tdiff = F.SECONDS_TO_MIDDLE;
% what does it mean for deployment to have nan for seconds to middle
% find(isnan(tdiff))

% remove bins that are all nans
% noting that there may be a number of deployment gaps
nt = length(tr);
inn = find(~isnan(nanmean(Ug,1)));
unn = ug(:,inn);
pcnan = sum(isnan(unn),2)/size(unn,2);
%igz = find(pcnan<0.5);
%nz = length(igz);
% initial choice of pcnan<0.5 is based on the combination of:
% the increasing pcnans at surface and bottom and variability of the mean
% next step - check the spread of nans & depth count through the ts

plot_check = 0;
if plot_check
% check for number of obs going into each hourly average along each z level
% sometimes when ADCP is sampling at 20min intervals there are only 2/hour 
%e.g. Yongala -pre 2011
% identify which z-levels are consistently getting <=1 obs/hour

vvs = {'Ndat';'Ug';'Vg'};
ttls = {'Number of records going to hourly value in LTSP'; ...
        'Eastward velocity'; ...
        'Northward velocity'};

for iv = 1:length(vvs)
    vnam = vvs{iv};
    eval(['VV = ' vnam ';'])
    figure('position',[0 0 1200 400],'color',[1 1 1])
    set(gca,'position',[0.05 0.09 0.9 0.84])
    xcolor(tr,-Zc,VV)
    colorbar
    hold on
    %plot(tr,repmat(-zc,[1 length(tr)]),'.k','markersize',1)
    datetick('x','yyyy','keeplimits')
    set(gca,'tickdir','in')
    title([lcn ' ' ttls{iv} ' - output velocity levels indicated'])
    ylabel('Depth')
    ax1 = axis;
end

vvs = {'ug';'vg'}; cax = [-1 1];
for iv = 1:length(vvs)
    vnam = vvs{iv};
    eval(['VV = ' vnam ';'])
    figure('position',[0 0 1200 400],'color',[1 1 1])
    set(gca,'position',[0.05 0.09 0.9 0.84])
    xcolor(tr,-zc,VV)
    colorbar
    hold on
    plot(tr,repmat(-zc,[1 length(tr)]),'.k','markersize',1)
    datetick('x','yyyy','keeplimits')
    caxis(cax)
    set(gca,'tickdir','in')
    title([lcn ' binned ' vnam])
    ylabel('Depth')
    axis(ax1)
end
end
nnn = sum(~isnan(nanmean(ug,1)));
pcgd = sum(~isnan(ug),2)/nnn*100;
igz = find(pcgd>.4);
nz = length(igz);

uu = ug(igz,:);
vv = vg(igz,:);
ub = nanmean(uu,1);
vb = nanmean(vv,1);
zc = zc(igz);

% before filtering fill gaps that are less than 6 hours big
uup = nan(size(uu));
vvp = nan(size(vv));

for iz = 1:length(zc)
  inn = find(~isnan(uu(iz,:)));
  gaps = find(diff(tr(inn))>1.5/24);
  ngaps = length(gaps);
  %disp(['iz: ' int2str(iz) '  ngaps: ' int2str(ngaps)])
  uup(iz,:) = gapfil(tr,uu(iz,:),6,24);
  vvp(iz,:) = gapfil(tr,vv(iz,:),6,24);
end

nnan1 = sum(isnan(uu),2);
nnan2 = sum(isnan(uup),2);
if plot_check
figure
plot(nnan1,-zc,'.','markersize',10)
hold on
plot(nnan2,-zc,'.','markersize',10)
end
pcnan = nnan2/length(tr);

% remove levels where after filling gaps there are too many nans
iok = find(pcnan<.40);  % more than 40% nans
uu = uup(iok,:);
vv = vvp(iok,:);
ub = nanmean(uup,1);
vb = nanmean(vvp,1);
zc = zc(iok);
nz = length(zc);

% time mean velocity at each depth
umn = nanmean(uup,2);
vmn = nanmean(vvp,2);
% direction of barotropic current
[mm,b0] = linear_fit(ub,vb);
thd = 90-atand(mm);
thz = nan(nz,1);
thzf = nan(nz,1);
uf = nan(nz,nt);
vf = nan(nz,nt);
ulp = uf; vlp = uf; uhp = uf; vhp = uf; 
nav = 40*2/median(diff(tr))/24+1;
for ii = 1:nz
  [mm,b0] = linear_fit(uu(ii,:),vv(ii,:));
  % get dirn of major axis of flow (tides included)
  % atand gives angle measured anticlockwise from x-axis
  % 90-atand gives angle measured clockwise from y-axis
  thz(ii) = 90-atand(mm);  
  uf(ii,:) = hanfilmc(uu(ii,:),nav,0);
  vf(ii,:) = hanfilmc(vv(ii,:),nav,0);
  % compare tsf3 wavelet filter 
  %[ulp(ii,:),uhp(ii,:),wave,period,ttf] = tsf3(uup(ii,:),tr,40,0,24);
  %[vlp(ii,:),vhp(ii,:),wave,period,ttf] = tsf3(vvp(ii,:),tr,40,0,24);
  % now get the direction of the low-pass velocity
  [mm,b0] = linear_fit(uf(ii,:),vf(ii,:));
  thzf(ii) = 90-atand(mm);  
end

% subsample filtered velocities
iss = 12:24:length(tr);
ufs = uf(:,iss);
vfs = vf(:,iss);
ttf = tr(iss);

% time mean of filtered velocity at each depth
 umnf = nanmean(ufs,2);
 vmnf = nanmean(vfs,2);
 % flip u&v in atan2d to get clockwise from 90
 % direction of the mean
 thz_mn = atan2d(umnf,vmnf);

% estimate local alongshore direction using filtered velocities
% in the middle of the watercolumn
  if median(znom) > 90
    imid = find(zc>20 & max(abs(znom))-abs(zc)>20);
  else 
    imid = 1:length(thzf);
  end
  th0 = nanmean(thz_mn(imid));
  th1 = nanmean(thzf(imid));
  
% rotate velocites to get across-shelf and alongshore
% vrotate rotates clockwise 
[ur,vr] = vrotate(ufs,vfs,-th0,'degrees');
umr = nanmean(ur,2); 
vmr = nanmean(vr,2);

vvs = {'ur';'vr'}; cax = [-.5 .5];
vlbl = 'rotated'
for iv = 1:length(vvs)
    vnam = vvs{iv};
    eval(['VV = ' vnam ';'])
    figure('position',[0 0 1200 400],'color',[1 1 1])
    set(gca,'position',[0.05 0.09 0.9 0.84])
    %xcolor(ttf,-zc,VV)
    %colorbar
    plot(ttf,VV,'.')
    hold on
    %plot(ttf,repmat(-zc,[1 length(ttf)]),'.k','markersize',1)
    datetick('x','yyyy','keeplimits')
    caxis(cax)
    set(gca,'tickdir','in')
    title([lcn ' ' vlbl ' ' vnam])
    ylabel('Depth')
    %axis(ax1)
end

  
keyboard
figure
ZI = zeros(size(zc));
quiver3(ZI,ZI,-zc,umr,vmr,ZI)

 keyboard
