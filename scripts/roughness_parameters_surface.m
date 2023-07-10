function param = roughness_parameters_surface(surf, z)

% Gadelmawla et al., 2002
% https://doi.org/10.1016/S0924-0136(02)00060-2

% preparation
% finds peaks and valleys based on values only
pks_locs = z > prctile(z,95);
vls_locs = z < prctile(z,5);
% intersection of mean line
edg = SurfStatEdg(surf);
n0_list = [0 0];
for ii = 1:length(edg)
    if sum(sign(z(edg(ii,:)) - mean(z)))==0
        edg_n0(ii) = 1;
    else
        edg_n0(ii) = 0;
    end
end
n0_list = edg(edg_n0==1,:);

% AMPLITUDE PARAMETERS
% arithmetic average height, also known as centre line average
Ra = mean(abs(z - mean(z)));
% RMS roughness - standard deviation of surface heights
Rq = std(z);
% Maximum height of peaks 
Rp = diff([mean(z); max(z)]);
% Maximum depth of valleys
Rv = diff([min(z); mean(z)]);
% Mean height of peaks
if sum(pks_locs)~=0
    Rpm = mean(z(pks_locs)) - mean(z);
else
    Rpm = 0;
end
% Mean depth of valleys
if sum(vls_locs)~=0
    Rvm = mean(z) - mean(z(vls_locs));
else
    Rvm = 0;
end
% Maximum height of the profile
Rmax = max(z);
% Profile solidarity factor
k = Rv/Rmax;
% Skewness
Rsk = skewness(z);
% Kurtosis
Rku = kurtosis(z);

% SPACING PARAMETERS
% Peak count - adapted to be unique peaks in the landscape
rand_data = rand(10, length(surf.coord));
slm = SurfStatLinMod(rand_data,1,surf);
slm = SurfStatT(slm, ones(10,1));
slm.t = pks_locs';
[~, ~, clusid] = SurfStatPeakClus(slm, ones(1,length(pks_locs)), 0.5);
if ~isempty(clusid)
    Pc = max(clusid);
else
    Pc = 0;
end
% Mean spacing of adjacent local peaks - adapted to avg. Euclidean distance
% between centroids
pks_cent = [];
pks_cent = grpstats(surf.coord',clusid,'mean');
S = mean(pdist(pks_cent(2:end,:)));
% Mean spacing at mean line - adapted to avg. Euclidean distance
tmp = [];
for ii = 1:length(n0_list)
    tmp(ii,:) = pdist(surf.coord(:,n0_list(ii,:))');
end
Sm = mean(tmp);
% Number of intersections of the mean line, relative to size of surface -
% adapted from length to sqrt(area)
n0 = length(n0_list);
% mean radius or asperities - adapted be average size of peaks
tmp = tabulate(clusid);
if ~size(tmp,1) > 1
    rp = mean(tmp(2:end,2));
else
    rp = 0;
end

% HYBRID PARAMETERS
% Slope at mean line  - adapted to absolute slope
tmp = [];
for ii = 1:length(n0_list)
    tmp(ii) = abs(diff(z(n0_list(ii,:))));
end
gamma = mean(tmp);
% Mean slope of surface - adapted to absolute slope
delta_a = mean(abs(diff(z)));
% RMS of slope of surface
delta_q = sqrt(mean(abs(diff(z)).^2));
% Average wavelength
lambda_a = 2*pi*Ra*delta_a;
% RMS wavelength
lambda_q = 2*pi*Rq*delta_q;
% Steepness factor of the surface
Sf = Ra*Sm;
% Waviness factor of the surface
Wf = numel(z)/Ra;

param = table(Ra, Rq, Rp, Rv, Rpm, Rvm, Rmax, k, Rsk, Rku, Pc, S, Sm, n0, rp, ...
    gamma, delta_a, delta_q, lambda_a, lambda_q, Sf, Wf, ...
    'VariableNames', {'Ra', 'Rq', 'Rp', 'Rv', 'Rpm', 'Rvm', 'Rmax', 'k', ...
    'Rsk', 'Rku', 'Pc', 'S', 'Sm', 'n0', 'rp', 'gamma', 'delta_a', ...
    'delta_q', 'lambda_a', 'lambda_q', 'Sf', 'Wf'});