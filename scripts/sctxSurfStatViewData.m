function [a,cb]=sctxSurfStatViewData(data, surf_L, surf_R, pos, clim);

% sctxSurfStatViewData is a simple viewer for subcortical surface data
% 
% Usage: [a, cb] = sctxSurfStatViewData(data, surf [,title [,background]]);
% 
% data        = 1 x v vector of data, v=# of subcortical structures (14)
%   one value per subcortical structure
%   order is L-accumbens, L-amygdala, L-caudate, L-hippocampus, L-pallidum
%               L-putamen, L-thalamus, R-accumbens, R-amygdala, R-caudate,
%                   R-hippocampus, R-pallidum, R-putamen, R-thalamus
%
% For both left and right subcortical surfaces:
%   surf.coord  = 3 x v matrix of coordinates.
%   surf.tri    = 3 x t matrix of triangle indices, 1-based, t=#triangles.
% title       = any string, data name by default.
% background  = background colour, any matlab ColorSpec, such as 
%   'white' (default), 'black'=='k', 'r'==[1 0 0], [1 0.4 0.6] (pink) etc.
%   Letter and line colours are inverted if background is dark (mean<0.5).
%
% a  = vector of handles to the axes, left to right, top to bottom. 
% cb = handle to the colorbar.
%
% Usage example:
%
% f = figure,
%   sctxSurfStatViewData(ones(14,1), sctx_template_L, sctx_template_R, 'ello!');
%
%
% Modified from BoSurfStatViewData.m
% SL | a rainy November night 2019


data = [ repmat(data(1), 867, 1); repmat(data(2), 1419, 1); ...
            repmat(data(3), 3012, 1); repmat(data(4), 3784, 1); ...
            repmat(data(5), 1446, 1); repmat(data(6), 4003, 1); ...
            repmat(data(7), 3726, 1); repmat(data(8), 838, 1); ...
            repmat(data(9), 1457, 1); repmat(data(10), 3208, 1); ...
            repmat(data(11), 3742, 1); repmat(data(12), 1373, 1); ...
            repmat(data(13), 3871, 1); repmat(data(14), 3699, 1) ];

vl   = 1:size(surf_L.coord, 2);
vr   = [1:size(surf_R.coord, 2)] + max(size(surf_L.coord, 2));

tl   = 1:size(surf_L.tri, 1);
tr   = [1:size(surf_R.tri, 1)] + max(size(surf_L.tri, 1));


% clf;
colormap(spectral(256));

a(1)=axes('position',pos(1,:));
trisurf(surf_L.tri,surf_L.coord(1,:),surf_L.coord(2,:),surf_L.coord(3,:),...
    double(data(vl)),'EdgeColor','none');
% view(0,90);
view(-90,0)
daspect([1 1 1]); axis tight; camlight; axis vis3d off;
lighting phong; material dull; shading flat;

a(2)=axes('position',pos(2,:));
trisurf(surf_L.tri,surf_L.coord(1,:),surf_L.coord(2,:),surf_L.coord(3,:),...
    double(data(vl)),'EdgeColor','none');
view(90,0); 
% view(90,-90)
daspect([1 1 1]); axis tight; camlight; axis vis3d off;
lighting phong; material dull; shading flat;

a(3)=axes('position',pos(3,:));
trisurf(surf_R.tri,surf_R.coord(1,:),surf_R.coord(2,:),surf_R.coord(3,:),...
    double(data(vr)),'EdgeColor','none');
view(-90,0); 
% view(-90,-90)
daspect([1 1 1]); axis tight; camlight; axis vis3d off;
lighting phong; material dull; shading flat;

a(4)=axes('position',pos(4,:));
trisurf(surf_R.tri,surf_R.coord(1,:),surf_R.coord(2,:),surf_R.coord(3,:),...
    double(data(vr)),'EdgeColor','none');
view(90,0);
% view(0,90);
daspect([1 1 1]); axis tight; camlight; axis vis3d off;
lighting phong; material dull; shading flat;


for i=1:length(a)
    set(a(i),'CLim',clim);
    set(a(i),'Tag',['SurfStatView ' num2str(i) ]);
end


cb=colorbar('location','South');
set(cb,'Position',[0.35 0.18 0.3 0.03]);
set(cb,'XAxisLocation','bottom');

dcm_obj=datacursormode(gcf);
set(dcm_obj,'UpdateFcn',@SurfStatDataCursor,'DisplayStyle','window');

%% colormaps 
% matlab original (boring)
colormap([.7 .7 .7; parula(256)])


return
end
