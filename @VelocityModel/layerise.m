function [ztop,zbot,vplay,vslay,rhlay] = layerise(Z,VP,VS,RH)
% We must estimate a layercake structure for our velocity model
% Here's a really simple way to do it, which probably breaks down
% in the presence of strong gradients.  If you're going to use this
% routine, it's better to parameterize discontinuities as steps rather
% than gradients.


    function [avg] = get_layer_avg(V,ipt)
        avg = (V(ipt) + V(ipt+1))/2.;
    end

nlayMax = 100;

%Initialize vars
ilay = 0;
ztop = zeros(nlayMax,1);
zbot = zeros(nlayMax,1);
vplay = zeros(nlayMax,1);
vslay = zeros(nlayMax,1);
rhlay = zeros(nlayMax,1);

npts = length(Z);

for ipt = 1:(npts-1);
    %Check for true discontinuity
    if Z(ipt+1) == Z(ipt); continue; end
    ilay = ilay + 1;
    
    zbot(ilay) = Z(ipt);
    ztop(ilay) = Z(ipt+1);
    vplay(ilay) = get_layer_avg(VP,ipt);
    vslay(ilay) = get_layer_avg(VS,ipt);
    rhlay(ilay) = get_layer_avg(RH,ipt);
    
    %%Kluge to put in a thick MLD layer
    %if ztop(ilay) > 180 & zbot(ilay) < 250;
    %    vslay(ilay) = vslay(ilay)*0.90;
    %    vplay(ilay) = vplay(ilay)*0.90;
    %end
    
    if vslay(ilay) < 0.5;
        fprintf('***Warning, oceanic layer detected. Extending upper crust to surface.\n')
        vslay(ilay) = 3.2;
        vplay(ilay) = 5.8;
        rhlay(ilay) = 2.6;
    end
    
    
end

ztop =  flipud(ztop(1:ilay));
zbot =  flipud(zbot(1:ilay));
vplay = flipud(vplay(1:ilay));
vslay = flipud(vslay(1:ilay));
rhlay = flipud(rhlay(1:ilay));

end

    