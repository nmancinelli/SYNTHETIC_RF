
function [R,T,Z,tt] = write_propmat_syn(zt,zb,vp,vs,rh,RayParam,PERIOD)

ANGL=asind(RayParam*vs(end));

fprintf('Generating synthetics for p, angle, period = %f, %f, %f\n',RayParam,ANGL,PERIOD);

nlay=length(zt);
LAYmodel.nlay = nlay;

LAYmodel.zlayt =zt;
LAYmodel.zlayb =zb;
LAYmodel.Vp    =vp;
LAYmodel.Vs    =vs;
LAYmodel.rho   =rh;
LAYmodel.VpAnis=zeros(nlay,1)';
LAYmodel.VsAnis=zeros(nlay,1)';
LAYmodel.faz   =zeros(nlay,1)';
LAYmodel.fpl   =zeros(nlay,1)';
LAYmodel.intstr=zeros(nlay,1)';
LAYmodel.intdip=zeros(nlay,1)';

%%Print velocity model
% FID=fopen('vmodel.txt','w');
% 
% fprintf('%10s %10s %10s %10s %10s\n','ZTOP','ZBOT','VP','VS','RHO');
% 
% for ii = 1:length(zt);
% 
% fprintf('%10f %10f %10f %10f %10f\n',zt(ii),zb(ii),vp(ii),vs(ii),rh(ii));
% fprintf(FID,'%10f %10f %10f %10f\n',zt(ii),vp(ii),vs(ii),rh(ii));
% fprintf(FID,'%10f %10f %10f %10f\n',zb(ii),vp(ii),vs(ii),rh(ii));
% 
% end
% 
% fclose(FID);
%%
%Step 2 - Generate synthetics
%cd('/Users/mancinelli/PROJECTS/ARRAY_STACK/ReceiverFunctions/Make_Receiver_Functions/CLEAN/Scattered_Waves/SYN/PROPMAT/SUBS/matlab_to_propmat')
[traces,tt,~,~] = PROPMAT.run_propmat(LAYmodel,'eg','Sp',10,ANGL,PERIOD);

R=traces(:,1);
T=traces(:,2);
Z=traces(:,3);

end
