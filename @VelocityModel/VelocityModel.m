classdef VelocityModel
    %VelocityModel Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        z
        vs
        vp
        rho
        zlayt
        zlayb
        vplay
        vslay
        rhlay
    end
    
    methods
        function VM = VelocityModel(path2cardfile)
            [VM.z, VM.vp, VM.vs, VM.rho] = VelocityModel.load_from_card_file(path2cardfile);
            [VM.zlayt,VM.zlayb,VM.vplay,VM.vslay,VM.rhlay] = VelocityModel.layerise(VM.z,VM.vp,VM.vs,VM.rho);
        end
        function plot(VM)
            hold on
            plot(VM.vs,VM.z);
            plot(VM.vp,VM.z);
            plot(VM.rho,VM.z);
            plot(VM.vslay,VM.zlayt,'o');
            plot(VM.vslay,VM.zlayb,'o');
            xlabel('Vp,Vs (km/s)');
            ylabel('Depth (km)')
            set(gca,'Ydir','reverse')
        end
        [Z,T] = migrator(VelocityModel,ray_parameter);
        
    end
    
    methods (Static)
       [depth, vpv, vsv, rho] = load_from_card_file ( path2cardfile )
       [ztop,zbot,vplay,vslay,rhlay] = layerise(Z,VP,VS,RH)
       
    end
    
end

