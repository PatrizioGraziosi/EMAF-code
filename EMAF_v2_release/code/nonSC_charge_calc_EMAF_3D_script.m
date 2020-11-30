% Copyright 2020 Patrizio Graziosi, Neophytos Neophytou                   %
% Written and developed by                                                %
% Patrizio Graziosi and Neophytos Neophytou                               %
% patrizio.graziosi@cnr.it and n.neophytou@warwick.ac.uk                  %
% under the                                                               %
% Marie Curie - Individual Fellowships -  GENESIS - project ID 788465     %
% Generic transport simulator for new generation thermoelectric materials %
% ----------------------------------------------------------------------- %
% This file is distributed under the terms of the GNU                     %
% General Public License. See the file `LICENSE' in  the root directory   %
% of the present distribution.                                            %
% ----------------------------------------------------------------------- %
%                                                                         %
% Please cite the code source and the paper [1] when publishing results   %
% obtained  using the present EMAF code                                   %
% [1] P. Graziosi, C. Kumarasinghe, N. Neophytou, ACS Appl. Energy Mater. 
% 3, 6, 5913 (2020) ; https://pubs.acs.org/doi/abs/10.1021/acsaem.0c00825 %
%                                                                         %
% ----------------------------------------------------------------------- %


%--------- Input parameters ------------  

%---- % Physical Constants ---------------------------------------
%-----------------------------------------------------------------
m0=9.1094e-31;            % [kg]
hbar=(6.6261e-34)/(2*pi); % [J-sec]
q0=1.609e-19;             % [col]
% T=300;                    % [K]
kB=1.38e-23;              % [J/K]
kT_q = kB*T/q0;           % [eV]   
%------------------------------------------------------------------------

EF_2_shift = 40*kT_q; % 1.5; % represents V_SD, that shall be large, eq. A8 in IEEE Trans. Electr. Dev. 50, 1853 (2003), "Theory of Ballistic Nanotransistors" 
dEF = 0.05;
%----------------------------------------------


%--------- Initializations ---------------------------------
% kT_q = kB*T/q0; % in eV
Ek_reduced_mult_factor=1;
%------------------------------------

E_edge = min(min(min(min(Ek)))); % the code shifting_bands_VOMBATO that runs in the EMAF_main shifts the bands to zero, this is for generality

%------------- Get the VELOCITIES of the k-states ------
clear Vx_transp;
clear Vy_transp;
clear Vz_transp;
nkx=size(Ek,1); nky=size(Ek,2); nkz=size(Ek,3);
Vx_transp=zeros(size(Ek,1),size(Ek,2),size(Ek,3),size(Ek,4)); Vy_transp=zeros(size(Ek,1),size(Ek,2),size(Ek,3),size(Ek,4)); Vz_transp=zeros(size(Ek,1),size(Ek,2),size(Ek,3),size(Ek,4));
for ii_bands=1:size(pseudoconduction_bands,2) %n_bands_transp
    % central differences with dummies variables ix, iy iz        
    for id_x = 1 : nkx
       for id_y = 1 : nky
           for id_z = 1 : nkz
               if id_x == 1 
                   xd = 2; xu = id_x;
               elseif id_x == nkx
                   xu = nkx-1; xd = id_x;
               else
                   xu = id_x; xd = id_x;
               end
               if id_y == 1
                   yu = id_y; yd = 2;
               elseif id_y == nky
                   yu = nky-1; yd = id_y;
               else
                   yu = id_y; yd = id_y;
               end
               if id_z == 1
                   zd = 2; zu = id_z;
               elseif id_z == nkz
                   zu = nkz-1; zd = id_z;
               else
                   zu = id_z; zd = id_z;
               end
                                  
                dkx_v = [ kx_matrix(xu+1,id_y,id_z) ky_matrix(xu+1,id_y,id_z) kz_matrix(xu+1,id_y,id_z) ] - [ kx_matrix(xd-1,id_y,id_z) ky_matrix(xd-1,id_y,id_z) kz_matrix(xd-1,id_y,id_z) ];
                dky_v = [ kx_matrix(id_x,yu+1,id_z) ky_matrix(id_x,yu+1,id_z) kz_matrix(id_x,yu+1,id_z) ] - [ kx_matrix(id_x,yd-1,id_z) ky_matrix(id_x,yd-1,id_z) kz_matrix(id_x,yd-1,id_z) ];
                dkz_v = [ kx_matrix(id_x,id_y,zu+1) ky_matrix(id_x,id_y,zu+1) kz_matrix(id_x,id_y,zu+1) ] - [ kx_matrix(id_x,id_y,zd-1) ky_matrix(id_x,id_y,zd-1) kz_matrix(id_x,id_y,zd-1) ];
                dkx_temp = norm(dkx_v); dky_temp = norm(dky_v); dkz_temp = norm(dkz_v);
                
                Vx_transp(id_x,id_y,id_z,ii_bands) = q0/hbar * ( Ek(xu+1,id_y,id_z,ii_bands ) - Ek(xd-1,id_y,id_z,ii_bands) ) / dkx_temp;
                Vy_transp(id_x,id_y,id_z,ii_bands) = q0/hbar * ( Ek(id_x,yu+1,id_z,ii_bands ) - Ek(id_x,yd-1,id_z,ii_bands) ) / dky_temp;
                Vz_transp(id_x,id_y,id_z,ii_bands) = q0/hbar * ( Ek(id_x,id_y,zu+1,ii_bands ) - Ek(id_x,id_y,zd-1,ii_bands) ) / dkz_temp;          
                                        
           end
       end
    end
end





%-----------Fermi---------------------------------
EF_scan = E_edge-0.5:dEF:0.7; % represents V_G

Charge_S_EF = zeros(1,size(EF_scan,2));
Charge_D_EF = zeros(1,size(EF_scan,2));
%------------------------------------------




%--------- Ef scan of Ek --------------------------------------------

for i_EF = size(EF_scan,2):-1:1
    EF_1_tmp = EF_scan(i_EF);
    EF_2_tmp = EF_1_tmp - EF_2_shift;

    
    Charge_source = Ek_reduced_mult_factor*0.5*(dVk/(8*pi^3)).*(1./(1+exp((Ek-EF_1_tmp)/kT_q))).*(strcmp(Spin_Orbit_flag,'off')+1);                   
    Charge_drain  = Ek_reduced_mult_factor*0.5*(dVk/(8*pi^3)).*(1./(1+exp((Ek-EF_2_tmp)/kT_q))).*(strcmp(Spin_Orbit_flag,'off')+1);
                 
 
    IStmp = q0*Charge_source .* abs(Vx_transp); IS_x = sum(sum(sum(IStmp)));
    IDtmp = q0*Charge_drain .* abs(Vx_transp);  ID_x = sum(sum(sum(IDtmp)));
    I_TB_tmp_x = sum( IS_x - ID_x ); 
    IStmp = q0 * Charge_source .* abs(Vy_transp); IS_y = sum(sum(sum(IStmp)));
    IDtmp = q0 * Charge_drain .* abs(Vy_transp);  ID_y = sum(sum(sum(IDtmp)));
    I_TB_tmp_y = sum( IS_y - ID_y );
    IStmp = q0 * Charge_source .* abs(Vz_transp); IS_z = sum(sum(sum(IStmp)));
    IDtmp = q0 * Charge_drain .* abs(Vz_transp);  ID_z = sum(sum(sum(IDtmp)));
    I_TB_tmp_z = sum( IS_z - ID_z );
    
    I_array(i_EF) = ( I_TB_tmp_x + I_TB_tmp_y + I_TB_tmp_z ) / 3;


    % The injection velocities
    v_inj_non_SC_x(i_EF) = I_TB_tmp_x./sum(sum(sum(sum((squeeze(Charge_source))))))/q0;  
    v_inj_non_SC_S_D_x(i_EF) = I_TB_tmp_x./sum(sum(sum(sum((squeeze(Charge_source-Charge_drain))))))/q0;
    v_inj_non_SC_y(i_EF) = I_TB_tmp_y./sum(sum(sum(sum((squeeze(Charge_source))))))/q0;  
    v_inj_non_SC_S_D_y(i_EF) = I_TB_tmp_y./sum(sum(sum(sum((squeeze(Charge_source-Charge_drain))))))/q0;  
    v_inj_non_SC_z(i_EF) = I_TB_tmp_z./sum(sum(sum(sum((squeeze(Charge_source))))))/q0;  
    v_inj_non_SC_S_D_z(i_EF) = I_TB_tmp_z./sum(sum(sum(sum((squeeze(Charge_source-Charge_drain))))))/q0;
%     v_inj_non_SC(i_EF) = I_TB_tmp ./ sum(sum(sum(sum((squeeze(Charge_source))))))/q0;  
%     v_inj_non_SC_S_D(i_EF) = I_TB_tmp ./ sum(sum(sum(sum((squeeze(Charge_source-Charge_drain))))))/q0;     

    Charge_S_EF(i_EF) = sum(sum(sum(sum(Charge_source))));  
    Charge_D_EF(i_EF) = sum(sum(sum(sum(Charge_drain))));
    
end  %
%--------- END Ef scan of Ek --------------------------------------------           

Charge_total=Charge_S_EF+Charge_D_EF;


% --- Extraction of the thermal velocities and the equivalent masses -----
clear v_inj_array

% what follows is for the conductivity effective mass, that does 
% reproduce the TDF
v_inj_min_x = min(v_inj_non_SC_x);
v_inj_array_x = [v_inj_min_x];
for i_v = 1:size(v_inj_non_SC_x,2)
    if (v_inj_non_SC_x(i_v)-v_inj_min_x)/v_inj_min_x < 0.003
        v_inj_array = [v_inj_array_x, v_inj_non_SC_x(i_v)];
    end
end
v_inj_thermal_x = mean(v_inj_array_x);
m_inj_x = 1/( v_inj_thermal_x^2 * pi / (2* kB * T) ) / m0 ;

v_inj_min_y = min(v_inj_non_SC_y);
v_inj_array_y = [v_inj_min_y];
for i_v = 1:size(v_inj_non_SC_y,2)
    if (v_inj_non_SC_y(i_v)-v_inj_min_y)/v_inj_min_y < 0.003
        v_inj_array = [v_inj_array_y, v_inj_non_SC_y(i_v)];
    end
end
v_inj_thermal_y = mean(v_inj_array_y);
m_inj_y = 1/( v_inj_thermal_y^2 * pi / (2* kB * T) ) / m0 ;

v_inj_min_z = min(v_inj_non_SC_z);
v_inj_array_z = [v_inj_min_z];
for i_v = 1:size(v_inj_non_SC_z,2)
    if (v_inj_non_SC_z(i_v)-v_inj_min_z)/v_inj_min_z < 0.003
        v_inj_array = [v_inj_array_z, v_inj_non_SC_z(i_v)];
    end
end
v_inj_thermal_z = mean(v_inj_array_z);
m_inj_z = 1/( v_inj_thermal_z^2 * pi / (2* kB * T) ) / m0 ;

m_equivalent_injection = 3 / (1/m_inj_x + 1/m_inj_y + 1/m_inj_z);

% -----------------------------------------------------------------------





if strcmp(plot_velocities,'yes')
        %------------------------------------------------------
        %---------- Plotting the velocities vs nL ----------
        %------------------------------------------------------
    if strcmp(carriers, 'electrons')
         %--- The charge vs. Ef --------------
        figure(1);
        h1 = plot(EF_scan-E_edge,Charge_total,'-b','linewidth',[lwpl]);hold on;
        set(gca,'Fontsize',[fsize],'linewidth',[lwbor]);
        xlabel('E_f [eV]');
        ylabel('n [m^{-3}]');

        %--- Plotting the injection velocities of TB
        figure(2);
        h1 = semilogx(Charge_S_EF+Charge_D_EF,v_inj_non_SC_S_D_x,'-b'); hold on;
        set(gca,'Fontsize',[fsize],'linewidth',[lwbor]);
        set(h1,'linewidth',[lwpl]);
        xlabel('n [m^{-3}]');
        ylabel('v_{inj} [m/s]');

        figure(3);
        h1 = plot(EF_scan-E_edge,v_inj_non_SC_S_D_x,'-b'); hold on;
        set(gca,'Fontsize',[fsize],'linewidth',[lwbor]);
        set(h1,'linewidth',[lwpl]);
        xlabel('E_F [eV]');
        ylabel('v_{inj} [m/s]');


        %--- The TB current --------------
        figure(4);
        h1 = semilogx(Charge_total,I_array,'-b','linewidth',[lwpl]);hold on;
        set(gca,'Fontsize',[fsize],'linewidth',[lwbor]);
        xlabel('n [m^{-^3}]');
        ylabel('I [A]')

        figure(5);
        h1 = plot(EF_scan-E_edge,I_array,'b-','linewidth',[lwpl]);hold on;
        set(gca,'Fontsize',[fsize],'linewidth',[lwbor]);
        xlabel('E_f [eV]');
        ylabel('I [A]')

    elseif strcmp(carriers, 'holes')
         %--- The charge vs. Ef --------------
        figure(1);
        h1 = plot(EF_scan-E_edge,Charge_total,'-r','linewidth',[lwpl]);hold on;
        set(gca,'Fontsize',[fsize],'linewidth',[lwbor]);
        xlabel('E_f [eV]');
        ylabel('n [m^{-3}]');

        %--- Plotting the injection velocities of TB
        figure(2);
        h1 = semilogx(Charge_S_EF+Charge_D_EF,v_inj_non_SC_S_D_x,'-r'); hold on;
        set(gca,'Fontsize',[fsize],'linewidth',[lwbor]);
        set(h1,'linewidth',[lwpl]);
        xlabel('n [m^{-3}]');
        ylabel('v_{inj} [m/s]');

        figure(3);
        h1 = plot(EF_scan-E_edge,v_inj_non_SC_S_D_x,'-r'); hold on;
        set(gca,'Fontsize',[fsize],'linewidth',[lwbor]);
        set(h1,'linewidth',[lwpl]);
        xlabel('E_F [eV]');
        ylabel('v_{inj} [m/s]');


        %--- The TB current --------------
        figure(4);
        h1 = semilogx(Charge_total,I_array,'-r','linewidth',[lwpl]);hold on;
        set(gca,'Fontsize',[fsize],'linewidth',[lwbor]);
        xlabel('n [m^{-^3}]');
        ylabel('I [A]')

        figure(5);
        h1 = plot(EF_scan-E_edge,I_array,'-r','linewidth',[lwpl]);hold on;
        set(gca,'Fontsize',[fsize],'linewidth',[lwbor]);
        xlabel('E_f [eV]');
        ylabel('I [A]')
    end
    
end
