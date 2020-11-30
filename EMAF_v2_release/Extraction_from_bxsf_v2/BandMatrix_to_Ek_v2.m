% Written and developed by                                                %
% Patrizio Graziosi, patrizio.graziosi@warwick.ac.uk, during the          %  
% Marie Curie - Individual Fellowships  GENESIS - project ID 788465       %
% Generic transport simulator for new generation thermoelectric materials %
% ----------------------------------------------------------------------- %
% This file is distributed under the terms of the GNU                     %
% General Public License. See the file `LICENSE' in  the root directory   %
% of the present distribution.                                            %
% ----------------------------------------------------------------------- %
%                                                                         %
% Please cite the code source when publishing results obtained            %
% using the present code                                                  %
%                                                                         %
% ----------------------------------------------------------------------- %

% BandMatrix has one column for each band
load('E_temp.mat')

Ek=zeros(points_in_axis_kx,points_in_axis_ky,points_in_axis_kz,num_of_bands);
for id_band = 1:num_of_bands
    id_k=1;
    for id_x=1:points_in_axis_kx
        for id_y=1:points_in_axis_ky
            for id_z=1:points_in_axis_kz
                Ek(id_x,id_y,id_z,id_band) = BandMatrix(id_k,id_band)-Fermi;                
                id_k=id_k+1;
            end
        end
    end
end


% alat shall be inputted      
if exist('blat','var') == 0 
    blat = alat;
end
if exist('clat','var') == 0 
    clat = alat;
end

% coordinates system change matrix, uc* vectors as rows
B_matrix = [ a*2*pi/(alat*1e-9) ; b*2*pi/(blat*1e-9); c*2*pi/(clat*1e-9) ] ;


    for id_x = (points_in_axis_kx - 1) : -1 : 0
        for id_y = (points_in_axis_ky - 1) : -1 : 0 
            for id_z = (points_in_axis_kz - 1) : -1 : 0
                
                k_vector_not_norm = [id_x id_y id_z]*B_matrix; % row vector * B matrix
                
                kx_matrix(id_x+1,id_y+1,id_z+1) = 1/points_in_axis_kx * k_vector_not_norm(1);
                ky_matrix(id_x+1,id_y+1,id_z+1) = 1/points_in_axis_ky * k_vector_not_norm(2);
                kz_matrix(id_x+1,id_y+1,id_z+1) = 1/points_in_axis_kz * k_vector_not_norm(3);
            end
        end
    end
    


clearvars -except Ek kx_matrix ky_matrix kz_matrix a b c kx_array ky_array kz_array bands_interpolation bands_centering nk_new material_name alat