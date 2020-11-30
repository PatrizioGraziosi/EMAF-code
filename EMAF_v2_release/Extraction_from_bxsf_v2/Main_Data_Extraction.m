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

% input information
fileName = 'Mg3Sb2_fs.bxsf'; % name of the bxsf file
material_name = 'Mg3Sb2_interp'; % material name that will be used to save the band structure as Ek_'material_name'

num_of_bands = 5; % write the number of bands in the file


points_in_axis_kx = 61; % write the number of points in the k space
points_in_axis_ky = 61;
points_in_axis_kz = 41;

num_elements_each_row = 6; % how many energy values in each row, 6 is typical for bsxf

Fermi = 7.2709; % in eV, value of the Fermi level. 
               % In the Ek file, it will be placed to zero, NOTE!
               
alat = 0.4597 ; 


% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% num_lines_to_skip = 19; % number of line to be skipped in order to get
% the first ilne of energy values, 19 is the typical values in bsxf files
% cell_vector_width = 1; 
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------


% extraction codes run

num_of_points = points_in_axis_kx * points_in_axis_ky * points_in_axis_kz  ; 

Data_from_bxsf(fileName,num_of_bands,num_of_points,num_elements_each_row)

BandMatrix_to_Ek_v2


% additional features

if strcmp(bands_interpolation,'yes')
    Interpolation_3D
end


clearvars -except Ek kx_matrix ky_matrix kz_matrix kx_array ky_array kz_array a b c alat material_name
save_filename=['..\Bands\Ek_',material_name,'.mat'];   
save(save_filename) 