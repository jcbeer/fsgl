%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fused Sparse Group Lasso ABIDE Application
% Assign each voxel to the independent component of maximal value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script used for analyses reported in the manuscript
% "Incorporating Prior Information with Fused Sparse Group Lasso:
% Application to Prediction of Clinical Measures from Neuroimages"
%%% INPUTS: 
% melodic_IC_3mm.nii
%%% OUTPUTS:
% ICgroups.nii
% ICgroups.csv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% go to data directory
cd('./data')

% read the resing state ICA maps (from Cerliani et al 2015)
V=spm_vol('/Users/beerj2/Box Sync/Box Sync/ABIDE/data/melodic_IC_3mm.nii');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we will use 19 ICs:
% 1, 3, 5, 8, 9
% 10, 13, 15, 16, 17
% 19, 21, 23, 24, 25
% 27, 29, 30, 33
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save each IC map as a vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1, 3, 5, 8, 9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% volume 1
IC01=spm_read_vols(V(1));
IC01vect = reshape(IC01, [1,271633]);
sum(IC01vect==0)
clear IC01
% volume 3
IC03=spm_read_vols(V(3));
IC03vect = reshape(IC03, [1,271633]);
sum(IC03vect==0)
clear IC03
% volume 5
IC05=spm_read_vols(V(5));
IC05vect = reshape(IC05, [1,271633]);
sum(IC05vect==0)
clear IC05
% volume 8
IC08=spm_read_vols(V(8));
IC08vect = reshape(IC08, [1,271633]);
sum(IC08vect==0)
clear IC08
% volume 9
IC09=spm_read_vols(V(9));
IC09vect = reshape(IC09, [1,271633]);
sum(IC09vect==0)
clear IC09

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10 13 15 16 17 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% volume 10
IC10=spm_read_vols(V(10));
IC10vect = reshape(IC10, [1,271633]);
sum(IC10vect==0)
clear IC10
% volume 13
IC13=spm_read_vols(V(13));
IC13vect = reshape(IC13, [1,271633]);
sum(IC13vect==0)
clear IC13
% volume 15
IC15=spm_read_vols(V(15));
IC15vect = reshape(IC15, [1,271633]);
sum(IC15vect==0)
clear IC15
% volume 16
IC16=spm_read_vols(V(16));
IC16vect = reshape(IC16, [1,271633]);
sum(IC16vect==0)
clear IC16
% volume 17
IC17=spm_read_vols(V(17));
IC17vect = reshape(IC17, [1,271633]);
sum(IC17vect==0)
clear IC17

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 19 21 23 24 25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% volume 19
IC19=spm_read_vols(V(19));
IC19vect = reshape(IC19, [1,271633]);
sum(IC19vect==0)
clear IC19
% volume 21
IC21=spm_read_vols(V(21));
IC21vect = reshape(IC21, [1,271633]);
sum(IC21vect==0)
clear IC21
% volume 23
IC23=spm_read_vols(V(23));
IC23vect = reshape(IC23, [1,271633]);
sum(IC23vect==0)
clear IC23
% volume 24
IC24=spm_read_vols(V(24));
IC24vect = reshape(IC24, [1,271633]);
sum(IC24vect==0)
clear IC24
% volume 25
IC25=spm_read_vols(V(25));
IC25vect = reshape(IC25, [1,271633]);
sum(IC25vect==0)
clear IC25

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 27 29 30 33
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% volume 27
IC27=spm_read_vols(V(27));
IC27vect = reshape(IC27, [1,271633]);
sum(IC27vect==0)
clear IC27
% volume 29
IC29=spm_read_vols(V(29));
IC29vect = reshape(IC29, [1,271633]);
sum(IC29vect==0)
clear IC29
% volume 30
IC30=spm_read_vols(V(30));
IC30vect = reshape(IC30, [1,271633]);
sum(IC30vect==0)
clear IC30
% volume 33
IC33=spm_read_vols(V(33));
IC33vect = reshape(IC33, [1,271633]);
sum(IC33vect==0)
clear IC33

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make a component that is all zeros 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IC00vect=zeros([1 271633]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% combine all IC vectors into a matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
icmatrix = [ IC00vect;
    IC01vect ; IC03vect ; IC05vect ; IC08vect ; IC09vect ; 
    IC10vect ; IC13vect ; IC15vect ; IC16vect ; IC17vect ; 
    IC19vect ; IC21vect ; IC23vect ; IC24vect ; IC25vect ; 
    IC27vect ; IC29vect ; IC30vect ; IC33vect];
clear IC00vect IC01vect IC03vect IC05vect IC08vect IC09vect IC10vect IC13vect IC15vect IC16vect IC17vect IC19vect IC21vect IC23vect IC24vect IC25vect IC27vect IC29vect IC30vect IC33vect

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get maximum component for each column
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxic = max(icmatrix);
sum(maxic==0)
[maxic3 idx] = max(icmatrix,[],1);
[row, col] = ind2sub(size(icmatrix), idx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tabulate IC groups
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: matching this coding of ICs 
% to the coding in Cerliani et al 2015
% this code     Cerliani et al code
%   1           None -- outside of brain
%   2           IC1           
%   3           IC3
%   4           IC5 -- dorsal sensorimotor
%   5           IC8
%   6           IC9
%   7           IC10
%   8           IC13
%   9           IC15
%   10          IC16
%   11          IC17 -- basal ganglia / thalamus
%   12          IC19
%   13          IC21
%   14          IC23
%   15          IC24
%   16          IC25
%   17          IC27
%   18          IC29 -- ventral sensorimotor
%   19          IC30
%   20          IC33
tabulate(row)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save a NIFTI file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reshape to make a volume
ICgpvol = reshape(row, [61 73 61]);
% take the header information from a previous file with similar dimensions 
% and voxel sizes and change the filename in the header
HeaderInfo=spm_vol('/data/MNI152maskbin_3mm.nii')
% fill in the new filename
HeaderInfo.fname = 'ICgroups.nii';  
% replace the old filename in another location within the header
HeaderInfo.private.dat.fname = HeaderInfo.fname;  
% use spm_write_vol to write out the new data
% give spm_write_vol the new header information and corresponding data matrix
% HeaderInfo is your header information for the new file
% ICgpvol is the image matrix corresponding to the image you'll be writing out
spm_write_vol(HeaderInfo,ICgpvol);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save a csv file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subtract one so that voxels not in mask are in group zero
icgps = row - 1;
% save group assignments to a csv file
csvwrite('ICgroups.csv', icgps);