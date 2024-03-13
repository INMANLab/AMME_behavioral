%Script to estimate the distance/relationship between amygdala
%stimulation and hippocampal structure.

%% fetch pertinent data:

%define directories:
project_dir='/Volumes/groups/WillieLabPatientData/AMME_BIDS'; %path to AMME data directory
figures_dir='~/Desktop'; %path to output directory
sessions_list=[project_dir '/stg-preproc/group-data/AMME_electrodeStimContacts.xlsx']; %table with list of sessions and stim electrode labels for each session
fs_lookuptable_fname='/Volumes/groups/WillieLabPatientData/ElectrodeLocalization_Anat/lookup-tables/FreeSurferColorLUT_v7-1-0.txt'; %path to freesurfer labels (lut) lookup table
ashs_lookuptable_fname = '.txt';%path to the ashs look up table file

%get freesurfer label-to-name correspondance:
ashs_lookuptable={};
fileID=fopen(ashs_lookuptable_fname,'r');
while ~feof(fileID) %while we have not reached end of file
    current_line=fgetl(fileID); %get next line
    if ~isempty(current_line) && ~strcmp(current_line(1),'#')
        ashs_lookuptable(end+1,:)=textscan(current_line,'%d%s%d%d%d%d');
    end
end
fclose(fileID);

%fetch stim electrode information for all sessions, and clean the table:
stim_electrodes_table=readtable(sessions_list,'ReadRowNames',true,'PreserveVariableNames',true); %fetch list of stim session and the stim electrode labels
stim_electrodes_table.Properties.RowNames=lower(stim_electrodes_table.Properties.RowNames); %subject IDs in lower case
stim_electrodes_table=stim_electrodes_table(:,contains(stim_electrodes_table.Properties.VariableNames,{'cathode','anode','d''','task','test_day'})); %only keep columns we need
if ~all(startsWith(stim_electrodes_table.Properties.RowNames,'amyg0'))
    error('Not all rows start with ''amyg0''');
end
stim_electrodes_table(isnan(stim_electrodes_table.('d''')),:)=[]; %remove sessions where we don't have d'
stim_electrodes_table(contains(stim_electrodes_table.Properties.RowNames,{'hc','ecx','cing'}),:)=[]; %remove sessions where stimulate hippocampal region or cingulate
[row,col]=find(ismember(table2array(stim_electrodes_table(:,1:end-3)),{'ground','x'})); %find rows that have these labels, and replace with empty cell (those indicate monopolar stimulation)
disp(['Total sessions with monopolar stimulation: ' num2str(length(row))]);
for i=1:length(row)
    stim_electrodes_table(row(i),col(i))={''};
end
stim_electrodes_table(stim_electrodes_table.test_day>1,:)=[]; %remove sessions for which we did 1 day test more than 1 day later
stim_electrodes_table(sum(arrayfun(@(x) length(x{:})~=0, table2array(stim_electrodes_table(:,1:end-3)))')==4,:)=[];

%define list of subject IDs:
all_subjectIDs=stim_electrodes_table.Properties.RowNames;
all_subjectIDs=regexp(all_subjectIDs,'amyg...','match'); %keep only the subject ID (not anything appended to it)
all_subjectIDs=[all_subjectIDs{:}];
all_subjectIDs=unique(all_subjectIDs); %only keep unique subject IDs

%fetch anatomical data for all subjects:
subject=struct(); %initialize structure that will keep track of anatomical information for each subject
subject_nber=0;
for subjectID=all_subjectIDs
    subject_nber=subject_nber+1;
    subject(subject_nber).subjectID=subjectID{:}; %define subjectID
    
    %import electrodes information, including transformation matrices
    subject(subject_nber).anat.electrodes_info=importdata([project_dir '/stg-preproc/sub-' subject(subject_nber).subjectID '/anat/electrode-models/sub-' subject(subject_nber).subjectID '_electrodes_info.mat'],'electrodes_info');
    
    %fetch mni world coordinates (already calculated):
    subject(subject_nber).anat.fs_coords_mni=subject(subject_nber).anat.electrodes_info.ecoords_mni_world;
    subject(subject_nber).anat.ashs_coords_mni=subject(subject_nber).anat.electrodes_info.ecoords_mni_world;

    %calculate subject world coordinates:
    temp_electrodes_info=subject(subject_nber).anat.electrodes_info; %fetch relevant transformation matrices information for this subject
    temp_coords=temp_electrodes_info.voxcoords_postopct; %get post op ct vox coords
    temp_coords(:,4)=1;
    temp_coords=temp_coords';
    temp_coords=temp_electrodes_info.postopct.vox2ras*temp_coords; %convert to postop ct RAS space
    temp_coords=temp_electrodes_info.postopct.ras2ras*temp_coords; %convert to preop T1 MRI space
    temp_coords=temp_coords';
    temp_coords=temp_coords(:,1:3);
    subject(subject_nber).anat.coords_nat=temp_coords;
    
    %fetch freesurfer segmentation for that subject:
    subject(subject_nber).anat.fs_aparcaseg=ft_read_mri([project_dir '/stg-preproc/sub-' subject(subject_nber).subjectID '/anat/segmentation/sub-' subject(subject_nber).subjectID '_freesurfer-output/mri/aparc+aseg.mgz']);
    %NEED TO FIGURE OUT HOW TO FETCH FS_HIPPOAMYG SEGMENTATION
    
    %convert atlas voxels to coordinates (each coordinate is center of 1 voxel):
    atlas_volume=subject(subject_nber).anat.fs_aparcaseg;
    atlas_voxel_coords=zeros([numel(atlas_volume.anatomy) 4]);
    row_nber=0;
    for x=1:size(atlas_volume.anatomy,1)
        for y=1:size(atlas_volume.anatomy,2)
            for z=1:size(atlas_volume.anatomy,3)
                row_nber=row_nber+1;
                atlas_voxel_coords(row_nber,:)=[x-0.5 y-0.5 z-0.5 atlas_volume.anatomy(x,y,z)];
            end
        end
    end
    %convert atlas coordinates to RAS space (BELOW CORRECT?):
    temp=atlas_voxel_coords(:,1:3);
    temp(:,4)=1;
    temp=temp';
    temp=atlas_volume.transform*temp;
    temp=temp';
    temp=temp(:,1:3);
    atlas_voxel_coords(:,1:3)=temp;
    subject(subject_nber).anat.fs_aparcaseg_coords=atlas_voxel_coords;
    
end


%% create table and store relevant information:

struct_label='Hippocampus'; %structure of interest
struct_label=[fs_lookuptable{endsWith([fs_lookuptable{:,2}],struct_label),1}]; %find freesurfer labels for that structure of interest
results_table=cell([0 5]); %initialize results table
row_nber=0;
for ses=stim_electrodes_table.Properties.RowNames' %for each AMME session
    subjectID=regexprep(ses{:},'_.+',''); %get subject ID
    sub_nber=find(strcmp({subject(:).subjectID},subjectID));
    temp_stimelectrodes=stim_electrodes_table{ses,1:4};
    for temp=temp_stimelectrodes %for each stim electrode
        if ~isempty(temp{:}) %if we do have a stim electrode
            row_nber=row_nber+1; %increment results_table row number
            results_table{row_nber,1}=ses{:}; %get session ID
            results_table{row_nber,2}=temp{:}; %get electrode label
            results_table{row_nber,3}=subject(sub_nber).anat.coords_nat(strcmp(subject(sub_nber).anat.electrodes_info.labels,temp),:);
            
            %get closest distance to hippocampus:
            hippo_coords=subject(sub_nber).anat.fs_aparcaseg_coords(ismember(subject(sub_nber).anat.fs_aparcaseg_coords(:,4),struct_label),1:3); %get all RAS hippo coords
            [temp_dist,temp_index]=min(sum(((hippo_coords(:,1:3)-results_table{row_nber,3}).^2)'));
            results_table{row_nber,4}=hippo_coords(temp_index,1:3); %RAS coordinate of hippocampal voxel closest to electrode
            results_table{row_nber,5}=temp_dist;
        end
    end
end

%save results table to excel sheet:
results_table=array2table(results_table,'VariableNames',{'sessionID','electrode_label','electrode_nat_coords','closest_struct_coords','closest_struct_dist'});
writetable(results_table,[figures_dir '/EC_dist_results_table.csv']);


