'''
function [matched,unmatched] = corr_mols(set1_pos,set2_pos,tform_start,match_radius)

% This function takes two input sets of positions and matches the points to
% each other based on the initial transform tform_start and the radial
% tolerance match_radius. set1 is the master set, and set2 is matched to
% it. Both input variables set1_pos and set2_pos must be structs with
% length 1. set1_pos and set2_pos should have fields .x and .y that contain
% vectors of the x and y positions. The input tform_start is the spatial
% transformation matrix (see cp2tform). The output matched and unmatched
% each have fields set1_inds and set2_inds.
%
% SYNTAX: [matched,unmatched] = corr_mols(set1_pos,set2_pos,tform_start,match_radius)

% initialize variables
num_mols_set1 = length(set1_pos.x);
matching_set2_inds = zeros(1,num_mols_set1);

% overlap set2 onto set1 (shift, or also warp, depending on tform_start)
[registered_set2_pos.x,registered_set2_pos.y] = tforminv(tform_start,set2_pos.x',set2_pos.y');

% correlate molecules in the two sets
for n = 1:num_mols_set1 % loop over each set1 index
    xdiff = set1_pos.x(n) - registered_set2_pos.x;
    ydiff = set1_pos.y(n) - registered_set2_pos.y;
    matching_set2_ind = find(sqrt(xdiff.^2 + ydiff.^2)<match_radius);
    if length(matching_set2_ind)==1  % no more/less than 1 matching molecule within match radius
        matching_set2_inds(n) = matching_set2_ind;
    end
end

% get only the matching molecules
matched.set1_inds = find(matching_set2_inds>0);
matched.set2_inds = matching_set2_inds(matched.set1_inds);
unmatched.set1_inds = find(matching_set2_inds==0);
unmatched.set2_inds = setdiff(1:length(set2_pos),matched.set2_inds);
'''
import numpy as np

def corr_mols(set1_pos_x, set1_pos_y, set2_pos_x, set2_pos_y, tform_start, match_radius): 
    
    #% initialize variables
    num_mols_set1 = len(set1_pos_x);
    
    #matching_set2_inds = np.zeros(num_mols_set1);
    matching_set2_inds=np.zeros((num_mols_set1,), dtype=np.int)-1
    
    #% overlap set2 onto set1 (shift, or also warp, depending on tform_start)
    
    '''
    set1_pos = np.column_stack([set1_pos_x, set1_pos_y]);
    registered_set2_pos_x = (tform_start_inverse(tform_start(set1_pos)))[:,0];
    registered_set2_pos_y = (tform_start_inverse(tform_start(set1_pos)))[:,1]; 
    '''
    
    set2_pos = np.column_stack([set2_pos_x, set2_pos_y]);
    registered_set2_pos_x = (tform_start.inverse(set2_pos))[:,0];
    registered_set2_pos_y = (tform_start.inverse(set2_pos))[:,1];   
    
    #% correlate molecules in the two sets
    for n in range(0, num_mols_set1): # % loop over each set1 index
        xdiff = set1_pos_x[n] - registered_set2_pos_x;
        ydiff = set1_pos_y[n] - registered_set2_pos_y;
        
        matching_set2_ind, = np.where((np.sqrt(np.square(xdiff) + np.square(ydiff)) < match_radius));
        
        if len(matching_set2_ind)==1:  #% no more/less than 1 matching molecule within match radius
            matching_set2_inds[n] = matching_set2_ind[0];
            
    #% get only the matching molecules
    matched_set1_inds, = np.where(matching_set2_inds != -1);    
    matched_set2_inds = matching_set2_inds[matched_set1_inds];
        
    #unmatched_set1_inds, = np.where(matching_set2_inds==0);
    
    unmatched_set1_inds, = np.where(matching_set2_inds==-1);
    
    #put start=0 and stop=1 to match results from matlab code. 1:length(set2_pos) is the same as 1:1 in matlab, not same in python
    temp_arr = np.arange(start=0, stop=len(set2_pos), step=1);
    #temp_arr = np.arange(start=0, stop=1, step=1);
    
    unmatched_set2_inds = np.setdiff1d(temp_arr, matched_set2_inds);
    
    return  matched_set1_inds, matched_set2_inds, unmatched_set1_inds, unmatched_set2_inds;