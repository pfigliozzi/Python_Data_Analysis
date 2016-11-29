import scipy
import numpy as np
import pandas as pd


def check_full_rotation(group):
    '''Returns the indices where a particle has gone from one side of a
    periodic boundary to the other. For example, a trajectory where a 
    particle goes from theta = 355 to theta = 2 degrees.
    
    :param group: The group of one track from the groupby operation, should
    have no duplicate entries w.r.t. frame and track id. Must have a boolean
    column 'start_region' and 'end_region'.
    '''
    # Find indicies where particle just came from end_region
    new_group = group.iloc[:-1]
    cur_bool = group.end_region.iloc[:-1]
    fut_bool = group.start_region.shift(-1).iloc[:-1].astype(np.bool)
    from_end_region_frwrd = new_group[(cur_bool & fut_bool)].index
    
    # Find indicies where particle just went to start_region
    new_group = group.iloc[1:]
    cur_bool = group.start_region.iloc[1:]
    past_bool = group.end_region.shift(1).iloc[1:].astype(np.bool)
    to_start_region_frwrd = new_group[(cur_bool & past_bool)].index
    
    # Find indicies where particle just came from end_region
    new_group = group.iloc[:-1]
    cur_bool = group.start_region.iloc[:-1]
    fut_bool = group.end_region.shift(-1).iloc[:-1].astype(np.bool)
    from_end_region_rev = new_group[(cur_bool & fut_bool)].index
    
    # Find indicies where particle just went to start_region
    new_group = group.iloc[1:]
    cur_bool = group.end_region.iloc[1:]
    past_bool = group.start_region.shift(1).iloc[1:].astype(np.bool)
    to_start_region_rev = new_group[(cur_bool & past_bool)].index
    
    #print from_end_region_frwrd, to_start_region_frwrd, from_end_region_rev, to_start_region_rev
    
    indicies = np.sort(np.concatenate([from_end_region_frwrd, to_start_region_frwrd, from_end_region_rev, to_start_region_rev]))
    
    #print group.loc[indicies]
    
    return group.loc[indicies].index

def find_full_trajs_around_ring(df, periodic_limit=360):
    '''Returns a DataFrame of trajectories where trajectories that traverse
    the entirety of the periodic region are labeled.
    
    This function accepts a DataFrame of rotational trajectories and labels
    the trajectories that traverse through the entirety of the periodic region.
    This function is not really robust and only has been tested with the single
    particle L ramping over the nanoplate. Trajectories are not necessarily
    consecutive. It works for both positive and negative L's.
    
    :param df: DataFrame of trajectories to be analyzed. Does not need to have
    unique values w.r.t. 'frame' and 'track id'.
    :param periodic_limit: The upper limit for the cut off of the periodic 
    boundary condition. The lower limit is assumed to be 0.
    '''
    
    # Designate if particle positions are in the start region with 'start_region" column
    df_unique = df.copy().drop_duplicates(subset=['frame', 'track id'], keep='first')
    df_in_region = (0 < df_unique.theta) & (df_unique.theta < periodic_limit * 0.3)
    df_unique.loc[df_in_region, 'start_region'] = True
    df_unique.loc[~df_in_region, 'start_region'] = False
    
    # Designate if particle positions are in the end region with 'end_region" column
    df_in_region = (periodic_limit * 0.6 < df_unique.theta) & (df_unique.theta < periodic_limit)
    df_unique.loc[df_in_region, 'end_region'] = True
    df_unique.loc[~df_in_region, 'end_region'] = False
    
    # Find positions where particles just enter or exit the region
    full_traj_indices = df_unique.groupby('track id', group_keys=False).apply(check_full_rotation)
    
    #print full_traj_indices
    
    df_unique['full_track_region'] = False
    df_unique['full_track_id'] = np.nan
    for track in full_traj_indices:
        for num in range(0, len(track), 2):
            num += 1
            if num == len(track)-1:
                continue
            start = track[num]
            end = track[num+1]
            df_unique.loc[start:end,'full_track_id'] = int((num - 1)/2)
            df_unique.loc[start:end, 'full_track_region'] = True
    
    return df_unique.drop(['start_region', 'end_region'], axis=1)