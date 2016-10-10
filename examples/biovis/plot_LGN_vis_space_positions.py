import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
matplotlib.rcParams.update({'font.size': 20})

type_colors = { 'transient_ON': 'green', 'transient_OFF': 'magenta', 'transient_ON_OFF': 'cyan' }


map_df = pd.read_csv('../build/ll2_inputs_from_LGN.csv', sep=' ')
gids = np.array(list(set(map_df['index'].values)))

filters_data_fname = '../../6-LGN_firing_rates_and_spikes/LGN_spike_trains/LGN2_visual_space_positions_and_cell_types.dat'

'''
# Find cells with a relatively large number of ON_OFF inputs.
ON_OFF_3 = []
for gid in gids:
    if (gid >= 8500):
        break
    if (gid % 100 == 0):
        print 'Processing gid %d.' % (gid)
    tmp_df = map_df.loc[map_df['index'] == gid]
    N_ON_OFF = len(tmp_df.loc[tmp_df['src_type'] == 'transient_ON_OFF'].index)
    if (N_ON_OFF >= 3):
        ON_OFF_3.append(gid)
print len(ON_OFF_3)
'''

# Plot positions of the inout LGN filters (in the visual space) for the selected cell.
filters_data = pd.read_csv(filters_data_fname, sep=' ')
filters_data.columns = ['LGN_type', 'x', 'y', 'x_offset', 'y_offset', 'sigma_c', 'sigma_s', 'r0', 'scaling_factor', 'k_alpha']
for cell_gid in [107]: #[107, 1232, 2102, 4621, 4731, 7192]: #ON_OFF_3: #[300, 690]: #xrange(0, 8500, 10): #Cells with 5 ON_OFF inputs: [622, 653, 7252]
    tmp_df = map_df.loc[map_df['index'] == cell_gid]
    print 'N_LGN_filters =', len(tmp_df.index)

    # Print out the parameters of LGN source filters associated with the cell_gid.
    src_gids = tmp_df['src_gid'].values
    print filters_data.ix[src_gids] # Use src_gids as indices in the data frame; gids are not saved in the original file, but entries are saved in order of gids.

    src_types = list(set(tmp_df['src_type'].values))
    for type in src_types:
        x = tmp_df[tmp_df['src_type'] == type]['src_vis_x'].values
        y = tmp_df[tmp_df['src_type'] == type]['src_vis_y'].values
        plt.scatter(x, y, s=500, c=type_colors[type], label=type)
    plt.gca().set_aspect('equal')
    plt.legend()
    plt.title('Cell gid %d' % (cell_gid))
    plt.savefig('plot_LGN_vis_space_positions/ll2_LGN_src.eps', format='eps')
    plt.show()

