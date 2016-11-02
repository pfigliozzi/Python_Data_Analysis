import scipy.io as sio

def add_axes_label_inches(ax, (right_left, down_up), string, corner='upper left', **kwargs):
    """A helper function to add a text label a specific distance from an axes corner in inches.
    
    This function was created so that when constructing multi-panel figures for a manuscript all the axes
    labels (e.g. a, b, c, d...) are a consistent distance from the corner of the axis. The relative
    coordinates used for text labels can result in a label being a different physical distance compared
    to subplots of different sizes.
    
    :param ax: A mpl.axes object that you want to add text to
    :param (right_left, down_up): Distance in inches you want the text label to be from the selected
    axes corner
    :param string: The text string you want displayed at that location.
    :param corner: The corner you want to have the text origin to be from the axes corner
    """
    fig = ax.get_figure()
    fig_size = fig.get_size_inches()
    ax_bbox = ax.get_position()
    ax_rect_inches = ax_bbox.x0*fig_size[0], ax_bbox.y0*fig_size[1], ax_bbox.x1*fig_size[0], ax_bbox.y1*fig_size[1]
    if corner == 'upper left':
        text_location_inches = [right_left, ax_rect_inches[3]-ax_rect_inches[1]-down_up]
        va = 'top'
        ha = 'left'
    if corner == 'upper right':
        text_location_inches = [ax_rect_inches[2]-ax_rect_inches[0] - right_left, ax_rect_inches[3]-ax_rect_inches[1]-down_up]
        va = 'top'
        ha = 'right'
    if corner == 'lower left':
        text_location_inches = [right_left, down_up]
        va = 'bottom'
        ha = 'left'
    if corner == 'lower right':
        text_location_inches = [ax_rect_inches[2]-ax_rect_inches[0] - right_left, down_up]
        va = 'bottom'
        ha = 'right'
    text_position_rel_coors = text_location_inches[0]/(ax_rect_inches[2]-ax_rect_inches[0]), text_location_inches[1]/(ax_rect_inches[3]-ax_rect_inches[1])
    return ax.text(text_position_rel_coors[0], text_position_rel_coors[1], string, transform=ax.transAxes, va=va, ha=ha, **kwargs)

def make_cbar_match_polar_axis_height(polar_ax, cbar_ax):
    """Makes a color bar that corresponds to a polar projection plot
    the same height as the polar plot.
    
    This function works by converting data coordinates to of the top 
    and bottom of the polar plot into figure coordinates. Once figure 
    coordinates are discovered for the polar plot the colorbar's height is
    adjusted to match the polar plot. This function only works if your x
    axis (theta axis) goes from 0 to 2pi. The key to get this to work is
    forcing matplotlib to draw the polar plot once so that the data 
    coordinates can be properly converted into figure coordinates.
    
    :param polar_ax: The mpl.axes object that contains the polar plot
    :param cbar_ax: The mpl.axes object that contains the color bar
    """
    plt.draw()
    fig = polar_ax.get_figure()
    ylim_max = polar_ax.get_ylim()[1]
    theta_offset = ax.get_theta_offset()
    
    x_90 = (np.pi/2) - theta_offset
    
    top_point_disp = polar_ax.transData.transform_point([x_90, ylim_max])
    top_point_fig = fig.transFigure.inverted().transform(top_point_disp)
    
    bot_point_disp =  polar_ax.transData.transform_point([x_90+np.pi, ylim_max])
    bot_point_fig = fig.transFigure.inverted().transform(bot_point_disp)
    
    cbar_pos = cbar_ax.ax.get_position()
    cbar_ax.ax.set_position([cbar_pos.x0, bot_point_fig[1],
                          cbar_pos.x1 - cbar_pos.x0,
