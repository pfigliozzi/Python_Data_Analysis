from skimage.viewer import CollectionViewer
from skimage.draw import circle_perimeter
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import numpy as np
import pims
import glob

class ParticlePositionViewer(CollectionViewer):
    """A subclass of the skimage viewer that will display
    the position of particles on the images in a pims 
    object.

    Parameters
    ----------
    image_collection : ND image collection
        
    df_positions : DataFrame
        needs to contain columns ['frame', 'x',
        'y'].
    """
    def __init__(self, image_collection, df_positions, pos_columns=['x','y']):
        CollectionViewer.__init__(self, image_collection)
        self.df_positions = df_positions
        self.artists = None
        self.circle_kwargs = dict(edgecolor='r', radius=6, facecolor=None, fill=False, lw=1)
        
        self.pos_columns = pos_columns

        image = self.image_collection[0]
        
        # Use the min/max of the first image for all other images
        self.vmin_vmax = (np.min(image), np.max(image))
        
        self.update_image(image)
        
        
        
    def update_image(self, image):
        """Update displayed image.
        This method was overwritten in the skimage viewer. This method will
        plot the particle positions on the image.
        """
        if self.artists != None:
            self.artists.remove()
        image_df = self.df_positions[self.df_positions['frame']==self.index]

        # Add the circles over each particle
        patches = [Circle((v[self.pos_columns[0]], v[self.pos_columns[1]]), **self.circle_kwargs) for idx,v in image_df.iterrows()]
        collection = PatchCollection(patches, color='r', facecolor='none')
        self.artists = self.ax.add_collection(collection)
        
        
        self._update_original_image(image)
        self._image_plot.set_clim(self.vmin_vmax)


def annotate(df_positions, pims_obj):
    """Function to display the localized positon of particles in
    an image viewer.
    """
    viewer = ParticlePositionViewer(pims_obj, df_positions)
    viewer.show()

def annotate_from_path(df_positions, image_path, **kwargs):
    """Function to display the localized positon of particles in
    an image viewer by loading images from a specific path.
    """
    image_filenames = glob.glob(image_path+"\*.tif")
    pims_seq = pims.ImageSequence(image_filenames)
    print kwargs
    viewer = ParticlePositionViewer(pims_seq, df_positions, **kwargs)
    viewer.show()

class TrajectoryViewer(CollectionViewer):
    """A subclass of the skimage viewer that will display
    the trajectories of particles on the images in a pims 
    object.

    Parameters
    ----------
    image_collection : ND image collection
        
    df_positions : DataFrame
        needs to contain columns ['frame', 'x',
        'y'].
    """
    def __init__(self, image_collection, df_positions, pos_columns=['x','y'], tail_length=3, track_column='track'):
        CollectionViewer.__init__(self, image_collection)
        
        self.artists = None
        self.circle_kwargs = dict(edgecolor='r', radius=6, facecolor=None, fill=False, lw=1)
        
        self.pos_columns = pos_columns
        self.track_column = track_column
        self.df_positions = df_positions.sort_values(['frame', self.track_column])

        image = self.image_collection[0]
        
        # Use the min/max of the first image for all other images
        self.vmin_vmax = (np.min(image), np.max(image))
        
        
        self.tail_length = tail_length
        

        self.traj_colors = ["#ffbfbf", "#ced936", "#3db6f2",
                            "#6c468c", "#f2553d", "#98b386",
                            "#566573", "#a336d9", "#7f4840",
                            "#00ff00", "#0057d9", "#ffbffb",
                            "#7f3300", "#79f279", "#bfd9ff",
                            "#f200c2", "#d97736", "#1d7328",
                            "#1d3473", "#73566d", "#e6cbac",
                            "#66ccaa", "#7999f2", "#73004d",
                            "#ffaa00", "#00a7b3", "#0000f2",
                            "#ff0066", "#73561d", "#bffbff",
                            "#0000e6", "#d96c98", "#b2a159",
                            "#004d73", "#3d0073", "#73000f"]

        self.traj_colors_len = len(self.traj_colors)
        
        self.update_image(image)
        
        
    def update_image(self, image):
        """Update displayed image.
        This method was overwritten in the skimage viewer. This method will
        plot the particle positions and tails showing particle trajectories
        on the image.
        """
        if self.artists != None:
            for artist in self.artists:
                artist.remove()

        upper_frame = self.df_positions['frame'] <= self.index
        lower_frame = self.df_positions['frame'] >= self.index - self.tail_length
        image_df = self.df_positions[upper_frame & lower_frame]

        # Add the circles over each particle
        self.artists = []
        for name, group in image_df.groupby(self.track_column):
            cur_pos = group[group['frame'] == self.index]
            # Only draw circles on particles in the frame
            if len(cur_pos) != 0:
                circle = Circle((cur_pos[self.pos_columns[0]], cur_pos[self.pos_columns[1]]), **self.circle_kwargs)
                self.artists.append(self.ax.add_patch(circle))
            color = self.traj_colors[int(name % self.traj_colors_len)]
            line = self.ax.plot(group[self.pos_columns[0]], group[self.pos_columns[1]], color)
            self.artists.append(line[0])

        self._update_original_image(image)
        self._image_plot.set_clim(self.vmin_vmax)
        # patches = [Circle((v[self.pos_columns[0]], v[self.pos_columns[1]]), **self.circle_kwargs) for idx,v in image_df.iterrows()]
        # collection = PatchCollection(patches, color='r', facecolor='none')
        # self.artists = self.ax.add_collection(collection)

def traj_view(df_positions, pims_obj, **kwargs):
    """Function to display the localized positon of particles in
    an image viewer.
    """
    viewer = TrajectoryViewer(pims_obj, df_positions, **kwargs)
    viewer.show()

def traj_view_from_path(df_positions, image_path, **kwargs):
    """Function to display the localized positon of particles in
    an image viewer.
    """
    image_filenames = glob.glob(image_path+"\*.tif")
    pims_seq = pims.ImageSequence(image_filenames)
    viewer = TrajectoryViewer(pims_seq, df_positions, **kwargs)
    viewer.show()

# class CompareParticlePositionViewer(CollectionViewer):
#     """A subclass of the skimage viewer that will display
#     the position of particles on the images in a pims 
#     object.

#     Parameters
#     ----------
#     image_collection : ND image collection
        
#     df_positions : DataFrame
#         needs to contain columns ['frame', 'x',
#         'y'].
#     """
#     def __init__(self, image_collection, df_list, pos_columns=['x','y']):
#         CollectionViewer.__init__(self, image_collection)
#         self.df_list = df_list
#         self.artists = None
#         self.circle_kwargs = dict(edgecolor='r', radius=6, facecolor=None, fill=False, lw=1)
        
#         self.pos_columns = pos_columns

#         image = self.image_collection[0]
        
#         # Use the min/max of the first image for all other images
#         self.vmin_vmax = (np.min(image), np.max(image))
        
        
        
#         self.traj_colors = ["#ffbfbf", "#ced936", "#3db6f2",
#                             "#6c468c", "#f2553d", "#98b386",
#                             "#566573", "#a336d9", "#7f4840",
#                             "#00ff00", "#0057d9", "#ffbffb",
#                             "#7f3300", "#79f279", "#bfd9ff",
#                             "#f200c2", "#d97736", "#1d7328",
#                             "#1d3473", "#73566d", "#e6cbac",
#                             "#66ccaa", "#7999f2", "#73004d",
#                             "#ffaa00", "#00a7b3", "#0000f2",
#                             "#ff0066", "#73561d", "#bffbff",
#                             "#0000e6", "#d96c98", "#b2a159",
#                             "#004d73", "#3d0073", "#73000f"]

#         self.traj_colors_len = len(self.traj_colors)

#         self.update_image(image)
        
#     def update_image(self, image):
#         """Update displayed image.
#         This method was overwritten in the skimage viewer. This method will
#         plot the particle positions on the image.
#         """
#         if self.artists != None:
#             self.artists.remove()
#         image_df = self.df_positions[self.df_positions['frame']==self.index]

#         # Add the circles over each particle
#         patches = [Circle((v[self.pos_columns[0]], v[self.pos_columns[1]]), **self.circle_kwargs) for idx,v in image_df.iterrows()]
#         collection = PatchCollection(patches, color='r', facecolor='none')
#         self.artists = self.ax.add_collection(collection)
        
        
#         self._update_original_image(image)
#         self._image_plot.set_clim(self.vmin_vmax)


# def annotate(df_positions, pims_obj):
#     """Function to display the localized positon of particles in
#     an image viewer.
#     """
#     viewer = ParticlePositionViewer(pims_obj, df_positions)
#     viewer.show()

# def annotate_from_path(df_positions, image_path, **kwargs):
#     """Function to display the localized positon of particles in
#     an image viewer.
#     """
#     image_filenames = glob.glob(image_path+"\*.tif")
#     pims_seq = pims.ImageSequence(image_filenames)
#     print kwargs
#     viewer = ParticlePositionViewer(pims_seq, df_positions, **kwargs)
#     viewer.show()
