from skimage.viewer import CollectionViewer
from skimage.draw import circle_perimeter
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import numpy as np
import pims

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
    def __init__(self, image_collection, df_positions):
        CollectionViewer.__init__(self, image_collection)
        self.df_positions = df_positions
        self.artists = None
        self.circle_kwargs = dict(edgecolor='r', radius=6, facecolor=None, fill=False, lw=1)
        
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
        patches = [Circle((v['x'], v['y']), **self.circle_kwargs) for idx,v in image_df.iterrows()]
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

def annotate_from_path(df_positions, image_path):
	"""Function to display the localized positon of particles in
	an image viewer.
	"""
	image_filenames = glob.glob(image_path+"\*.tif")
	pims_seq = pims.ImageSequence(image_files)
	viewer = ParticlePositionViewer(pims_seq, df_positions)
	viewer.show()
