import sys
import re
import yaml
from pathlib import Path
import pandas as pd
import numpy as np
from pathlib import Path
from tifffile import imread
from tifffile import imwrite



def buildTrackObjects(segmentation, image):
    """A function to build btrack objects """
    PROPERTIES = ('area', )
    
    return btrack.utils.segmentation_to_objects(segmentation, 
                                                intensity_image=image, 
                                                use_weighted_centroid=False,
                                                properties=PROPERTIES)



def tracks_to_dataframe(tracks):
    """A function to convert tracklets to dataframe """
    return pd.concat([pd.DataFrame(x.to_dict()) for x in tracks])


def track(objects, trackConfigPath):
    
    FEATURES = [
        "area",
    ]

    TRACKING_UPDATES = [
        "motion",
        "visual",
    ]
    
    tracker = btrack.BayesianTracker()
    tracker.configure(trackConfigPath)
    
    # set up the features to use as a list
    tracker.features = FEATURES

    # append the objects to be tracked
    tracker.append(objects)

    # perform tracking
    tracker.track()
    tracker.optimize()

    # get tracks
    tracks = tracker.tracks
    return tracks_to_dataframe(tracks)
    



if __name__ == "__main__":
    btrackPath = sys.argv[1]
    imgPath = sys.argv[2]
    segPath = sys.argv[3]
    trackConfigPath = sys.argv[4] # tracker config
    outpath = sys.argv[5] # tracker config
    
    # add btrack to path
    sys.path.append(btrackPath)
    import btrack
    
    # load images
    segmentation = imread(segPath)
    img = imread(imgPath)
    
    # # reformat img to channel last
    img = np.moveaxis(img, 0, 3)
    
    # build objects
    objects = buildTrackObjects(segmentation, img)
    
    # track objects and save
    tracks = track(objects, trackConfigPath)
    tracks.to_csv(outpath, index=False)
    
    
    
    
    
    
    