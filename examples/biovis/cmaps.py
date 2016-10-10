# color maps:
import numpy as np

#*** BASIC COLOUR DICTIONARY ******
RED = 255, 0, 0
ORANGE = 255, 127, 0
YELLOW = 255, 255, 0
GREEN = 0, 255, 0
CYAN = 0, 255, 255
LIGHTBLUE = 0, 127, 255
BLUE = 0, 0, 255
VIOLET = 127, 0, 255
MAGENTA = 255, 0, 255  
GREY = 85, 85, 85
WHITE = 255, 255, 255
DARK_GREY = 30, 30, 30
PURPLE = 154, 44, 209
CMAP_BASIC = [RED, ORANGE, YELLOW, GREEN, CYAN, LIGHTBLUE, BLUE, VIOLET, MAGENTA, GREY, WHITE, PURPLE]



from matplotlib import colors
import six



def get_named_colors():

    colors_ = list(six.iteritems(colors.cnames))

    return colors_


def construct_colors():

    '''
    Returns
    -------
    color_lib: dict
               color_name: (r,g,b) for all available named colors                
    '''
    colors_ = list(six.iteritems(colors.cnames))
    colors_tup = [(color[0],colors.hex2color(color[1])) for color in colors_]
    colors_dict = dict(colors_tup) # name-to-rgb conversion
    colors_lib = {}

    for col,rgb_float in colors_dict.items():
       rgb_int = int(255*rgb_float[0]),int(255*rgb_float[1]),int(255*rgb_float[2])
       colors_lib[col]=rgb_int
    
    return colors_lib


def convert_to_rgb(mymap): 
    '''
    Convert names from a specific colormap to rgb expressed using 0..255

    Parameters:
    ----------
    mymap: dict
           label:color_name
    Returns
    -------
    color_lib: dict
              label:(r,g,b)
    '''

    cmap_rgb = {}
    for label,col_name in mymap.items(): # convert my color_names to rgb
       cmap_rgb[label] = colors_lib[col_name]

    return cmap_rgb


def make_color_map(labels,colors):

    '''
    Make map for a list of labels given a list of colors. This helps when you do not want to manually define the colormap
    Parameters:
    ----------
    labels: np.array
           unique labels from a particular column
    colors: list 
           to associate with the labels
    Returns
    -------
    cmap: dict
          label:(r,g,b)
    '''

    cmap = {}

    for color,label in zip(colors,labels):
        cmap[label] = color
    return cmap
 
# run on import to have the library of named colors
colors_lib = construct_colors()


