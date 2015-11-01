import os
import sys
import numpy as np

from PIL import Image

import net
from utils import depth_montage, normals_montage

def load_depth_nn(model_name = 'depthnormals_nyud_vgg', nn_path):
    #model_name = 'depthnormals_nyud_alexnet' # an alternative model_name

    # add path of the neural network
    sys.path.append(nn_path)

    depth_nn = dict()

    # location of depth module, config and parameters
    depth_nn['model_name'] = model_name
    depth_nn['module_fn']  = 'models/iccv15/%s.py' % model_name
    depth_nn['config_fn'] = 'models/iccv15/%s.conf' % model_name
    depth_nn['params_dir'] = 'weights/iccv15/%s' % model_name
    

    # load depth network
    depth_nn['machine'] = net.create_machine( 
        depth_nn['module_fn'],  depth_nn['config_fn'],  depth_nn['params_dir'])

    return depth_nn

def infer_depth_and_normals(depth_nn, rgb_img)
    

    (pred_depths, pred_normals) = depth_nn['machine'].infer_depth_and_normals(rgb_img)

    return (pred_depths, pred_normals)

