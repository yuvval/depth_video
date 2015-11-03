import os
import sys
import numpy as np

from PIL import Image


def test(n = 5):
    print 'Test\n'

def load_depth_nn(model_name = 'depthnormals_nyud_vgg', nn_path = '/afs/csail.mit.edu/u/y/yuvval/externals/dnl-depthnormals/'):
    #model_name = 'depthnormals_nyud_alexnet' # an alternative model_name

    # add path of the neural network
    sys.path.append(nn_path)
    from utils import depth_montage, normals_montage
    import net
    depth_nn = dict()
    # location of depth module, config and parameters
    depth_nn['model_name'] = model_name
    depth_nn['module_fn']  = '%smodels/iccv15/%s.py' %(nn_path, model_name)
    depth_nn['config_fn'] = '%smodels/iccv15/%s.conf' %(nn_path, model_name)
    depth_nn['params_dir'] = '%sweights/iccv15/%s' %(nn_path, model_name)


    # load depth network
    depth_nn['machine'] = net.create_machine(
        depth_nn['module_fn'],  depth_nn['config_fn'],  depth_nn['params_dir'])

    return depth_nn

def infer_depth_and_normals(depth_nn, rgb_img):

    (pred_depths, pred_normals) = depth_nn['machine'].infer_depth_and_normals(rgb_img)

    return (pred_depths, pred_normals)


def infer_depth_and_normals_frames_seq(
        frames, 
        model_name = 'depthnormals_nyud_vgg', 
        nn_path = '/afs/csail.mit.edu/u/y/yuvval/externals/dnl-depthnormals/'):

    nn = load_depth_nn(model_name, nn_path)

    frames = frames.tolist() # this originated from a matlab cell array. convert it to a list
    Nframes = len(frames)

    est_depth_frames, est_normals_frames = [None]*Nframes, [None]*Nframes

    for k in range(Nframes):
        img  = np.asarray(frames[k]).reshape((1, 240, 320, 3))
        est_depth_frames[k], est_normals_frames[k] = infer_depth_and_normals(nn, img)

    est_depth_frames, est_normals_frames

        
