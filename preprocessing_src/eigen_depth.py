import os
import sys
from numpy import *

from PIL import Image
import ipdb
ipdb.launch_ipdb_on_exception()

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

    imgs = frames.tolist() # this originated from a matlab cell array. convert it to a list
    Nimgs = len(imgs[0][0][0])
    # print 'Nimgs = {Nimgs}\n'.format(**locals())
    est_depth_imgs, est_normals_imgs = [None]*Nimgs, [None]*Nimgs

    for k in range(Nimgs):
        img  = frames[0][0][0][k].reshape(1,240,320,3)*255
        est_depth_imgs[k], est_normals_imgs[k] = infer_depth_and_normals(nn, img.astype(uint8))

    # est_depth_imgs, est_normals_imgs = infer_depth_and_normals(nn, imgs)
    return est_depth_imgs, est_normals_imgs

        
