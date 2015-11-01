# depth_video
inference of depth from a video sequence

# Installation instructions for using theano with a GPU

virtualenv --system-site-packages -p python2.7 theano-env
source theano-env/bin/activate
pip install theano
pip install ipdb

# Add to .bashrc (or just set it manually)
export PATH="/usr/local/cuda-7.5/bin":$PATH
export CUDA_HOME="/usr/local/cuda-7.5"
export LD_LIBRARY_PATH="/usr/local/cuda-7.5/lib64"


# Run eigen depth normals demo (from where it is installed. e.g. ~/externals/dnl-depthnormals$
THEANO_FLAGS=floatX=float32,device=gpu0 python demo_depthnormals.py
