import numpy
from _tvdenoising import tvdenoise

def im3ddenoise(img,iter=50,lmbda=15.0):
    imgmin = min(img.ravel())
    img = img-imgmin
    imgmax = max(img.ravel())
    img = img/imgmax

    return tvdenoise(img,iter,lmbda)

