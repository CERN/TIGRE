from __future__ import division
import numpy as np
import tigre
import copy
from tigre.utilities.im_3d_denoise import im3ddenoise
from tigre.demos.Test_data import data_loader
import time

def fista(y, geo, angles, niter, trueimg):
    lmbda = 0.1
    y_rec = np.zeros((geo.nVoxel), dtype=np.float32)
    x_rec = copy.deepcopy(y_rec)
    t = 1
    #L = 2.e7
    L = 2.e8
    bm = 1 / L
    relativeError = []
    avgtime = []
    for i in range(niter):
        avgtic = time.clock()
        y_rec += bm * 2 * (tigre.Atb((y - tigre.Ax(y_rec, geo, angles, 'interpolated')),
                                     geo, angles, 'matched'))
        if i ==0:
            print('FISTA algorithm in progress.')
            toc = time.clock()
        if i ==1:
            tic = time.clock()
            print('Estimated time until completion: ' +str(niter*(tic-toc))+ '(s)')

        lambdaForTv = 2 * bm * lmbda

        #y_rec = tvdenoise(y_rec, 20, 1./lambdaForTv)




        x_rec_old = copy.deepcopy(x_rec)
        x_rec = im3ddenoise(y_rec, 20, 1./lambdaForTv)
        t_old = t
        t = (1 + np.sqrt(1 + 4 * t ** 2)) / 2
        y_rec = x_rec + (t_old - 1) / t * (x_rec - x_rec_old)
        avgtoc = time.clock()
        avgtime.append(abs(avgtic-avgtoc))
        relativeError.append(np.linalg.norm((y_rec-trueimg).ravel(),2)/np.linalg.norm(trueimg.ravel(),2))
    print('Average time taken for each iteration for FISTA:' + str(sum(avgtime)/len(avgtime)) + '(s)')
    return y_rec ,relativeError

def ista(y, geo, angles, niter, trueimg):
    lmbda = 0.1
    y_rec = np.zeros((geo.nVoxel), dtype=np.float32)
    x_rec = copy.deepcopy(y_rec)
    t = 1
    #L = 2.e7
    L = 2.e8
    bm = 1 / L
    relativeError = []
    avgtime = []
    for i in range(niter):
        avgtic = time.clock()
        y_rec += bm * 2 * (tigre.Atb((y - tigre.Ax(y_rec, geo, angles, 'interpolated')),
                                     geo, angles, 'matched'))
        if i ==0:
            print('ISTA algorithm in progress.')
            toc = time.clock()
        if i ==1:
            tic = time.clock()
            print('Estimated time until completion: ' +str(niter*(tic-toc))+ '(s)')

        lambdaForTv = 2 * bm * lmbda

        y_rec = im3ddenoise(y_rec, 20, 1./lambdaForTv)
        avgtoc = time.clock()
        avgtime.append(abs(avgtic-avgtoc))
        relativeError.append(np.linalg.norm((y_rec - trueimg).ravel(), 2) / np.linalg.norm(trueimg.ravel(), 2))
    print('Average time taken for each iteration for ISTA:' + str(sum(avgtime)/len(avgtime)) + '(s)')
    return y_rec ,relativeError

def find_nan(img,line,iteration):
    for i in img.ravel():
        if str(i) == 'nan':
            print(img)
            print('NAN DETECTED at line ' +str(line))
            print('iteration ' +str(iteration))

            raise SystemExit
angles = np.linspace(0,2*np.pi,200)


geo = tigre.geometry(mode='cone', nVoxel=np.array([256, 256, 256]), default_geo=True)

src_img = data_loader.load_head_phantom(geo.nVoxel)
proj = tigre.Ax(src_img,geo,angles)
numit = 400


from tigre.algorithms.conjugate_gradient_algorithms import CGLS

class CGLS_relativeerror(CGLS):
    def __init__(self,proj,geo,angles,niter,**kwargs):
        CGLS.__init__(self,proj,geo,angles,niter,**kwargs)

    def run_main_iter(self,trueimg):
        self.l2l = np.zeros([self.niter], dtype=np.float32)
        relativeError = []
        avgtime = []
        for i in range(self.niter):

            if i == 0:
                print("CGLS Algorithm in progress.")
                toc = time.clock()
            if i == 1:
                tic = time.clock()
                print('Esitmated time until completetion (s): ' + str((self.niter - 1) * (tic - toc)))
            avgtic = time.clock()
            q = tigre.Ax(self.__p__, self.geo, self.angles, 'ray-voxel')
            q_norm = np.linalg.norm(q.ravel(), 2)
            alpha = self.__gamma__ / (q_norm * q_norm)
            self.res += alpha * self.__p__
            avgtoc = time.clock()
            avgtime.append(abs(avgtic - avgtoc))
            for item in self.__dict__:
                if type(getattr(self, item)) == np.ndarray:
                    if np.isnan(getattr(self, item)).any():
                        raise ValueError('nan found for ' + item + ' at iteraton ' + str(i))

            aux = self.proj - tigre.Ax(self.res, self.geo, self.angles, 'ray-voxel')
            self.l2l[i] = np.linalg.norm(aux.ravel(), 2)
            if i > 0 and self.l2l[i] > self.l2l[i - 1]:
                print('re-initialization was called at iter:' + str(i))
                self.res -= alpha * self.__p__
                self.reinitialise_cgls()

            self.__r__ -= alpha * q
            s = tigre.Atb(self.__r__, self.geo, self.angles)
            s_norm = np.linalg.norm(s.ravel(), 2)

            gamma1 = s_norm * s_norm
            beta = gamma1 / self.__gamma__
            if self.log_parameters:
                self.parameter_history['alpha'][i] = alpha
                self.parameter_history['beta'][i] = beta
                self.parameter_history['gamma'][i] = self.__gamma__
                self.parameter_history['q_norm'][i] = q_norm
                self.parameter_history['s_norm'][i] = s_norm

            self.__gamma__ = gamma1
            self.__p__ = s + beta * self.__p__

            relativeError.append(np.linalg.norm((self.res - trueimg).ravel(), 2) / np.linalg.norm(trueimg.ravel(), 2))
        print('Average time taken for each iteration for CGLS:' + str(sum(avgtime) / len(avgtime)) + '(s)')
        return self.__r__, relativeError

from tigre.algorithms.art_family_algorithms import SIRT
class SIRT_relativeerror(SIRT):
    def __init__(self,proj,geo,angles,niter,**kwargs):
        SIRT.__init__(proj,geo,angles,niter,**kwargs)

    def run_main_iter(self,trueimg):
        """
        Goes through the main iteration for the given configuration.
        :return: None
        """
        Quameasopts = self.Quameasopts
        relativeError = []
        avgtime = []
        for i in range(self.niter):

            res_prev = None
            if Quameasopts is not None:
                res_prev = copy.deepcopy(self.res)
            if self.verbose:
                if i == 0:
                    print(str(self.name).upper() + ' ' + "algorithm in progress.")
                    toc = time.clock()
                if i == 1:
                    tic = time.clock()
                    print('Esitmated time until completetion (s): ' + str((self.niter - 1) * (tic - toc)))
            avgtic = time.clock()
            getattr(self, self.dataminimizing)()
            avgtoc = time.clock()
            avgtime.append(abs(avgtic-avgtoc))
            self.error_measurement(res_prev, i)
            relativeError.append(np.linalg.norm((self.res - trueimg).ravel(), 2) / np.linalg.norm(trueimg.ravel(), 2))
        print('Average time taken for each iteration for SIRT:' + str(sum(avgtime) / len(avgtime)) + '(s)')
        return self.res, relativeError

output_ista,relativeError_ista = ista(proj,geo,angles,numit,src_img)
output_fista,relativeError_fista = fista(proj,geo,angles,numit,src_img)
output_cgls,relativeError_cgls = CGLS_relativeerror(proj,geo,angles,numit).run_main_iter(src_img)
output_sirt,relativeError_sirt = SIRT_relativeerror(proj,geo,angles,numit).run_main_iter(src_img)

np.save(time.asctime()+'.npy', dict(output_fista= output_fista,
                                    relativeError_fista = relativeError_fista,
                                    output_ista = output_ista,
                                    relativeError_ista = relativeError_ista,
                                    output_cgls=output_cgls,
                                    relativeError_cgls=relativeError_cgls,
                                    output_sirt = output_sirt,
                                    relativeError_sirt = relativeError_sirt,
                                    ))

"""
from matplotlib import pyplot as plt
plt.plot(relativeError_fista,label='FISTA')
plt.plot(relativeError_cgls,label='CGLS')
plt.legend()
plt.show()
"""