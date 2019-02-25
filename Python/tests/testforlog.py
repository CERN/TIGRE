from __future__ import print_function
import tigre
import sys, traceback
import os
import numpy as np
from tigre.demos.Test_data import data_loader
from matplotlib import pyplot as plt
import time
import tigre.algorithms as algs
import warnings
from matplotlib import pyplot as plt
warnings.filterwarnings("error")
from git import Repo
def testandlog(alglist, geo, angles, niter, logfilename=None, saveresult=False, newsubdirectory=True,**kwargs):
    nangles = angles.shape[0]
    # createlogfile

    dirname = os.path.dirname(__file__)

    # create top level folder if this does not exist

    if 'logs' not in os.listdir(dirname):
        os.system('mkdir ' + 'logs')
    logdirectory = os.path.join(dirname, 'logs')
    logdate = time.strftime("%a_%d_%Y")
    if logdate not in os.listdir(logdirectory):
        os.system('mkdir ' + os.path.join(logdirectory,logdate))
    subdirectory = os.path.join(logdirectory,logdate)
    subdirectoryname = 'test' + str(len(os.listdir(subdirectory))-1)

    if newsubdirectory:
        subdirectoryname = 'test' +str(len(os.listdir(subdirectory)))

    os.system('mkdir ' + os.path.join(subdirectory, subdirectoryname))
    subdirectory = os.path.join(logdirectory,logdate, subdirectoryname)
    timestamp = str(time.asctime()).replace(' ', '_')
    if logfilename is not None:
        logf = open(os.path.join(subdirectory, logfilename) + '.log', 'w')
    else:
        logf = open(os.path.join(subdirectory, str(geo.mode) + '_' + timestamp + '.log'), 'w')

    # GIT header for logfile


    repo_path=os.path.join(dirname,'..','..','.git/')
    repo = Repo(repo_path)
    if not repo.bare:
        print('Repo at {} successfully loaded.'.format(repo_path))
        logf.write('------------------------------------------------\n')
        logf.write(make_repository_str(repo))

        # create list of commits then print some of them to stdout
        commit = list(repo.iter_commits(repo.active_branch))[0]
        logf.write(make_commit_str(commit))
        logf.write('\n------------------------------------------------\n')
    else:
        print('Could not load repository at {} :('.format(repo_path))


    source_img = data_loader.load_head_phantom(number_of_voxels=geo.nVoxel)

    proj = tigre.Ax(source_img, geo, angles)

    logf.write('GEOMETRY used for instance of testandlog: \n')
    for item in geo.__dict__:
        logf.write(item + ': ' + str(getattr(geo, item)) + '\n')
    logf.write('------------------------------------------------\n')
    for alg in alglist:
        logf.write(str(alg).upper() + ' ' + time.asctime() + '\n')
        logf.write('nproj: ' + str(nangles) + ' niter: ' + str(niter) + '\n')
        if np.sum(angles[:,1])!=0 or np.sum(angles[:,2])!=0:
            spherical = True
        else:
            spherical = False
        logf.write('spherical projection: ' + str(spherical) + '\n')
        algpassed = False
        try:
            tic = time.clock()
            res = getattr(algs, alg)(proj, geo, angles, niter=niter,**kwargs)
            algpassed = True
        except Exception:
            formatedlines = traceback.format_exc()
            logf.write('ERROR at ' + str(time.strftime("%H:%M:%S")+'\n'))
            logf.write(''.join(formatedlines) + '\n')
        finally:
            toc = time.clock()
            pass
        if algpassed:
            logf.write('total time taken: ' + str(abs(tic - toc)) + '\n')
        logf.write('------------------------------------------------\n')

        if saveresult and algpassed:
            warnings.filterwarnings("ignore")
            plt.figure()
            plt.subplot(3, 1, 1)
            plt.imshow(res[geo.nVoxel[0]/2])
            plt.title('results for '+alg)
            plt.ylabel('dim 0')

            plt.subplot(3, 1, 2)
            plt.imshow(res[:,geo.nVoxel[1]/2])
            plt.ylabel('dim 1')

            plt.subplot(3, 1, 3)
            plt.imshow(res[:,:,geo.nVoxel[2]/2])
            plt.ylabel('dim 2')
            plt.savefig(os.path.join(subdirectory, alg+ '_'+geo.mode+'_' +timestamp+'.png'))
            warnings.filterwarnings("error")
def make_commit_str(commit):
    l = []
    l.append(str(commit.hexsha))
    l.append("\"{}\" by {} ({})".format(commit.summary,
                                     commit.author.name,
                                     commit.author.email))
    l.append(str(commit.authored_datetime))
    l.append(str("count: {} and size: {}".format(commit.count(),
                                              commit.size)))
    commit_str=('\n'.join(l))
    return commit_str
def make_repository_str(repo):
    l = []
    l.append('Repo description: {}'.format(repo.description))
    l.append('Repo active branch is {}'.format(repo.active_branch))
    for remote in repo.remotes:
        l.append('Remote named "{}" with URL "{}"'.format(remote, remote.url))
    l.append('Last commit for repo is {}.'.format(str(repo.head.commit.hexsha)))
    repo_str = ('\n'.join(l))
    return repo_str