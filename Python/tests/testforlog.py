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


def testandlog(alglist, geo, angles, niter, saveresult=False, newsubdirectory=False, **kwargs):
    nangles = angles.shape[0]
    # createlogfile

    dirname = os.path.dirname(__file__)
    kwargs.update(dict(blocksize=20))
    # create top level folder if this does not exist

    if 'logs' not in os.listdir(dirname):
        os.system('mkdir ' + 'logs')
    logdirectory = os.path.join(dirname, 'logs')
    logdate = time.strftime("%a_%b_%d")
    if logdate not in os.listdir(logdirectory):
        os.system('mkdir ' + os.path.join(logdirectory, logdate))
    subdirectory = os.path.join(logdirectory, logdate)
    if 'subdirectoryname' in kwargs:
        subdirectoryname = kwargs['subdirectoryname']
    else:
        subdirectoryname = 'test' + str(len(os.listdir(subdirectory)) - 1)

    if newsubdirectory:
        subdirectoryname = 'test' + str(len(os.listdir(subdirectory)))

    # create subdirectory for test if this does not exist

    if subdirectoryname not in os.listdir(subdirectory):
        os.system('mkdir ' + os.path.join(subdirectory, subdirectoryname))
    subdirectory = os.path.join(logdirectory, logdate, subdirectoryname)
    timestamp = str(time.asctime()).replace(' ', '_')

    # create/append to log file

    logfilename = None
    for filename in os.listdir(subdirectory):
        if os.path.splitext(filename)[-1] == '.log':
            logfilename = os.path.split(filename)[-1]
    logflist = []
    if logfilename is not None:
        logflist.extend(open(os.path.join(subdirectory, logfilename), 'r').readlines())
    else:
        logfilename = timestamp + '.log'
    try:
        algsuccess = np.load(os.path.join(subdirectory, os.path.splitext(logfilename)[0]) + '.npy').item()
    except Exception:
        algsuccess = dict()
    # TODO: need to fix directory path here.
    repo_path = os.path.join(dirname, '..', '..', '.git/')
    repo = Repo(repo_path)
    commit = list(repo.iter_commits(repo.active_branch))[0]
    uname = str(os.uname())
    current_config = dict(uname=uname, commithash=commit.hexsha)
    try:
        prev_config = open(os.path.join(subdirectory, os.path.splitext(logfilename)[0] + '_config') + '.npy').item()
    except Exception:
        prev_config = dict()
        prev_config.update(current_config)

    for item in current_config:
        if prev_config[item] != current_config[item]:
            logflist.append('------------------------------------------------\n')
            logflist.append('Configuration changed for' + item + ':\n')
            logflist.append('From: ' + str(prev_config[item]))
            logflist.append('To  : ' + str(current_config[item]))

    if len(logflist) == 0:  # If file is empty
        # MACHINE

        logflist.append(uname + '\n')

        # GIT header for logfile

        if not repo.bare:
            logflist.append('------------------------------------------------\n')
            logflist.append(make_repository_str(repo))

            # create list of commits then print some of them to stdout

            logflist.append(make_commit_str(commit))
            logflist.append('\n------------------------------------------------\n')
        else:
            print('Could not load repository at {} :('.format(repo_path))

        logflist.append('GEOMETRY used for instance of testandlog: \n')
        for item in geo.__dict__:
            logflist.append(item + ': ' + str(getattr(geo, item)) + '\n')
        logflist.append('------------------------------------------------\n')

    # Run the algorithms

    source_img = data_loader.load_head_phantom(number_of_voxels=geo.nVoxel)
    proj = tigre.Ax(source_img, geo, angles)

    for alg in alglist:
        logflist.append(str(alg).upper() + ' ' + time.asctime() + '\n')
        logflist.append('nproj: ' + str(nangles) + ' niter: ' + str(niter) + '\n')
        if np.sum(angles[:, 1]) != 0 or np.sum(angles[:, 2]) != 0:
            spherical = True
        else:
            spherical = False
        logflist.append('spherical projection: ' + str(spherical) + '\n')
        algpassed = False
        try:
            tic = time.clock()
            res = getattr(algs, alg)(proj, geo, angles, niter=niter, **kwargs)
            algpassed = True
        except Exception:
            formatedlines = traceback.format_exc()
            logflist.append('ERROR at ' + str(time.strftime("%H:%M:%S") + '\n'))
            logflist.append(''.join(formatedlines) + '\n')
        finally:
            toc = time.clock()
            algsuccess.update({alg: algpassed})
            np.save(os.path.join(subdirectory, os.path.splitext(logfilename)[0]), algsuccess)
            np.save(os.path.join(subdirectory,os.path.splitext(logfilename)[0])+'_config',current_config)
            pass
        if algpassed:
            logflist.append('total time taken: ' + str(abs(tic - toc)) + '\n')
        logflist.append('------------------------------------------------\n')
        logf = open(os.path.join(subdirectory, logfilename), 'w')
        logf.write(''.join(logflist))
        logf.close()
        if saveresult and algpassed:
            warnings.filterwarnings("ignore")
            plt.figure()
            plt.subplot(3, 1, 1)
            plt.imshow(res[geo.nVoxel[0] / 2])
            plt.title('results for ' + alg)
            plt.ylabel('dim 0')

            plt.subplot(3, 1, 2)
            plt.imshow(res[:, geo.nVoxel[1] / 2])
            plt.ylabel('dim 1')

            plt.subplot(3, 1, 3)
            plt.imshow(res[:, :, geo.nVoxel[2] / 2])
            plt.ylabel('dim 2')
            plt.savefig(os.path.join(subdirectory, alg + '_' + geo.mode + '_' + timestamp + '.png'))
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
    commit_str = ('\n'.join(l))
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


def log_summary(logdate, tobool=False):
    dirname = os.path.dirname(__file__)
    if 'logs' not in os.listdir(dirname):
        raise Exception('Logs not found under ' + dirname)
    logdirectory = os.path.join(dirname, 'logs')

    if logdate not in os.listdir(logdirectory):
        raise Exception('Date directory ' + logdate + ' not found.')
    logdatedirectory = os.path.join(logdirectory, logdate)
    logflist = []
    logfilename = ''
    for filename in os.listdir(logdirectory):
        if os.path.splitext(filename)[-1] == '.log':
            logfilename = os.path.split(filename)[-1]

    if logdate in logfilename:
        logflist.extend(open(os.path.join(logdirectory, logfilename), 'r').readlines())
    else:
        logfilename = logdate + '_testsummary.log'
        logflist.append(str(os.uname()) + '\n')
        logflist.append('------------------------------------------------\n')
        repo_path = os.path.join(dirname, '..', '..', '.git/')
        repo = Repo(repo_path)
        if not repo.bare:
            logflist.append('------------------------------------------------\n')
            logflist.append(make_repository_str(repo))

            # create list of commits then print some of them to stdout
            commit = list(repo.iter_commits(repo.active_branch))[0]
            logflist.append(make_commit_str(commit))
            logflist.append('\n------------------------------------------------\n')
        else:
            print('Could not load repository at {} :('.format(repo_path))

    toboolist = []
    logflist.append('summary carried out on ' + time.asctime() + '\n')
    for subdirectory in os.listdir(logdatedirectory):
        for testfile in os.listdir(os.path.join(logdatedirectory, subdirectory)):
            if os.path.splitext(testfile)[-1] == '.log':
                timestamp = os.path.splitext(testfile)[0]

                try:
                    testpassed = np.load(os.path.join(logdatedirectory, subdirectory, timestamp) + '.npy').item()
                    testconfig = np.load(os.path.join(logdatedirectory,subdirectory,timestamp)+'_config.npy').item()
                    formattedlist = []
                    for alg in testpassed:
                        toboolist.append(testpassed[alg])
                        if testpassed[alg]:
                            formattedlist.append('      ' + alg + ' : PASSED')
                        else:
                            formattedlist.append('      ' + alg + ' : ERROR')
                    logflist.append('\n------------------------------------------------\n')
                    logflist.append(subdirectory + ': \n' + str('\n'.join(formattedlist)) + '\n')
                    logflist.extend([item + ': ' + str(testconfig[item]) + '\n' for item in testconfig])
                except IOError:
                    logflist.append(subdirectory + ' no .npy file found. \n')

    logflist.append('\n All tests passed:' + str(all(toboolist)) + '\n')
    logflist.append('\n------------------------------------------------\n')

    with open(os.path.join(logdirectory, logfilename), 'w') as logf:
        logf.write(''.join(logflist))

    if tobool:
        return all(toboolist)


if __name__ == '__main__':
    method_name = sys.argv[1]
    param_name = sys.argv[2]
    tobool = False
    if len(sys.argv) == 4:
        tobool = sys.argv[3]

    getattr(sys.modules[__name__], method_name)(param_name, tobool)
