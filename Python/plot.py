
import numpy as np
# import matplotlib as mt
# pgf_with_rc_fonts = {"pgf.texsystem": "xelatex"}
# mt.rcParams.update(pgf_with_rc_fonts)
from matplotlib import pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=14)
results = np.load('resultsforconvergence.npy').item()
xx =np.linspace(1,100,100)
plt.semilogy(xx,results['relativeError_fista'],label=r"FISTA")
plt.semilogy(xx,results['relativeError_cgls'],label=r"CGLS")
plt.semilogy(xx,results['relativeError_ista'],label=r"ISTA")
plt.semilogy(xx,results['relativeError_sirt'],label=r"SIRT")
plt.xlabel(r'\textit{Iterations}')
plt.ylabel(r'$||x^* - x_k||_2/||x^*||_2$')
plt.grid(True,which='minor',linestyle='--')
plt.grid(True)
plt.yticks([0.1,1],('0.1','1.0'),fontsize=10)
plt.xticks(fontsize=10)
plt.xlim(0,100)
plt.legend(fontsize=10)
plt.show()
