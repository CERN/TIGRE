
import numpy as np
import matplotlib as mt
pgf_with_rc_fonts = {"pgf.texsystem": "xelatex"}
mt.rcParams.update(pgf_with_rc_fonts)
from mt import pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
results = np.load('resultsforconvergence.npy').item()
xx =np.linspace(1,100,100)
plt.semilogy(xx,results['relativeError_fista'],label=r"\textit{Time/Iteration}=21.6")
plt.ylabel(r'\textbf{Iteration}')
plt.legend()
plt.show()
