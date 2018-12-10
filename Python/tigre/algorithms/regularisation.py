from _minTV import minTV
from _AwminTV import AwminTV
class Regularisation(object):

    def minimizeTV(self,res_prev,dtvg):
        return minTV(res_prev,dtvg,self.numiter_tv)

    def minimizeAwTV(self,res_prev,dtvg):
        return AwminTV(res_prev,dtvg,self.numiter_tv,self.delta)
