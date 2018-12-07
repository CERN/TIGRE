from _minTV import minTV
class Regularisation(object):

    def minimizeTV(self,res_prev,dtvg):
        return minTV(res_prev,dtvg,self.numiter_tv)


