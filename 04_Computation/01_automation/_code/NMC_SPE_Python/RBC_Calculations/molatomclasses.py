#Copyright Lieven Bekaert 2021,2022,2023 all rights reserved - lieven.bekaert@vub.be

class atom(object):
    def __init__(self,id,symbol,pos):
        self.id = id
        self.symbol = symbol
        self.pos = pos
        self.iselectrode=False

    def get_symbol(self):
        return self.symbol

    def __str__(self):
        return "ATOM id="+str(self.id)+" symbol="+self.symbol+" pos="+str(self.pos)+" iselectrode="+str(self.iselectrode)
