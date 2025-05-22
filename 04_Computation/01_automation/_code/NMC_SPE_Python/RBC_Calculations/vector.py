#Copyright Lieven Bekaert 2021,2022,2023 all rights reserved - lieven.bekaert@vub.be

#Vector file
import math
class vector:
    def __init__(self, tupleorlistinput):
        self.x=tupleorlistinput[0]
        self.y=tupleorlistinput[1]
        self.z=tupleorlistinput[2]

    def getlength(self):
        return math.sqrt((self.x**2)+(self.y**2)+(self.z**2))

    def __str__(self):
        return "vector x="+str(self.x)+" y="+str(self.y)+" z="+str(self.z)+" length="+str(self.getlength())

    def __mul__(self,other):
        return vector((self.x*other,self.y*other,self.z*other))

    def __div__(self,other):
        return vector((self.x/other,self.y/other,self.z/other))

    def __add__(self,other):
        return vector((self.x+other.x,self.y+other.y,self.z+other.z))

    def __sub__(self,other):
        return vector((self.x-other.x,self.y-other.y,self.z-other.z))

    def totuple(self):
        return (self.x,self.y,self.z)

    def finddotproduct(self,other):
        return self.x*other.x + self.y*other.y + self.z*other.z

    def __pow__(self,other):
        return vector((self.enablezerosupport_for_power(self.x,other),self.enablezerosupport_for_power(self.y,other),self.enablezerosupport_for_power(self.z,other)))

    def enablezerosupport_for_power(self,toverify,other):
        if toverify<=0:
            return 0
        else:
            return toverify**other

    def absolute(self):
        self.x=abs(self.x)
        self.y=abs(self.y)
        self.z=abs(self.z)
        return self

    def __le__(self,other):
        #IS <=
        return (self.x<=other.x)and(self.y<=other.y)and(self.z<=other.z)

    def __ge__(self,other):
        #IS >=
        return (self.x>=other.x)and(self.y>=other.y)and(self.z>=other.z)
