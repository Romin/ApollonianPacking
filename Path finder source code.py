from __future__ import division
from visual import*
from decimal import*

getcontext().prec = 100

##############################################################################

#Defining R, S and T

#R: doing permutation (13) or swapping which rational approx you're adding to

def R(array):

    temp = array[0]
    array[0] = array[2]
    array[2] = temp

    return array

#S: doing permutation (23) or doing the inversion in the rational approx

def S(array):

    temp = array[1]
    array[1] = array[2]
    array[2] = temp

    return array

#T: replacing the smallest circel that we care about with the biggest circle
#that is tangent to all of them or adding 1 to the rationa approx

def T(array):

    temp = array[2]
    array[2] = Decimal(2)*(array[0]+array[1]+array[2])-array[3]
    array[3] = temp

    return array

################################################################################

#S1*S2*S3

##epsilonOne = 1/8+(3/8)*(9+4*sqrt(5))
##epsilonTwo = 3/8 + (1/8)*(9+4*sqrt(5))
##epsilonThree = 1
##epsilonFour = 0

##epsilonOne = 1/8+(3/8)*(9-4*sqrt(5))
##epsilonTwo = 3/8 + (1/8)*(9-4*sqrt(5))
##epsilonThree = 1
##epsilonFour = 0

#S1*S2*S3*S2

##epsilonOne = 1/6+ (1/6)*(17+12*sqrt(2))
##epsilonTwo = 2
##epsilonThree = 1
##epsilonFour = 0

##epsilonOne = 1/6+ (1/6)*(17-12*sqrt(2))
##epsilonTwo = 2
##epsilonThree = 1
##epsilonFour = 0

#S1*S2*S3*S2*S1*S2

##epsilonOne = 6
##epsilonTwo = 119/10+(1/10)*(-49-20*sqrt(6))
##epsilonThree = 1
##epsilonFour = 0

##epsilonOne = 6
##epsilonTwo = 119/10+(1/10)*(-49+20*sqrt(6))
##epsilonThree = 1
##epsilonFour = 0

#S1*S2*S3*S4

##epsilonOne = -0.938326
##epsilonTwo = -0.324674
##epsilonThree = -0.112342
##epsilonFour = -0.0388719

##epsilonOne = 0.938326
##epsilonTwo = 0.324674
##epsilonThree = 0.112342
##epsilonFour = 0.0388719

#S1*S2*S3*S4*S2

epsilonOne = Decimal(3)/Decimal(56) + (Decimal(9)/Decimal(56))*(Decimal(73)+Decimal(12)*sqrt(Decimal(37)))
epsilonTwo = Decimal(75)/Decimal(28) + (Decimal(1)/Decimal(28))*(Decimal(73)+Decimal(12)*sqrt(Decimal(37)))
epsilonThree = Decimal(19)/Decimal(56) + (Decimal(1)/Decimal(56))*(Decimal(73)+Decimal(12)*sqrt(Decimal(37)))
epsilonFour = Decimal(1)

##epsilonOne = 3/56 + (9/56)*(73-12*sqrt(37))
##epsilonTwo = 75/28 + (1/28)*(73-12*sqrt(37))
##epsilonThree = 19/56 + (1/56)*(73-12*sqrt(37))
##epsilonFour = 1

#S1*S2*S3*S2*S4*S1*S3

##epsilonOne = 47/16 + (1/16)*(385+16*sqrt(579))
##epsilonTwo = 5/16 + (1/48)*(385+16*sqrt(579)) 
##epsilonThree = 25/3
##epsilonFour = 1

##epsilonOne = 47/16 + (1/16)*(385-16*sqrt(579))
##epsilonTwo = 5/16 + (1/48)*(385-16*sqrt(579))
##epsilonThree = 25/3
##epsilonFour = 1

#S1*S2*S3*S4*S2*S3*S2*S3

##epsilonOne = -0.942596
##epsilonTwo = -0.31012
##epsilonThree = -0.11464
##epsilonFour = -0.0468571

#S1*S2*S1*S2*S3*S1*S2*S3*S4

##epsilonOne = Decimal(str(0.8515888153473002510311442077037374683179444773007433976295426175604623792447276687147746355987063198))
##epsilonTwo = Decimal(str(0.5227364663847090603124620567210819088996352301076624941534953230779830127427502039560025923298769774))
##epsilonThree = Decimal(str(0.03927725725442969897925177366279809302956181242405485191923250847235195065304534242379270323551061834))
##epsilonFour = Decimal(str(0.0006110233866560384788306861220806083225303555444873294744092888799351358215838237558176031635476571400))

#S1*S2*S1*S2*S3*S1*S2*S3*S4

##epsilonOne = Decimal(81)/Decimal(1736) + (Decimal(523)/Decimal(1736))*(Decimal(2313) + Decimal(68)*sqrt(Decimal(1157)))
##epsilonTwo = Decimal(219)/Decimal(1736) + (Decimal(321)/Decimal(1736))*(Decimal(2313) + Decimal(68)*sqrt(Decimal(1157)))
##epsilonThree = Decimal(71)/Decimal(217) + (Decimal(3)/Decimal(217))*(Decimal(2313) + Decimal(68)*sqrt(Decimal(1157)))
##epsilonFour = Decimal(1)

epsilons = [Decimal(str(epsilonOne)), Decimal(str(epsilonTwo)), Decimal(str(epsilonThree)), Decimal(str(epsilonFour))]
path = []

print(epsilons)

unsortedEpsilons = epsilons + []

epsilons.sort()

for i in range(len(epsilons)):
    for j in range(len(epsilons)):
        if epsilons[i] == unsortedEpsilons[j]:
            print("The sort sends", j+1, "to", i+1)
        

##print(epsilons)

##norm = epsilons[2]
##
####for i in range(len(epsilons)):
####    epsilons[i] = epsilons[i]/norm
##
##epsilons.sort()

##print"The sorted, normalized epsilons are: \n", epsilons

while 1:

    actAgain = raw_input("Press enter to move again")

    if epsilons[2] <  epsilons[0]:

        R(epsilons)
        S(epsilons)
        path.append("R")

    elif epsilons[2] < epsilons[1]:

        S(epsilons)
        path.append("S")

    elif epsilons[2] >= epsilons[1]:

        T(epsilons)
        path.append("T")

    else:
        print("Huh something else is the case")

##    print(epsilons)

    print( "".join(path))
    
