from __future__ import division
from visual.controls import*
import cmath

################################################################################

#Example of an unbounded packing. Copy and paste this into the input

##0.851588815347300251031144207703737468317944477300743397629542617560462379244727668714774635598706,0.522736466384709060312462056721081908899635230107662494153495323077983012742750203956002592329876,0.039277257254429698979251773662798093029561812424054851919232508472351950653045342423792703235510,0.000611023386656038478830686122080608322530355544487329474409288879935135821583823755817603163547

#########################################################################################################################

#Changable parameters

keepOldCircles = True


#########################################################################################################################

#Useful functions

def areClose(numberOne, numberTwo):
    
    if abs(numberOne-numberTwo) < 0.1:
        return True
    else:
        return False

def complexMultiply(vectorOne, vectorTwo):

    temp = vector(vectorOne[0]*vectorTwo[0] - vectorOne[1]*vectorTwo[1], vectorOne[0]*vectorTwo[1] + vectorOne[1]*vectorTwo[0],0)
    return temp

def lawOfCosines(sOne, sTwo, sThree):

    '''returns the cos(theta) angle between sOne and sTwo'''

    temp = (sOne**2 + sTwo**2 - sThree**2)/(2*sOne*sTwo)
##    print(temp)
    return temp

def thirdCircleCenter(cOne, rOne, rTwo, rThree):

    '''takes the center of the first circle and the radii of the three, and
       determines the position of the center of the third circle.'''

    cosTheta = lawOfCosines(rOne + rTwo, rOne + rThree, rTwo + rThree)
    base = (rOne+rThree)*cosTheta
    height = (rOne+rThree)* sin(acos(cosTheta)) #(sqrt(1-((cosTheta)**2)))

##    print(base, height)

    position = cOne + vector(base, height, 0)

##    print(position)
    return position

eps = 0.0001
def appEq(u, v): # approximately equal excluding signs
    return abs(abs(u) - abs(v)) <= eps

# test if 2 circles are approximately equal
def appEqCir(a, b):
    return (abs(a[0] - b[0]) <= eps and appEq(a[1], b[1]))

# test if 3 circles are tangent to each other
def isTangent(a, b, c):
    flag = True
    if appEqCir(a, b) or appEqCir(a, c) or appEqCir(b, c): # all circles must be different
        flag = False
    if not isTanCir(a, b):
        flag = False
    if not isTanCir(a, c):
        flag = False
    if not isTanCir(b, c):
        flag = False
    return flag

# test if 2 circles are tangent
def isTanCir(a, b):
    return (appEq(a[0] - b[0], abs(a[1]) + abs(b[1])) or appEq(a[0] - b[0], abs(a[1]) - abs(b[1])))

# test if 2 circles are intersecting
def isIntCir(a, b):
    dist = abs(a[0] - b[0]) # distance between 2 centers
    rmin = min(abs(a[1]), abs(b[1])) # min radius
    rmax = max(abs(a[1]), abs(b[1])) # max radius
    return (dist + eps < rmin + rmax and dist + rmin > rmax + eps)


# test if the 4th circle is tangent to previous 3
def isTan4(a, b, c, d):
    return isTangent(a, b, d) and isTangent(a, c, d) and isTangent(b, c, d)

# input: 3 circles w/ each as a tuple: (vector(x,y,0),r)
# output: 4 tangent circles as a tuple
def cif(a, b, c):
    # curvatures of the circles
    k1 = 1.0 / a[1]
    k2 = 1.0 / b[1]
    k3 = 1.0 / c[1]
    # curvatures of tangent circles
    temp0 = 2.0 * math.sqrt(k1 * k2 + k1 * k3 + k2 * k3)
    temp1 = k1 + k2 + k3
    k40 = temp1 + temp0
    k41 = temp1 - temp0
    # centers of tangent circles
##    temp0 = 2.0 * math.sqrt(k1*k2*complexMultiply(a[0],b[0])+k1*k3*complexMultiply(a[0],c[0])+k2*k3*complexMultiply(b[0],c[0]))
    temp0 = 2.0 * cmath.sqrt(k1*k2*a[0]*b[0]+k1*k3*a[0]*c[0]+k2*k3*b[0]*c[0])
    temp1 = k1 * a[0] + k2 * b[0] + k3 * c[0]
    c40 = (temp1 + temp0) / k40
    c41 = (temp1 - temp0) / k41
    c42 = (temp1 - temp0) / k40
    c43 = (temp1 + temp0) / k41
    # radius of tangent circles
    r0 = 1.0 / k40
    r1 = 1.0 / k41
    # there are 4 solutions and only 2 are correct (tangent) circles
    # those are c0 and c1 or c2 and c3
    c0 = (c40, r0)
    c1 = (c41, r1)
    c2 = (c42, r0)
    c3 = (c43, r1)
    if isTan4(a, b, c, c0) and isTan4(a, b, c, c1):
##        print(c0, c1)
        return (c0, c1)
    else:
##        print(c2, c3)
        return (c2, c3)

########################################################################################################################

#Getting the starting data

startingCircles = (input('What are the curvatures?'))
oldCircles = []

curvatures = list(startingCircles)
curvatures.sort()
for i in range(len(curvatures)):
    curvatures[i] = 1/curvatures[i]

radii = []
radii.extend(curvatures)

##radii.sort()
##print("the radii are", radii)

for i in range(len(curvatures)):
    curvatures[i] = 1/curvatures[i]

##print("the curcatures are", curvatures)
##print("the radii are", radii)

#########################################################################################################################



radiusOne = radii[0]
radiusTwo = radii[1]
radiusThree = radii[2]
radiusFour = radii[3]

circleThickness =(radiusOne)/100

centerOne = vector(0,0,0)
centerTwo = vector(radiusOne+radiusTwo,0,0)
centerThree = thirdCircleCenter(centerOne,radiusOne,radiusTwo,radiusThree)

##line = curve(pos=[centerOne, centerTwo, centerThree, centerOne])

circleOne = ring(pos=centerOne, radius = radiusOne, axis = (0,0,1), thickness = circleThickness)
circleTwo = ring(pos=centerTwo, radius = radiusTwo, axis = (0,0,1), thickness = circleThickness)
circleThree = ring(pos=centerThree, radius = radiusThree, axis = (0,0,1), thickness = circleThickness)

##print(abs(centerOne - centerTwo))
##print(abs(centerOne - centerThree))
##print(abs(centerTwo - centerThree))

desQuads = cif((complex(centerOne[0],centerOne[1]), radiusOne), (complex(centerTwo[0],centerTwo[1]), radiusTwo), (complex(centerThree[0],centerThree[1]), radiusThree))

if radiusFour == desQuads[0][1]:
    print("Awesome! That's a valid quadruple")
    circleFour = ring(pos=vector(real(desQuads[0][0]), imag(desQuads[0][0]),0), radius = desQuads[0][1], axis = (0,0,1), thickness = circleThickness, color=color.white)
    centerFour = vector(real(desQuads[0][0]), imag(desQuads[0][0]),0)
    radiusFour = desQuads[0][1]
    
elif radiusFour == desQuads[1][1]:
    print("Awesome! That's a valid quadruple")
    circleFour = ring(pos=vector(real(desQuads[1][0]), imag(desQuads[1][0]),0), radius = desQuads[1][1], axis = (0,0,1), thickness = circleThickness, color=color.white)
    centerFour = vector(real(desQuads[1][0]), imag(desQuads[1][0]),0)
    radiusFour = desQuads[1][1]
    
else:
    print("oops, that 4th curvature isn't a valid quadruple. The valid quadruples are", 1.0/desQuads[0][1], "or", 1.0/desQuads[1][1])
    circleFour = ring(pos=vector(real(desQuads[0][0]), imag(desQuads[0][0]),0), radius = desQuads[0][1], axis = (0,0,1), thickness = circleThickness, color=color.red)
    circleFour2 = ring(pos=vector(real(desQuads[1][0]), imag(desQuads[1][0]),0), radius = desQuads[1][1], axis = (0,0,1), thickness = circleThickness, color=color.red)
    centerFour = vector(real(desQuads[0][0]), imag(desQuads[0][0]),0)
    radiusFour = radius = desQuads[0][1]



#########################################################################################################################################################################

#########################################################################################################################################################################


#Button interactions

def S1():

    global circleOne, centerOne, radiusOne, circleTwo, centerTwo, radiusTwo, circleThree, centerThree, radiusThree, circleFour, centerFour, radiusFour
    
    curvatures[0] = 2*(curvatures[1]+curvatures[2]+curvatures[3])-curvatures[0]
##    print("changed curvatures")
    print(curvatures[0])
    
    newDesQuads = cif((complex(centerTwo[0],centerTwo[1]), radiusTwo), (complex(centerThree[0],centerThree[1]), radiusThree),(complex(centerFour[0],centerFour[1]), radiusFour))

    print(1/newDesQuads[0][1], 1/newDesQuads[1][1])

    if areClose(curvatures[0],1/newDesQuads[0][1]):
        tempCircle = circleOne
        tempCircle.visible = keepOldCircles
        oldCircles.append(tempCircle)
        circleOne = ring(pos=vector(real(newDesQuads[0][0]), imag(newDesQuads[0][0]),0), radius = newDesQuads[0][1], axis = (0,0,1), thickness = circleThickness, color=color.red)
        centerOne = vector(real(newDesQuads[0][0]), imag(newDesQuads[0][0]),0)
        radiusOne = newDesQuads[0][1]
        print("New Circle")
    elif areClose(curvatures[0],1/newDesQuads[1][1]):
        tempCircle = circleOne
        tempCircle.visible = keepOldCircles
        oldCircles.append(tempCircle)
        circleOne = ring(pos=vector(real(newDesQuads[1][0]), imag(newDesQuads[1][0]),0), radius = newDesQuads[1][1], axis = (0,0,1), thickness = circleThickness, color=color.red)
        centerOne = vector(real(newDesQuads[1][0]), imag(newDesQuads[1][0]),0)
        radiusOne = newDesQuads[1][1]
        print("New Circle")

def S2():

    global circleOne, centerOne, radiusOne, circleTwo, centerTwo, radiusTwo, circleThree, centerThree, radiusThree, circleFour, centerFour, radiusFour
    
    curvatures[1] = 2*(curvatures[0]+curvatures[2]+curvatures[3])-curvatures[1]
##    print("changed curvatures")
    print(curvatures[1])
    
    newDesQuads = cif((complex(centerOne[0],centerOne[1]), radiusOne), (complex(centerThree[0],centerThree[1]), radiusThree),(complex(centerFour[0],centerFour[1]), radiusFour))

    print(1/newDesQuads[0][1], 1/newDesQuads[1][1])

    if areClose(curvatures[1],1/newDesQuads[1][1]):
        tempCircle = circleTwo
        tempCircle.visible = keepOldCircles
        oldCircles.append(tempCircle)
        circleTwo = ring(pos=vector(real(newDesQuads[1][0]), imag(newDesQuads[1][0]),0), radius = newDesQuads[1][1], axis = (0,0,1), thickness = circleThickness, color=color.blue)
        centerTwo = vector(real(newDesQuads[1][0]), imag(newDesQuads[1][0]),0)
        radiusTwo = newDesQuads[1][1]
        print("New Circle")

    elif areClose(curvatures[1],1/newDesQuads[0][1]):
        tempCircle = circleTwo
        tempCircle.visible = keepOldCircles
        oldCircles.append(tempCircle)
        circleTwo = ring(pos=vector(real(newDesQuads[0][0]), imag(newDesQuads[0][0]),0), radius = newDesQuads[0][1], axis = (0,0,1), thickness = circleThickness, color=color.blue)
        centerTwo = vector(real(newDesQuads[0][0]), imag(newDesQuads[0][0]),0)
        radiusTwo = newDesQuads[0][1]
        print("New Circle")

def S3():

    global circleOne, centerOne, radiusOne, circleTwo, centerTwo, radiusTwo, circleThree, centerThree, radiusThree, circleFour, centerFour, radiusFour
    
    curvatures[2] = 2*(curvatures[0]+curvatures[1]+curvatures[3])-curvatures[2]
##    print("changed curvatures")
    print(curvatures[2])

    newDesQuads = cif((complex(centerOne[0],centerOne[1]), radiusOne), (complex(centerTwo[0],centerTwo[1]), radiusTwo),(complex(centerFour[0],centerFour[1]), radiusFour))

    print(1/newDesQuads[0][1], 1/newDesQuads[1][1])

    if areClose(curvatures[2],1/newDesQuads[1][1]):
        tempCircle = circleThree
        tempCircle.visible = keepOldCircles
        oldCircles.append(tempCircle)
        circleThree = ring(pos=vector(real(newDesQuads[1][0]), imag(newDesQuads[1][0]),0), radius = newDesQuads[1][1], axis = (0,0,1), thickness = circleThickness, color=color.green)
        centerThree = vector(real(newDesQuads[1][0]), imag(newDesQuads[1][0]),0)
        radiusThree = newDesQuads[1][1]
        print("New Circle")

    elif areClose(curvatures[2],1/newDesQuads[0][1]):
        tempCircle = circleThree
        tempCircle.visible = keepOldCircles
        oldCircles.append(tempCircle)
        circleThree= ring(pos=vector(real(newDesQuads[0][0]), imag(newDesQuads[0][0]),0), radius = newDesQuads[0][1], axis = (0,0,1), thickness = circleThickness, color=color.green)
        centerThree = vector(real(newDesQuads[0][0]), imag(newDesQuads[0][0]),0)
        radiusThree = newDesQuads[0][1]
        print("New Circle")

def S4():

    global circleOne, centerOne, radiusOne, circleTwo, centerTwo, radiusTwo, circleThree, centerThree, radiusThree, circleFour, centerFour, radiusFour
    
    curvatures[3] = 2*(curvatures[0]+curvatures[1]+curvatures[2])-curvatures[3]
##    print("changed curvatures")
    print(curvatures[3])
    
    newDesQuads = cif((complex(centerOne[0],centerOne[1]), radiusOne), (complex(centerTwo[0],centerTwo[1]), radiusTwo),(complex(centerThree[0],centerThree[1]), radiusThree))

    print(1/newDesQuads[0][1], 1/newDesQuads[1][1])

    if areClose(curvatures[3],1/newDesQuads[1][1]):
        tempCircle = circleFour
        tempCircle.visible = keepOldCircles
        oldCircles.append(tempCircle)
        circleFour = ring(pos=vector(real(newDesQuads[1][0]), imag(newDesQuads[1][0]),0), radius = newDesQuads[1][1], axis = (0,0,1), thickness = circleThickness, color=color.orange)
        centerFour = vector(real(newDesQuads[1][0]), imag(newDesQuads[1][0]),0)
        radiusFour = newDesQuads[1][1]
        print("New Circle")

    elif areClose(curvatures[3],1/newDesQuads[0][1]):
        tempCircle = circleFour
        tempCircle.visible = keepOldCircles
        oldCircles.append(tempCircle)
        circleFour= ring(pos=vector(real(newDesQuads[0][0]), imag(newDesQuads[0][0]),0), radius = newDesQuads[0][1], axis = (0,0,1), thickness = circleThickness, color=color.orange)
        centerFour = vector(real(newDesQuads[0][0]), imag(newDesQuads[0][0]),0)
        radiusFour = newDesQuads[0][1]
        print("New Circle")

def SGrow():

    par = curvatures[0] + curvatures[1] + curvatures[2] + curvatures[3]

    testOne = sum([2*(curvatures[1]+curvatures[2]+curvatures[3])-curvatures[0], curvatures[1], curvatures[2], curvatures[3]])
    testTwo = sum([curvatures[0], 2*(curvatures[0]+curvatures[2]+curvatures[3])-curvatures[1], curvatures[2], curvatures[3]])
    testThree= sum([curvatures[0], curvatures[1], 2*(curvatures[0]+curvatures[1]+curvatures[3])-curvatures[2], curvatures[3]])
    testFour = sum([curvatures[0], curvatures[1], curvatures[2], 2*(curvatures[0]+curvatures[1]+curvatures[2])-curvatures[3]])

    if testOne < par:
        print("Do S1")
        S1()
        return
    if testTwo < par:
        print("Do S2")
        S2()
        return
    if testThree < par:
        print("Do S3")
        S3()
        return
    if testFour < par:
        print("Do S4")
        S4()
        return
    else:
        print('You have the maximal configuration here')
        return

##########################################################################################################################################################################################3

#Adding the buttons

controlBox = controls(x=500)

S1Button = button(text='S1', action=lambda: S1(), pos=(-60,0))
S2Button = button(text='S2', action=lambda: S2(), pos=(-20,0))
S3Button = button(text='S3', action=lambda: S3(), pos=(20,0))
S4Button = button(text='S4', action=lambda: S4(), pos=(60,0))
SGrowButton = button(text='Grow', action=lambda: SGrow(), pos=(0,-50,0), width = 100)

while 1:
    controlBox.interact()

    if scene.mouse.clicked:
        m = scene.mouse.getclick()
        scene.center = m.pos

##    circleThickness = 0.01 ##1/(radiusOne*1000)





        

    
