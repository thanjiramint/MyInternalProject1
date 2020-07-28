import math
import numpy as np
import matplotlib.pyplot as plt
import types



def plot(f,start,stop, space=100):
    ### Plot ##
    x = np.linspace(start,stop,space)
    y = f(x)
    # plot the function
    plt.plot(x,y)
    plt.grid()
    # show the plot
    plt.show()

def swap (a,b):
    temp=a
    a=b
    b=temp
    return a,b

def goldensectionsearch(f, x1, x2, tolerance=10 ** -5, plot1 = False):

    
    
    '''
    Golden Section Search algorithm.    
    Implements the golden section search, which is an algorithm that used for
    finding a minimum of a unimodal continuous function f(x) on the interval [x1,x2].
    The idea of this method is recursively narrowing the bracketing interval until
    it reaches the smallest value or specified tolerance.
    
    Parameters
    ----------
    function : The function f should be unimodal on the interval [x1,x2] 
        which recommended declare by the lambda function.
        Eg : f = lambda x: (x+3)*(x-1)**2

    x1 : A numeric specifying the lower interval lies on the function.
    x2 : A numeric specifying the upper interval lies on the function.
    tolerance : A float specifying the required accuracy which the default value is assigned to 10^-5
    plot1 : A boolean
        Perform graphical plot for visualize the minimum output of the algorithm.
        If True, the graphical plot will be appear.
        If False, the graphical plot will be not appear. Default is False.
    
    '''
    
    #Check inputs interval : Must be integer
    if (isinstance(x1, int) == False or isinstance(x2, int) == False):
        return print("Interval input must be integer")
    
    #Swap if lower bound is more than upper bound
    if (x1 > x2):
        x1,x2  = swap(x1,x2)
    
    ##Check input Function
    if (isinstance(f, types.FunctionType) == False):
        return print("Input function is incorrect.")
             
    ## Declar for plotting
    lower1 = x1
    upper2 =x2
    
    phi = (math.sqrt(5)-1)/2  ## 0.618
    xl = x1 + phi*(x2-x1)
    xr = x2 - phi*(x2-x1)
    f1 = f(xl)
    f2 = f(xr)

    while math.fabs(x2 - x1) > tolerance :
        if(f2>f1):
            x1=xr
            xr=xl
            f2=f1
            xl= x1+phi*(x2-x1)
            f1 = f(xl)
        else:
            x2=xl
            xl=xr
            f1=f2
            xr= x2-phi*(x2-x1)
            f2 = f(xr)
    if(plot1 == True):
        plot(f,lower1,upper2)
        
    #return "Golden Section Search gives minimum value equals to {}".format((lower+upper)/2)
    return (x1+x2)/2

def bisection(f,x1, x2,tolerance=10 ** -5, plot1 = False, print_i= False):

    
    '''
    Bisection algorithm.    
    Implements Bisection search is known as numerical method, which is an algorithm that used for finding
    the root of a nonlinear equation f(x) = 0. This method is one of bracketing method, 
    which requires the two guesses. The root can be defined onthe interval with specific 
    lower and upper bound [x1,x2] where f(x1) and f(x2)have opposite signs.
    
    Parameters
    ----------
    function : The function f should be unimodal on the interval [x1,x2] 
        which recommended declare by the lambda function.
        Eg : f = lambda x: (x+3)*(x-1)**2

    lower : An integer specifying the lower interval lies on the function.
    upper : An integer specifying the upper interval lies on the function.
    tolerance : A float specifying the required accuracy which the default value is assigned to 10^-5
    plot1 : bool 
        Perform graphical plot for visualize the minimum output of the algorithm.
        If True, the graphical plot will be appear.
        If False, the graphical plot will be not appear. Default is False.
    ploti : bool 
        Perform printing the iteration number, whcih it the total steps used in one process 
        If True, prints the iteration number.
        If False, not prints the iteration number. Default is False.
    
    '''
    xl = x1
    xr = x2
    i=0
    #check sign f(a), f(b) must have difference sign
    #means two values a and b are chosen for which f(xr) > 0 and f(xl) < 0
    if f(x1)*f(x2) >0 :
        mid = 0
        print("The input intervals is not meet the root")
        
    else:
        
        #the interval magnitude is less than tolerance, means the value is close enough to be the root of the function
        while (abs(xl - xr) >= tolerance):
            i=i+1
            
            #find the mid point
            mid = (xl+xr)/2.0
            
            #if mid =0 , means we found the root of the function, which is mid point
            if f(mid) == 0:
                break;
                
            #if f(mid) has same sign as f(xl), root is on between  mid-xr, new xl assigns to mid point.
            if np.sign(f(xl)) == np.sign(f(mid)):
                xl = mid
             #if f(mid) has same sign as f(xr), root is on between  xl-mid, new xr assigns to mid point.
            else:
                xr = mid

    if(print_i == True):
        print("iteration = {}" .format(i))
    if(plot1 == True):
        plot(f,x1,x2)
        #print("Root is equals to {}" .format(mid))

       
    #return "Bisection Method gives root at x = {}".format(mid)
    return mid

## Xn+1 = Xn – f(Xn)/f'(Xn)
def newton_raphson(f,fprime, x, tolerance = 10 ** -5 ,iteration =1000, plot1 = False):
    '''
    Newton-Rapshon method.    
    Implements the Newton-Raphson method , which is a way to quickly find a good approximation
    for the root of a real-valued function f(x) = 0. It uses the idea that 
    a continuous and differentiable function can be approximated by a straight line tangent to it.
    
    Parameters
    ----------
    function : The function f should be unimodal on the interval [a,b] 
        which recommended declare by the lambda function.
        Eg : f = lambda x: (x+3)*(x-1)**2
    function : The derivative function of initial function.
    
    x : An integer specifying random intial number lies on the function
    tolerance : A float specifying the required accuracy which the default value is assigned to 10^-5
    iterlation : An integer specifying the maximum step required on running one process.
    plot1 : bool 
        Perform graphical plot for visualize the minimum output of the algorithm.
        If True, the graphical plot will be appear.
        If False, the graphical plot will be not appear. Default is False.
   '''
    

    for i in range(iteration):
        h=f(x)/fprime(x);
        xnew= x-h ;
        #print(" At Iteration no. {}, x = {} \n".format(i, xnew))
        
        if(abs(h) < tolerance):
            break;    
        x=xnew
        
    if(plot1 == True):
        plot(f,-x,x)
    
    #return "After {} iterations, root = {}".format(i, x)
    return x

def secant(f,x1, x2,tolerance=10 ** -5,iteration =1000, plot1 = False):
    '''
    Secant Method.    
    Implements the secant method, which is an algorithm used to approximate the roots of
    a given function f. The method is based on approximating f using secant lines.
    The secant method algorithm requires the selection of two initial approximations 
    x1 and x2, which may or may not bracket the desired root, but which are 
    chosen reasonably close to the exact root.
    
    Parameters
    ----------
    function : The function f should be unimodal on the interval [a,b] 
        which recommended declare by the lambda function.
        Eg : f = lambda x: (x+3)*(x-1)**2

    x1 : An integer specifying the lower interval lies on the function.
    x2 : An integer specifying the upper interval lies on the function.
    tolerance : A float specifying the required accuracy which the default value is assigned to 10^-5
    iterlation : An integer specifying the maximum step required on running one process.
    plot1 : bool 
        Perform graphical plot for visualize the minimum output of the algorithm.
        If True, the graphical plot will be appear.
        If False, the graphical plot will be not appear. Default is False.

    
    '''
    upper = x1
    lower = x2
    #check sign f(a), f(b) must have difference sign
    #means two values a and b are chosen for which f(xr) > 0 and f(xl) < 0
    if f(x1)*f(x2) >0 :
        xnew = 0
        print("The input intervals is not meets the root")
    
    else:
        for i in range(iteration):
            xnew = x2 - f(x2)*(x2-x1)/(f(x2)-f(x1))
            #xnew = (x1*f(x2) - x2*f(x1))/ (f(x2) - f(x1))
            if abs(xnew-x2) < tolerance:
                break
            else:
                x1 = x2
                x2 = xnew
    
    if(plot1 == True):
        plot(f,upper,lower)
    
    #return "Secant Method gives root at x = {} , at iterations {}".format(xnew, i)
    return xnew




    
def brent(f,x1, x2,tolerance=10 ** -11):
      
    '''
    Brent's method
    Implement the brent method, which is a root finding hybrid algorithm
    which combines the bracketing methods and the open method 
    where the open methods are inverse quadratic interpolation, 
    secant method and bracketing method is bisection method.
    
    Parameters
    ----------
    function : The function f should be unimodal on the interval [a,b] 
        which recommended declare by the lambda function.
        Eg : f = lambda x: (x+3)*(x-1)**2

    x1 : An integer specifying the lower interval lies on the function.
    x2 : An integer specifying the upper interval lies on the function.
    tolerance : A float specifying the required accuracy which the default value is assigned to 10^-5
    
    '''
    
    
    
    #check sign of functions
    if (f(x1)*f(x2) > 0):
        return print("error : Root is not found on this function.")
    
    if (abs(f(x1)) < abs(f(x2))):
        x1,x2  = swap(x1,x2)
        
    c=x1
    mflag = 1
    d=0
    i=0
    s=0
    while( f(x2) == 0 or abs(x1-x2) >= tolerance ):
        if (f(x1) != f(x2) and f(x2) != f(c) and (f(x1)-f(x2)) !=0 and (f(x1)-f(c)) !=0  and (f(x2)-f(c)) != 0 ):
            #Inverse Quadratic Interpolation
            s = ( x1*f(x2)*f(c)/( (f(x1)-f(x2))*(f(x1)-f(c)) ) ) + ( x2*f(x1)*f(c)/( (f(x2)-f(x1))*(f(x2)-f(c)) ) ) + ( c*f(x1)*f(x2)/( (f(c)-f(x1))*(f(c)-f(x2)) ) )    
        
        else :
            #secant
            s = x2 - (x2-x1)/(f(x2)-f(x1)) * f(x2)
            
        #Condition 1&2
        if ( not(s > (3*x1+x2)/4 and s < x2) or (mflag ==1 and abs(s-x2) >= abs(x2-c)/2 )):
            #bisection
            s = (x1+x2)/2
            mflag = 1
        #Condition 3&4&5
        elif ( i != 0  and ( (mflag ==0 and abs(s-x2) >= abs(c-d)/2) or (mflag ==1 and abs(x2-c)< tolerance ) or (mflag ==0 and abs(c-d)/2) < tolerance) ):
            #bisection        
            s = (x1+x2)/2
            mflag = 1
            
        else:
            mflag = 0
        
        d=c
        c=x2
        if (f(x1)*f(s) < 0 ):
            x2 = s
        else:
            x1 = s
            
        if(abs(f(x1)) < abs(f(x2))):
            swap(x1,x2)
        i+=1
 
    return s

