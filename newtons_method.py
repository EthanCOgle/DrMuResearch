def f (x):
  return x**3 + 4*x**2 - 10

def df (x):
   return 3*x**2 + 8*x 

def newtons_method(f, p0, n, tol):
    p = p0 - f(p0)/df(p0)

    for i in range(n):
      if abs(p - p0) < tol:
          print('We have found our root in ' + str(i) + ' iterations.')
          return p 
      p0 = p
      p = p0 - f(p0)/df(p0)
    
    print("We have not found our root.")
    return p

zero = newtons_method(f, 1.5, 100, 0.00001)
print (zero)
print (f(zero))