import numpy as np

def g1(x,y):
    return (2*x+x**2-y)/2

def g2(x,y):
    return (2*x-x**2+8)/9 + (4*y-y**2)/4

def fixedpoint_method_multivariable(g1, g2, x0, y0, n, tol):
    x = g1(x0, y0)
    y = g2(x0, y0)

    for i in range(n):
      print (f"Iteration {i}: x = {x:.3f}, y = {y:.3f}")
      print (f'residual: {abs(x - x0):.3e}, {abs(y - y0):.3e}')
      if abs(x - x0) < tol and abs(y - y0) < tol:
          print (f"Converged after {i} iterations.")
          print('We have found our root.')
          return x, y
      x0 = x
      y0 = y
      x = g1(x0, y0)
      y = g2(x0, y0)
    
    print("We have not found our root.")
    return x, y

zero = fixedpoint_method_multivariable(g1, g2, -1, 1, 100, 1e-9)
print (zero)
    