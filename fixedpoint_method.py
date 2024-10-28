import numpy as np

def g(x):
    return (10/x -4*x)**0.5

def fixedpoint_method(f, p0, n, tol):
    p = g(p0)

    for i in range(n):
      print (f"Iteration {i}: x = {p:.6f}, f(x) = {f(p):.6f}, g(x) = {abs(p - f(p)):.6e}, ")
      if abs(p - p0) < tol:
          print (f"Converged after {i} iterations.")
          print('We have found our root.')
          return p 
      p0 = p
      p = g(p0)
    
    print("We have not found our root.")
    return p

zero = fixedpoint_method(g, 1.5, 100, 1e-9)
print (zero)
print (g(zero))
    