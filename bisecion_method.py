def f(x):
    return x**3+4*x**2-10

def bisection_method(f, a, b, n, tol):
    if f(a) * f(b) > 0:
        print('Invalid a and b')
        return
    
    for i in range(n):
      p = (a + b)/2
      if f(p) == 0 or abs((b - a)/2) < tol:
          print('We have found our root.')
          return p 
      
      if f(p) * f(a) > 0:
          a = p
      else:
          b = p
    
    print("We have not found our root.")
    return p

zero = bisection_method(f, 1, 2, 100, 0.00001)
print (zero)
print (f(zero))
    