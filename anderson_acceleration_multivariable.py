import numpy as np

# Define the function f(x)
def g1(x,y):
    return (2*x+x**2-y)/2

def g2(x,y):
    return (2*x-x**2+8)/9 + (4*y-y**2)/4

def anderson_acceleration_multivariable(g1, g2, x0, y0, tol_res, k_max, m):
   # Initialize the vectors x and g
   x1 = g1(x0, y0)
   y1 = g2(x0, y0)
   x = np.array([[x0, x1], 
                 [y0, y1]]).T
   #print (np.shape(x))
   #print (x)
   g = np.array([[g1(x0, y0) - x0, g1(x1, y1) - x1], 
                 [g2(x0, y0) - y0, g2(x1, y1) - y1]]).T
   #print (np.shape(g))
   #print (f"g:{g}")
   #print ("G_K:")

   # Initialize matrices for increments in residuals and x values
   G_k = np.array([[g[0][1] - g[0][0]], 
                   [g[1][1] - g[1][0]]])  # Start with the first residual
   #print (G_k)
   #print ("X_K:") 
   
   X_k = np.array([[x[0][1] - x[0][0]], 
                   [x[1][1] - x[1][0]]])  # Start with the first difference in x
   #print (X_k)

   k = 1
   while k < k_max and (abs(g[0][k]) > tol_res or abs(g[1][k]) > tol_res):
      if (k > 3):
         m_k = 3
      else:
         m_k = 1
      #m_k = min(k, m)

      # Perform QR decomposition of G_k
      Q, R = np.linalg.qr(G_k, 'complete')
      #print (f'Q.T:{Q.T}')
      #print (f'R:{R}')
      #print (f'G_k:{G_k}')
      #print (f'shape:{np.shape(G_k)}')
      # print (f'g:{g}')
      # print (f'g_mk:{np.array(g[:,-m_k:])}')
      gamma_k, residuals, rank, s = np.linalg.lstsq(R, Q.T @ np.array(g[:,-m_k:]), rcond=None) 
      #print (f'gamma_k:{gamma_k}')

      # Compute new iterate and new residual
      #print (f'x[k]:{x[k]}')
      #print (f'g[k]:{g[k]}')
      #print (f'X_k:{X_k}')
      #print (f'G_k:{G_k}')
      #print (f'x[k] + g[k]:{np.add(x[k],g[k]).reshape(2,-1)}')
      #print (f'X_k + G_k:{np.add(X_k,G_k)}')
      #print (f'first:{np.add(np.add(x[k],g[k]).reshape(2,-1), -1*(np.add(X_k,G_k)))}')
      #print (f'second:{np.array([gamma_k]).T}')
      print (f"x_k + g_k:{np.add(x[:,-m_k:],g[:,-m_k:])}")
      print (f"X_k + G_k:{np.add(X_k,G_k)}")
      print (f"gamma_k:{np.array([gamma_k[0]]).T}")
      print (f"x_k + g_k:{np.shape(np.add(x[:,-m_k:],g[:,-m_k:]))}")
      print (f"X_k + G_k:{np.shape(np.add(X_k,G_k))}")
      print (f"gamma_k:{np.shape(np.array([gamma_k[0]]).T)}")
      new_x = np.add(np.add(x[:,-m_k:],g[:,-m_k:]), -1*(np.add(X_k,G_k))) @ np.array([gamma_k[0]]).T
      print (f'new_x:{new_x}')
      #print (f'new_x:{new_x[0][0]}')
      new_g= np.add(g1(new_x[0][0],new_x[0][1]), -1*new_x[0])
      #print (f'new_g:{new_g}')
      #new_g_y = f(new_x[0]) - new_x[0]

      # Update x and g vectors with new elements
      #print (f'x:{x}')
      #print (f'g:{g}')
      x = np.hstack((x, new_x[0]))
      g = np.hstack((g, new_g))
      #print (f'x:{x}')
      #print (f'g:{g}')

      # Update increment matrices with new elements
      #print (f'X_k:{X_k}')
      #print (f'G_k:{G_k}')
      #print (f'x[k]:{x[0][k]}')
      #print (f'new_x - x[k]:{np.add(new_x[0], -1*x[0][k])}')
      X_k = np.hstack((X_k, new_x[0] - x[0][k]))
      G_k = np.hstack((G_k, new_g - g[0][k]))
      #print (f'X_k:{X_k}')
      #print (f'G_k:{G_k}')


      k += 1
      # Keep only the last m_k elements
      # print (f'X_k:{X_k}')
      # print (f'G_k:{G_k}')
      # print (np.shape(X_k)[1])
      if np.shape(X_k)[1] > m_k:
         X_k = X_k[:,-m_k:]
         G_k = G_k[:,-m_k:]
      # print (f'X_k:{X_k}')
      # print (f'G_k:{G_k}')

      
      print (f"Iteration {k}")

   print(f"Computed fixed point {x[-1]:.6f} after {k} iterations")
   return x[-1]

x0 = -1
y0 = 1
N_max = 100
tol = 1e-9
m = 3

fixed_point = anderson_acceleration_multivariable(g1, g2, x0, y0, tol, N_max, m)
#print(f(fixed_point))
