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
                 [y0, y1]])
   g = np.array([[g1(x0, y0) - x0, g1(x1, y1) - x1], 
                 [g2(x0, y0) - y0, g2(x1, y1) - y1]])

   # Initialize matrices for increments in residuals and x values
   G_k = np.array([[g[0][1] - g[0][0]], 
                   [g[1][1] - g[1][0]]])  # Start with the first residual
   
   X_k = np.array([[x[0][1] - x[0][0]], 
                   [x[1][1] - x[1][0]]])  # Start with the first difference in x
   
   k = 1
   while k < k_max and (abs(g[0][k]) > tol_res or abs(g[1][k]) > tol_res):
      m_k = min(k, m)

      # Perform QR decomposition of G_k
      Q, R = np.linalg.qr(G_k, 'complete')
      gamma_k, residuals, rank, s = np.linalg.lstsq(R, Q.T @ np.array(g[:,k]).reshape(2,1), rcond=None) 
      
      new_x = (np.add(np.add(x[:,k],g[:,k]).reshape(2,1), -1*(np.add(X_k,G_k) @ np.array(gamma_k)))).reshape(2,1)
      new_g= np.array([g1(new_x[0][0],new_x[1][0]) - new_x[0][0], g2(new_x[0][0],new_x[1][0]) -new_x[1][0]]).reshape(2,1)

      x = np.hstack((x, new_x))
      g = np.hstack((g, new_g))
      

      # Update increment matrices with new elements
      X_k = np.hstack((X_k, np.add(new_x, -1*x[:,k].reshape(2,1))))
      G_k = np.hstack((G_k, np.add(new_g, -1*g[:,k].reshape(2,1))))

      k += 1
      # Keep only the last m_k elements
      
      if np.shape(X_k)[1] > m_k:
         X_k = X_k[:,-(m_k+1):]
         G_k = G_k[:,-(m_k+1):]

      print (f"Iteration {k}")

   print(f"Computed fixed point {x[:,-1]} after {k} iterations")
   return x[:,-1]

x0 = -1
y0 = 1
N_max = 10
tol = 1e-9
m = 3

fixed_point = anderson_acceleration_multivariable(g1, g2, x0, y0, tol, N_max, m)
#print(f(fixed_point))
