import numpy as np

# Define the function f(x)
def f(x):
   return (10/x -4*x)**0.5

def anderson_acceleration(f, x0, tol_res, k_max, m):
   # Initialize the vectors x and g
   x1 = f(x0)
   x = np.array([x0, x1])
   g = np.array([f(x0) - x0, f(x1) - x1])

   # Initialize matrices for increments in residuals and x values
   G_k = np.array([[g[1] - g[0]]])  # Start with the first residual
   X_k = np.array([[x[1] - x[0]]])  # Start with the first difference in x

   k = 1
   while k < k_max and abs(g[k]) > tol_res:
      if (k > 3):
         m_k = 3
      else:
         m_k = 1
      #m_k = min(k, m)

      # Perform QR decomposition of G_k
      Q, R = np.linalg.qr(G_k)
      print (f'Q.T:{Q.T}')
      #print (f'R:{R}')
      #print (f'G_k:{G_k}')
      #print (f'shape:{np.shape(G_k)}')
      print (f'g:{g}')
      print (f'g[k]:{g[k]}')
      #print (f'g_mk:{np.array(g[:,-m_k:])}')
      gamma_k, residuals, rank, s = np.linalg.lstsq(R, Q.T @ np.array(g[k]).reshape(-1, 1), rcond=None) 

      # Compute new iterate and new residual
      new_x = x[k] + g[k] - (X_k + G_k) @ np.array([gamma_k]).T
      new_g = f(new_x[0]) - new_x[0]

      # Update x and g vectors with new elements
      x = np.append(x, new_x[0][0][0])
      g = np.append(g, new_g[0][0])

      # Update increment matrices with new elements
      X_k = np.append(X_k, new_x - x[k]).reshape(1,-1)
      G_k = np.append(G_k, new_g - g[k]).reshape(1,-1)


      
      # Keep only the last m_k elements]
      if np.size(X_k) > m_k:
         #print("in if statement")
         X_k = np.array(X_k[0][-m_k:]).reshape(1,-1)
         G_k = np.array(G_k[0][-m_k:]).reshape(1,-1)
         

      
      #print (f"Iteration {k}: x = {x[k]:.10f}, f(x) = {f(x[k]):.10f}, g(x) = {abs(f(x[k]) - x[k]):.6e}")
      #print(f"X_k: {X_k}")
      k += 1

   print(f"Computed fixed point {x[-1]:.6f} after {k} iterations")
   return x[-1]

x0 = 1.5
N_max = 100
tol = 1e-9
m = 3

fixed_point = anderson_acceleration(f, x0, tol, N_max, m)
print(f(fixed_point))
