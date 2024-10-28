import numpy as np

def f(x):
  return np.sin(x) + np.arctan(x)


def anderson_acceleration(f, x0, tol, k_max, m):
  # Initialize iterates and residuals
  x = [x0, f(x0)]
  g = [f(x[0]) - x[0], f(x[1]) - x[1]]

  G_k = [g[1] - g[0]] # Matrix of increments in residuals
  X_k = [x[1] - x[0]] # Matrix of increments in x

  k = 1
  
  while k < k_max and abs(g[-1]) > tol:
    m_k = min(k, m)

    if k > 1:
      G_k_matrix = np.array(G_k[-m_k:]).reshape(-1,1) # Reshape to column vector
      print (G_k_matrix)
      Q, R = np.linalg.qr(G_k_matrix)
      print (Q)
      print (R)
      reshaped = np.array(g[-m_k:]).reshape(-1,1)
      print (reshaped)
      print (np.array(g[-m_k:]).reshape(2,1))
      gamma_k = np.linalg.solve(R, Q.T @ np.array(g[-m_k:]))
    else:
      gamma_k = 0

    # Compute new iterate and new residual
    if k > 1:
      x_new = x[-1] + g[-1] - (np.array(X_k[-m_k:]) + np.array(G_k[-m_k:])) @ gamma_k
    else:
      x_new = x[-1] + g[-1]

    g_new = f(x_new) - x_new

    # Append new iterate and new residual
    x.append(x_new)
    g.append(g_new)

    # Update increment matricies with new elements
    if k > 1:
      X_k.append(x[-1] - x[-2])
      G_k.append(g[-1] - g[-2])

    # Keep only the last m_k elements
    if len(X_k) > m_k:
      X_k = X_k[-m_k:]
      G_k = G_k[-m_k:]
    
    k += 1

  print(f"Computed fixed point {x[-1]:.6f} after {k} iterations")
  return x[-1] # Return the fixed point

x0 = 1
N_max = 100
tol = 1e-6
m = 3

fixed_point = anderson_acceleration(f, x0, tol, N_max, m)