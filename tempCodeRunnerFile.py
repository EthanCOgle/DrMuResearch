
      # print (f'G_k:{G_k}')
      # print (np.shape(X_k)[1])
      if np.shape(X_k)[1] > m_k:
         X_k = X_k[:,-m_k:]
         G_k = G_k[:,-m_k:]
      # print (f'X_k:{X_k}')
      # print (f'G_k:{G_k}')