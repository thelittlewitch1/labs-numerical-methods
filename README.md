# labs-numerical-methods
Realisation of some numerical methods with scilab, c++, java, python.

---

Shot description of reposytory content:

1. Solving systems of linear equations.
    
    A. Gaus method:
      https://en.wikipedia.org/wiki/Gaussian_elimination
    
      relise in: 
    
      /c++/matrix/, function Gaus
    
      /java/matrix/, class gauss_method
    
      /scilab/gaus.sce
    
    B. Gaus Seydel method:
      https://en.wikipedia.org/wiki/Gauss–Seidel_method
    
      relise in: 
    
      /c++/matrix/, function GaussSeidel
    
      /java/matrix/, class gauss_seidel_method
    
      /scilab/gaus-zeydel.sce
    
    C. Jacobi method:
      https://en.wikipedia.org/wiki/Jacobi_method
    
      relise in: 
    
      /c++/matrix/, function Jacobi
    
      /java/matrix/, class jacobi_method
    
      /scilab/jacobi.sce
      
2. Solving systems of nonlinear equations: 
  https://en.wikipedia.org/wiki/Newton%27s_method

    relise in:

    /c++/Newton_method
  
    /python/Newton
  
    /scilab/Newton
  
3. Finding max and min eigen values of matrix by power method
https://en.wikipedia.org/wiki/Power_iteration

    relise in:

    /c++/Eigenvalue
  
    /python/Power_method
  
    /scilab/Eigenvalue
  
4. Solving systems of ordinary differential equations
 
    A. Forward euler method:
      https://en.wikipedia.org/wiki/Euler_method
    
      relise in: 
    
      /python/SoDE.py, function forward_euler_method
    
    B. Implicit euler method:
      https://en.wikipedia.org/wiki/Backward_Euler_method
    
      relise in: 
    
      /python/SoDE.py, function implicit_euler_method
    
    C. shikhman method:
    
      relise in: 
    
      /python/SoDE.py, function shikhman_method
      
5. Solving a partial differential equation (will be added)
     
    A. Forward euler method:
      description
    
      relise in: 
    
      path
    
    B. Implicit euler method:
      description 
    
      relise in: 
    
      path
