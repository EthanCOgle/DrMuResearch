# Numerical Optimization Research

By Ethan Ogle with mentorship from Dr. Lin Mu at the University of Georgia.

## Anderson's Acceleration Method

Anderson's Acceleration is similar to fixed-point method, but instead of just looking at the current iteration to make the next guess, it has a parameter m = depth that looks at the past m guesses to make the next one.

## Advantages to Anderson's Acceleration

Anderson's Acceleration is faster and more reliable than fixed-point method. It typically has a faster convergence rate of super-linear and converges to the solution when fixed-point method would diverge instead. 

Compared to Newton's method, Anderson's Acceleration, at some points, can parallel the quadratic convergence rate of Newton's method while still converging when Newton's method would diverge. Additionally, Newton's method is reliant on a good initial guess while Anderson's Method is more flexible.

## What's Next?

Currently, we are researching ways to have a more adaptable depth that is variable from iteration to iteration. For example, starting out with a smaller depth and then increasing with the iterations.
