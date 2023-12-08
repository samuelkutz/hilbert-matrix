# hilbert-matrix
Various numerical analysis methods for solving a linear system of equations involving the Hilbert Matrix.

## The problem

Basically, we simply want to solve for $x$ the following system:

$$
Hx=b
$$

Where $H$ is the $n \times n$ Hilbert matrix given by the rule

$$
h_{ij} = \frac{1}{i + j - 1}
$$

And $b \in \mathbb{R}^{n}$ the vector given by 

$$
b = (1 \quad 1 \quad 1 \quad ... \quad 1)^{T}
$$

***In this implementation, we will fix the order of $n$
to be 100.***

The main metric used in this project is the **Relative residual norm** given by:

$$
    || r_{rel} || = \frac{|| b - Hx ||}{|| b ||} 
$$

## Methods

In this project, we implemented the following methods for aproximating the solution:

- LU decomposition
- Cholesky decomposition
- Jacobi Over-Relaxation (JOR)
- Successive Over-Relaxation (SOR)
- Gradient Ascent
- Conjugate Gradient

***All methods are implemented from scratch.***

## Results

<p style="font-style: italic; font-weight: lighter;">
    Don't mind the portuguese plots 
</p>

### [LU decomposition]()

Since LU is a relatively stable direct method for Symmetric Positive Definite (SDP) matrices, we can check the relative residual norm at each order $n$ compared to the desired tolerance

<p align="center">
    <img src="./plots/LU/residual_norm_LU.png" width="400">
</p>

### [Cholesky decomposition]()

Overall, this is the best method for SDP matrices. However, due to its bad conditioning, the floating-point errors acumulate so much it literally makes the $H$ decomposition lose its "positiveness".

In order to prevent this, we can create a new aproximated matrix. We just need to add small values $\lambda$ to de diagonal of $H$ and see what reduces the relative residual norm the most:

<p align="center">
    <img src="./plots/Cholesky/residual_norm_vs_lambda.png" width="400">
</p>
