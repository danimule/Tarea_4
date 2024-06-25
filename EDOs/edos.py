# EDOs/edos.py

"""
Contiene multiples formas de resolver ecuaciones diferenciales oridnarias.

Este modulo permite al usuario resolver sus ODEs mediante los metodos de Euler, RK2 y RK4.

"""


def Euler(f,x0,t0,tf,N):
    """
    Metodo de Euler para resolver ODE.

    Examples:
        >>> import numpy as np
        >>> def f(x,t):
        >>>     return (-1*x)**3 + np.sin(t)
        >>> x,t = Euler(f,0,0,10,20)

    Args:
        f (function): Ecuación a resolver.
        x0 (int): Valor inicial para la variable dependiente.
        t0 (int): Valor inicial para la variable independiente.
        tf (int): Valor final para la variable independiente.
        N (int): Número de particiones del rango de la variable independiente.

    Returns:
        x (vector): resultados de la variable dependiente.
        t (vector): resultados de la variable independiente.
    """
    t = np.linspace(t0,tf,N)
    x = np.zeros(t.size)
    x[0] = x0
    h = t[1] - t[0]
    
    for i in range(t.size - 1):
        x[i+1] = x[i] + h*f(x[i],t[i])
    
    return x,t


def RK2(f,x0,t0,tf,N):

    """
    Método de Runge-Kutta de Segundo Orden para resolver ODE.

    Examples:
        >>> import numpy as np
        >>> def f(x,t):
        >>>     return -x**3 + np.sin(t)
        >>> x,t = RK2(f,0,0,10,20)

    Args:
        f (function): Ecuación a resolver.
        x0 (int): Valor inicial de la variable dependiente.
        t0 (int): Valor inicial de la variable independiente.
        tf (int): valor final de la variable independiente.
        N (int): Número de particiones del rango de la variable independiente.

    returns:
        x (vector): resultados de la variable dependiente.
        t (vector): resultados de la variable independiente.

    """
    
    t = np.linspace(t0,tf,N)
    x = np.zeros(t.size)
    x[0] = x0
    h = t[1] - t[0]
    
    for i in range(t.size-1):
        k1 = h*f(x[i],t[i])
        k2 = h*f(x[i] + (k1/2), t[i] + (h/2))
        n = x[i] + (h/2)*f(x[i],t[i])
        x[i+1] = x[i] + h*f(n, t[i] + h/2)
    
    return x,t



def RK4(f,x0,t0,tf,N):
    """
    Método de Runge-Kutta de Cuarto Orden para resolver ODE.

    Examples:
        >>> import numpy as np
        >>> def f(x,t):
        >>>     return -x**3 + np.sin(t)
        >>> x,t = RK4(f,0,0,10,20)

    Args:
        f (func): Ecuación a resolver.
        x0 (int): Valor inicial para la variable dependiente.
        t0 (int): Valor inicial para la variable independiente.
        tf (int): Valor final para la variable independiente.
        N (int): Número de particiones del rango de la variable independiente.


    returns:
        x (vector): resultados de la variable dependiente.
        t (vector): resultados de la variable independiente.

     """
    t = np.linspace(t0,tf,N)
    x = np.zeros(t.size)
    x[0] = x0
    h = t[1] - t[0]
    
    for i in range(t.size - 1):
        k1 = h*f(x[i],t[i])
        k2 = h*f(x[i] + (k1/2), t[i] + (h/2))
        k3 = h*f(x[i] + (k2/2), t[i] + (h/2))
        k4 = h*f(x[i] + k3, t[i] + h)
        x[i+1] = x[i] + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
    
    return x,t