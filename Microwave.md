# Microwave Communications
# 10/11/23 Transmission lines and Telegrapher's Equation
- ## Striplines: velocity $v = \frac{c}{\sqrt{\epsilon_r}}$
- ## Microstrip lines: fields almost in board, while the wider the more field in board $v = \frac{c}{\epsilon_r}$
# Equivalent Dielectric Constant
- ## for finite-width microstrip lines:
$$\epsilon_e = \frac{\epsilon_r+1}{2}+\frac{\epsilon_r-1}{2}(1+12\frac{h}{W})^{-\frac{1}{2}}$$ 
# Phase velocity and Propagation Constant
$$v_{pd} = \frac{1}{\sqrt{\mu\epsilon}} = \frac{1}{\sqrt{\mu\epsilon_0}}\frac{1}{\sqrt{\epsilon_e}} = \frac{c}{\sqrt{\epsilon_e}}$$
$$v_{pa} = \frac{1}{\mu\epsilon_0} = 3\times10^8 m/s$$
$\beta$: It determines the sinusoidal amplitude / phase of the signal along a transmission line, at a constant time
$$k_0 = \omega\sqrt{\mu\epsilon_0} = \frac{\omega}{c} = \frac{2\pi{f}}{c} = \frac{2\pi}{\lambda}$$
$$\beta = \omega\sqrt{\mu\epsilon} = k_0\sqrt{\epsilon_e}$$
# Characteristic Impedance
$$\LARGE Z_0 = \left\{
    \begin{matrix}
    \frac{60}{\sqrt{\epsilon_e}}\ln({\frac{8h}{W}+\frac{W}{4h}}) & W/h \leq 1\\
    \frac{120\pi}{\sqrt{\epsilon_e}[W/h+1.393+0.667\ln(W/h+1.444)]} & W/h \geq 1
    \end{matrix}
    \right.$$
## for given $Z_0$ and $\epsilon_r$
$$\large \frac{W}{h} = \left\{
    \begin{matrix}
    \frac{8e^A}{e^{2A}-2} & W/h \leq 2 \\
    \frac{2}{\pi}[B-1-\ln(2B-1)+\frac{\epsilon_r-1}{2\epsilon_r}\{\ln(B-1)+0.39-\frac{0.61}{\epsilon_r}\}] & W/h \geq 2
    \end{matrix}
    \right. $$
$$ where\space A = \frac{Z_0}{60}\sqrt{\epsilon_r+1}{2}+\frac{\epsilon_r-1}{\epsilon_r+1}(0.23+\frac{0.11}{\epsilon_r})$$
$$ B = \frac{377\pi}{2Z_0\sqrt{\epsilon_r}}$$
# Transmission Line Theory
- Characteristic Impedance
$$Z_0 = \sqrt{\frac{R+j\omega{L}}{G+j\omega{C}}}$$
