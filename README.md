# Robust proportional multiple resonant controller

A code to compute robust proportional multiple resonant (PMR) controller and current feedback gain in order to ensure UPS stability and performance. Considering the state-space AVERAGED model, PMR controller and current gains are tuned by means of regional pole placement with LMI regions. Output to linear and non-linear load are plotted without the need for graphical environment circuits, which makes it simple for learning and practical design applications. Although the inverter's switching effects are not considered, it allows a brief simulation to evaluate the performance of controller obtained.

Only the following specifications are required:
- Number of resonant structures in the controller (h)
- Damping factor of the resonant modes (xi)
- The desired real part of the smallest pole of the closed-loop system (s)

Consider an example with h = 3, xi = 0, s = 100. The following output voltage response is obtained for linear load:
<img width="583" height="455" alt="ex_load" src="https://github.com/user-attachments/assets/98883c8c-bd2f-4382-aa9b-d29993de32da" />

The following output voltage response is obtained for non-linear load:
<img width="583" height="455" alt="ex_nl_load" src="https://github.com/user-attachments/assets/7535fdc3-b74c-4b72-a3eb-ed64aa999908" />

Note how the steady-state output voltage perfectly follows the reference for linear load and can reject disturbances caused by the nonlinear load, depending on the number of modes considered. The transient response, however, is related with the smallest pole constraint specified in the design. Due to the robustness of the controller, similar performance is ensured over the entire load range, from no-load to nominal load conditions.
