# Proportional multiple resonant controller design

A repository of algorithms for computing robust proportional multiple resonant (PMR) controller and current feedback gain for UPS systems. Each algorithm has its own folder, containing the main file and its function file.

The following algorithms are available:
- [PMR robust design](./PMR_robust_design)
- [PMR design](./PMR_design)
- [DT PMR design](./DT_PMR_design)

##

For a robust PMR desgin, the following output voltage response is obtained for linear load:
<img width="583" height="455" alt="ex_load" src="https://github.com/user-attachments/assets/98883c8c-bd2f-4382-aa9b-d29993de32da" />

The following output voltage response is obtained for non-linear load:
<img width="583" height="455" alt="ex_nl_load" src="https://github.com/user-attachments/assets/7535fdc3-b74c-4b72-a3eb-ed64aa999908" />

Note how the steady-state output voltage perfectly follows the reference for linear load and can reject disturbances caused by the nonlinear load, depending on the number of modes considered.